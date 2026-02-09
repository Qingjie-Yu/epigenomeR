# Apply Partial Dependence and PDP Analysis
# Post: Train Random Forest and Linear Regression models on RNA-seq and HiPlex weighted gene coverage data,
#       then compute Partial Dependence Plots (PDP) and Individual Conditional Expectation (ICE) curves
#       for feature importance analysis.
# Supported models:
#   - "rf": Random Forest with configurable hyperparameters (n_estimators, max_depth, min_samples_split)
#   - "lr": Linear Regression with default parameters
#
# Parameters:
#   gene_select_dir_filename: Path to gene selection dictionary JSON file
#   gene_select_name: Name of gene selection list to use from dictionary. Default: "coding_all"
#   rnaseq_dir_filename: Path to RNA-seq data CSV file
#   model_design: Model design identifier for selecting appropriate cutoff threshold
#   wgc_dir_filename: Path to weighted gene coverage matrix (.feather file)
#   features_path: Path to JSON file containing feature names
#   rnaseq_cutoffs_path: Path to JSON file containing RNA-seq cutoff values per design
#   model_params_file: Path to JSON file containing model hyperparameters and best params
#   out_dir: Directory path for saving output files
#   random_seed: Random seed for reproducibility. Default: 42
#   frag_type: Fragment type for analysis. Default: "mixed"
#   n_cores: Number of cores for parallel processing. Default: NULL (uses 1 core)
#   remove_TV: Whether to remove features prefixed with "T_" or "V_". Default: FALSE
#   test_size: Proportion of data for test split. Default: 0.2
#
# Output: Saves to out_dir/<random_seed>/:
#         - model_performance_results.csv: RMSE and Pearson^2 for each model
#         - top_features.csv: Normalized feature importance from Random Forest
#         - top_feature_importance.json: Ranked feature names
#         - top_features_pdp.rds: PDP average, ICE individual curves, and grid values per feature
#         - rnaseq_wgc_all_X.csv: Processed feature matrix used for training
apply_par_and_pdp <- function(gene_select_dir_filename, rnaseq_dir_filename, wgc_dir_filename, features_path, rnaseq_cutoffs_path, model_params_file, out_dir, model_design = "rnaseq_vs_hiplex_rm_outlier_log", gene_select_name = "coding_all", random_seed = 42, frag_type = "mixed", n_cores = NULL, remove_TV = FALSE, test_size = 0.2) {
    # Load Libraries
    suppressPackageStartupMessages({
        library(arrow)
        library(caret)
        library(jsonlite)
        library(randomForest)
        library(ggplot2)
        library(MASS)  
        library(viridis)
        library(Metrics) 
        library(R6)
        library(ranger)
        library(parallel)
        library(doParallel)
        library(foreach)
        library(pdp)  
    }) 

    # Create Folder
    save_dir <- paste0(out_dir, "/", random_seed)
    if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
    }

    # Detect Number of Cores Used
    if (is.null(n_cores)) {
        n_cores <- 1
        
        max_available <- parallel::detectCores()
        if (interactive()) { 
            message(sprintf(
                "Using 1 core. Set n_cores argument (max %d available) for parallel processing.",
                max_available
            ))
        }
    } else {
        cat("Using", n_cores, "cores\n")
    }

    # Load Genes
    gene_select_dict <- fromJSON(gene_select_dir_filename) # list of category of genes
    # Question: only consider coding_all? then what the input expect to be? four category and combine to coding_all?
    gene_select_dict[["coding_all"]] <- c(gene_select_dict[["coding_cpg"]], gene_select_dict[["coding_non_cpg"]])
    print(gene_select_name) # input variable, keep?
    gene_select_list <- gene_select_dict[[gene_select_name]]

    # Load RNA Seq for Gene
    rnaseq <- read.csv(rnaseq_dir_filename, header = TRUE, row.names = 1)
    rnaseq$sqrt_V <- log10(rnaseq$V1V2 + 1)

    # Load Weighted Gene Coverage
    wgc_raw <- read_feather(wgc_dir_filename)

    # Name of Features
    features <- fromJSON(features_path)

    wgc_vals <- wgc_raw[, features]
    zero_hiplex_genes <- wgc_raw$pos[rowSums(wgc_vals) == 0]
    wgc_raw[, features] <- log10(wgc_raw[, features] + 1)
    wgc_raw[, features] <- apply(wgc_raw[, features], 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
    })

    features <- c("pos", features) # add one col name
    wgc_raw <- wgc_raw[, features, drop = FALSE]
    wgc_raw_list <- wgc_raw$pos
    overlap_gene_list <- intersect(wgc_raw_list, gene_select_dict[["coding_all"]]) # # input variable, keep?
    rnaseq_avail <- rnaseq[overlap_gene_list, c("gene_id", "sqrt_V"), drop = FALSE]

    # Load Cutoff Points
    rnaseq_cutoffs <- fromJSON(rnaseq_cutoffs_path)
    q99 <- rnaseq_cutoffs[[model_design]]
    rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, , drop = FALSE]
    zero_rnaseq_genes <- rnaseq_avail$gene_id[rnaseq_avail$sqrt_V == 0]
    zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)

    rnaseq_wgc_raw <- merge(rnaseq_avail, wgc_raw, by.x = "gene_id", by.y = "pos", all = FALSE)  # inner join
    rownames(rnaseq_wgc_raw) <- rnaseq_wgc_raw$gene_id
    rnaseq_wgc_raw <- rnaseq_wgc_raw[!rownames(rnaseq_wgc_raw) %in% zero_all_genes, , drop = FALSE]
    rnaseq_wgc <- rnaseq_wgc_raw[, !names(rnaseq_wgc_raw) %in% c("pos", "gene_id"), drop = FALSE]

    # Split Train And Test Data
    set.seed(random_seed)
    n <- nrow(rnaseq_wgc)
    test_indices <- sample(1:n, size = floor(n * test_size), replace = FALSE)
    train_indices <- setdiff(1:n, test_indices)

    rnaseq_wgc_train <- rnaseq_wgc[train_indices, , drop = FALSE]
    rnaseq_wgc_test <- rnaseq_wgc[test_indices, , drop = FALSE]

    rnaseq_wgc_train_X <- rnaseq_wgc_train[, !names(rnaseq_wgc_train) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_train_y <- rnaseq_wgc_train$sqrt_V
    rnaseq_wgc_test_X <- rnaseq_wgc_test[, !names(rnaseq_wgc_test) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_test_y <- rnaseq_wgc_test$sqrt_V
    rnaseq_wgc_all_X <- rnaseq_wgc[, !names(rnaseq_wgc) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_all_y <- rnaseq_wgc$sqrt_V

    rnaseq_wgc_all_X_out <- cbind(pos = rownames(rnaseq_wgc_all_X), rnaseq_wgc_all_X)
    write.csv(rnaseq_wgc_all_X_out, paste0(save_dir, "/rnaseq_wgc_all_X.csv"), row.names = FALSE)

    # Train Models
    all_models <- list()
    all_models[["coding_all"]] <- list()
    model_params <- fromJSON(model_params_file, simplifyDataFrame = FALSE)

    ## Define RandomForestRegressor class 
    RandomForestRegressor <- R6Class("RandomForestRegressor",
        public = list(
            random_state = NULL,
            params = list(),
            model = NULL,
            feature_names = NULL,
            
            initialize = function(random_state = NULL) {
            self$random_state <- random_state
            },
            
            set_params = function(...) {
            new_params <- list(...)
            self$params <- c(self$params, new_params)
            invisible(self)
            },
            
            fit = function(X, y) {
            self$feature_names <- colnames(X)
            
            # use ranger's x/y
            ranger_params <- list(
                x = as.data.frame(X),
                y = y,
                num.threads = parallel::detectCores() - 1,
                seed = self$random_state,
                importance = "impurity"  # sklearn's feature_importances_
            )
            
            # save parameters
            if (!is.null(self$params$n_estimators)) {
                ranger_params$num.trees <- self$params$n_estimators
            }
            if (!is.null(self$params$max_depth)) {
                ranger_params$max.depth <- self$params$max_depth
            }
            if (!is.null(self$params$min_samples_split)) {
                ranger_params$min.node.size <- self$params$min_samples_split
            }
            
            self$model <- do.call(ranger, ranger_params)
            invisible(self)
            },
            
            predict = function(X) {
            predict(self$model, as.data.frame(X))$predictions
            }
        )
    )

    ## Define LinearRegression Class
    LinearRegression <- R6Class("LinearRegression",
        public = list(
            params = list(),
            model = NULL,
            
            set_params = function(...) {
                new_params <- list(...)
                self$params <- c(self$params, new_params)
                invisible(self)
            },
            
            fit = function(X, y) {
                data <- data.frame(y = y, X, check.names = FALSE)
                safe_names <- paste0("`", names(X), "`")
                formula <- as.formula(paste("y ~", paste(safe_names, collapse = " + ")))
                all_params <- c(list(formula = formula, data = data), self$params)
                self$model <- do.call(lm, all_params)
                invisible(self)
            },
            
            predict = function(X) {
                predict(self$model, X)
            }
        )
    )

    for (i in 1:length(model_params)) {
        model_name <- model_params[[i]][['model_name']]
        
        if (model_name == "rf") {
            model <- RandomForestRegressor$new(random_state = random_seed)
        } else {
            model <- LinearRegression$new()
        }
        
        params <- model_params[[i]][['best params']]
        
        if (length(params) > 0) {
            params <- remove_model_in_param_grid(params)
            print(params)
            do.call(model$set_params, params)
        }
        
        all_models[["coding_all"]][[model_name]] <- model
    }

    cat("Training data dimensions:\n")
    cat("  Features (X):", paste(dim(rnaseq_wgc_all_X), collapse = " x "), "\n")
    cat("  Targets (y):", length(rnaseq_wgc_all_y), "\n")

    ## Begin Train
    results <- list()
    for (frag_type in c("mixed")) {
        print(frag_type)
        
        models <- list(
            rf = all_models[[gene_select_name]][["rf"]],
            lr = all_models[[gene_select_name]][["lr"]]
        )
        
        for (model_name in names(models)) {
            model <- models[[model_name]]
            
            # RF Parallel
            if (model_name == "rf") {
            model$params$num.threads <- n_cores
            }
            
            model$fit(rnaseq_wgc_all_X, rnaseq_wgc_all_y)
            rnaseq_wgc_pred_y <- model$predict(rnaseq_wgc_all_X)
            
            rmse <- sqrt(mean((rnaseq_wgc_all_y - rnaseq_wgc_pred_y)^2))
            cor_val <- cor(rnaseq_wgc_all_y, rnaseq_wgc_pred_y)^2
            
            cat("Model Name:\n")
            print(model_name)
            print(paste("Root Mean Squared Error:", rmse))
            print(paste("Pearson^2 for all data:", cor_val))
            
            all_models[[gene_select_name]][[model_name]] <- model
            
            results <- c(results, list(data.frame(
                model_name = model_name,
                RMSE = rmse,
                `Pearson^2` = cor_val,
                check.names = FALSE
            )))
        }
    }
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(save_dir, "/model_performance_results.csv"), row.names = FALSE)

    # Features
    top_features_dict <- list()
    top_features_dict[[gene_select_name]] <- list()

    print(gene_select_name)

    model <- all_models[[gene_select_name]][["rf"]]

    # Get feature importance from the fitted model
    feature_importance <- data.frame(
        feature = names(model$model$variable.importance),
        importance = model$model$variable.importance
    )
    ## Normalize
    feature_importance$importance <- feature_importance$importance / sum(feature_importance$importance)
    feature_importance <- feature_importance[order(-feature_importance$importance), ]
    head(feature_importance)
    write.csv(feature_importance, paste0(save_dir, "/top_features.csv"), row.names = FALSE)

    # Top Features
    top_features_dict[[gene_select_name]][[model_design]] <- feature_importance$feature
    write_json(top_features_dict, paste0(save_dir, "/top_feature_importance.json"), pretty = TRUE, auto_unbox = TRUE)
    pdp_lines_result <- list()
    pdp_lines_result[[gene_select_name]] <- list()

    model <- all_models[[gene_select_name]][["rf"]]
    pdp_lines_dict <- list()
    target_pairs <- top_features_dict[[gene_select_name]][[model_design]]

    if (remove_TV) {
        target_pairs <- target_pairs[!grepl("^(T_|V_)", target_pairs)]
    }
    head(target_pairs)

    # Calculate PDP - Parallel
    if (n_cores > 1 && length(target_pairs) > 1) {
        cat("Computing PDP in parallel using", n_cores, "cores...\n")
        
        cl <- makeCluster(min(n_cores, length(target_pairs)))
        registerDoParallel(cl)
        
        clusterExport(cl, c("model", "rnaseq_wgc_all_X"), envir = environment())
        
        pdp_results_list <- foreach(
            target_pair = target_pairs,
            .packages = c("pdp", "ranger"),
            .combine = 'c'
        ) %dopar% {
            
            # PDP average
            pdp_result <- partial(
                object = model$model,
                pred.var = target_pair,
                train = rnaseq_wgc_all_X,
                grid.resolution = 85
            )
            
            # ICE (individual)
            ice_result <- partial(
                object = model$model,
                pred.var = target_pair,
                train = rnaseq_wgc_all_X,
                grid.resolution = 85,
                ice = TRUE
            )
            
            grid_values <- unique(pdp_result[[target_pair]])
            average_values <- pdp_result$yhat
            
            n_grid <- length(grid_values)
            n_samples <- nrow(rnaseq_wgc_all_X)
            
            individual_matrix <- matrix(ice_result$yhat, nrow = n_samples, ncol = n_grid, byrow = TRUE)
            
            result <- list(
                average = list(average_values),
                individual = list(list(individual_matrix)),
                grid_values = list(grid_values)
            )
            
            setNames(list(result), target_pair)
        }
        
        stopCluster(cl)
        
        pdp_lines_result[[gene_select_name]] <- pdp_results_list
        
    } else {
        # Loop Version
        for (target_pair in target_pairs) {
            # PDP average
            pdp_result <- partial(
                object = model$model,
                pred.var = target_pair,
                train = rnaseq_wgc_all_X,
                grid.resolution = 85
            )
            
            # ICE (individual)
            ice_result <- partial(
                object = model$model,
                pred.var = target_pair,
                train = rnaseq_wgc_all_X,
                grid.resolution = 85,
                ice = TRUE
            )
            
            grid_values <- unique(pdp_result[[target_pair]])
            average_values <- pdp_result$yhat
            
            # ICE Resuts 
            n_grid <- length(grid_values)
            n_samples <- nrow(rnaseq_wgc_all_X)
            
            individual_matrix <- matrix(ice_result$yhat, nrow = n_samples, ncol = n_grid, byrow = TRUE)
            
            pdp_lines_result[[gene_select_name]][[target_pair]] <- list(
                average = list(average_values),  
                individual = list(list(individual_matrix)), 
                grid_values = list(grid_values)
            )
        }
    }

    saveRDS(pdp_lines_result, paste0(save_dir, "/top_features_pdp.rds"))
}
