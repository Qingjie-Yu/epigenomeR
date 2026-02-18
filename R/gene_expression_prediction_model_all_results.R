# Main Function
# Performs comprehensive PDP (Partial Dependence Plot) analysis pipeline for gene expression prediction models.
# Description: This function executes a complete machine learning workflow including:
#              1) Data preprocessing and feature normalization
#              2) Random Forest and Linear Regression model training with cross-validation
#              3) K-means clustering of gene features
#              4) Feature importance calculation
#              5) Generation of multiple visualization outputs (bar plots, PDP plots, summary plots)
#              6) Effect size calculation for model interpretation
# 
# Parameters:  gene_select_dir_filename: Path to JSON file containing gene selection dictionaries
#                                        (must include 'coding_cpg' and 'coding_non_cpg' keys).
#              target_pair_mapping_file_36x: Path to space-delimited file mapping target pairs.
#              target_pairs_selected_file: Path to JSON file with selected target pairs.
#              cluster_result: Path to CSV file containing pre-computed cluster assignments.
#              model_design: String identifier for model design (used for RNA-seq cutoff selection).
#              frag_type: Fragment type identifier (e.g., "mixed").
#              gene_select_name: Gene selection category name (e.g., "coding_all").
#              rnaseq_dir_filename: Path to CSV file containing RNA-seq expression data.
#              wgc_dir_filename: Path to `.feather` file containing WGC (whole genome coverage) features.
#              feature_path: Path to JSON file listing feature names to use.
#              rnaseq_cutoff_path: Path to JSON file with RNA-seq quantile cutoffs per model design.
#              model_params_file: Path to JSON file containing model hyperparameters.
#              save_dir: Output directory for all results and plots.
#              feature_colors_path: Path to CSV file mapping features to colors for visualization.
#              apply_par_and_pdp_output_dir: Directory for PDP computation outputs.
#              color_list: Named list of colors for cluster visualization.
#              target_pairs_list: List of target pairs for single cluster PDP plots.
#              effect_size_method: Method for effect size calculation (e.g., "cohen_d").
#              random_seed: Random seed for reproducibility (default: 42).
#              all_seeds: Vector of seeds for ensemble model training 
#                         (default: c(0, 1, 7, 42, 123, 999, 1234, 1337, 2021, 31415)).
#              top_num_features: Number of top features to display in importance plots (default: 20).
#              n_cores: Number of CPU cores for parallel processing (default: NULL, uses 1 core).
#              kmeans_random_seed: Random seed specifically for k-means clustering (default: NULL).
#              max_clusters: Maximum number of clusters for k-means (default: NULL).
#              summary_pdp_width: Width of summary PDP plot in inches (default: NULL).
#              summary_pdp_height: Height of summary PDP plot in inches (default: NULL).
# 
# Output: Generates multiple outputs in `save_dir`:
#          - "cv_plots/": Cross-validation scatter plots with density coloring.
#          - "cv_plots/cv_summary_*.csv": Cross-validation metrics summary.
#          - "clusters/": K-means clustering results.
#          - Feature importance bar plot (Figure 1).
#          - PDP grid plot (Figure 2).
#          - "sqrt_coding_all_centered_same_range/": Single cluster PDP plots (Figure 3).
#          - Summary PDP plot (Figure 4).
#          - Effect size CSV file (Data output).
# Dependencies: jsonlite, readr, arrow, dplyr, caret, ggplot2, MASS, viridis, 
#               Metrics, R6, randomForest, ranger, parallel, foreach, doParallel
draw_pdp_results <- function(gene_select_dir_filename, target_pair_mapping_file_36x, target_pairs_selected_file, cluster_result,
                    model_design, frag_type, gene_select_name, rnaseq_dir_filename, wgc_dir_filename, feature_path, 
                    rnaseq_cutoff_path, model_params_file, save_dir, feature_colors_path, apply_par_and_pdp_output_dir, 
                    color_list, target_pairs_list, effect_size_method,
                    random_seed = 42, all_seeds = c(0, 1, 7, 42, 123, 999, 1234, 1337, 2021, 31415), top_num_features = 20,
                    n_cores = NULL, kmeans_random_seed = NULL,
                    max_clusters = NULL, summary_pdp_width = NULL, summary_pdp_height = NULL) {
    # Load Libraries
    suppressPackageStartupMessages({
        library(jsonlite)
        library(readr)
        library(arrow)
        library(dplyr)
        library(caret)
        library(ggplot2)
        library(MASS)  # for kde2d
        library(viridis)
        library(Metrics)  # for rmse
        library(R6)
        library(randomForest)
        library(ranger)
        library(parallel)
        library(foreach)
        library(doParallel)
    }) 

    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

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

    gene_select_dict <- fromJSON(gene_select_dir_filename)
    gene_select_dict[["coding_all"]] <- c(gene_select_dict[["coding_cpg"]], gene_select_dict[["coding_non_cpg"]])

    target_pair_mapping_df <- read.csv(target_pair_mapping_file_36x, sep = " ")
    head(target_pair_mapping_df)
    # target_pair_mapping_df <- read_delim(target_pair_mapping_file_36x, delim = " ")

    target_pairs <- fromJSON(target_pairs_selected_file)

    cluster_result <- read.csv(cluster_result, row.names = 1)
    rownames(cluster_result) <- cluster_result$gene_id
    cluster_result <- cluster_result[, !names(cluster_result) %in% "gene_id", drop = FALSE]

    data <- list()
    rnaseq <- read.csv(rnaseq_dir_filename, sep = ",", header = TRUE, row.names = 1)
    rnaseq$sqrt_V <- log10(rnaseq$V1V2 + 1)

    wgc_raw <- read_feather(wgc_dir_filename)

    features <- fromJSON(feature_path)
    wgc_vals <- wgc_raw[, features, drop = FALSE]  

    zero_hiplex_genes <- wgc_raw$pos[rowSums(wgc_vals) == 0]
    wgc_raw[, features] <- log10(wgc_raw[, features] + 1)

    for (feature in features) {
        min_val <- min(wgc_raw[[feature]], na.rm = TRUE)
        max_val <- max(wgc_raw[[feature]], na.rm = TRUE)
        wgc_raw[[feature]] <- (wgc_raw[[feature]] - min_val) / (max_val - min_val)
    }

    features <- c("pos", features)
    wgc_raw <- wgc_raw[, features, drop = FALSE]
    wgc_raw_list <- wgc_raw$pos
    overlap_gene_list <- intersect(wgc_raw_list, gene_select_dict[["coding_all"]])
    rnaseq_avail <- rnaseq[overlap_gene_list, c("gene_id", "sqrt_V"), drop = FALSE]

    rnaseq_cutoffs <- fromJSON(rnaseq_cutoff_path)
    q99 <- rnaseq_cutoffs[[model_design]]
    rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, , drop = FALSE]
    zero_rnaseq_genes <- rnaseq_avail$gene_id[rnaseq_avail$sqrt_V == 0]
    zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)

    rnaseq_wgc_raw <- merge(rnaseq_avail, wgc_raw, by.x = "gene_id", by.y = "pos", all = FALSE)  # inner join
    rownames(rnaseq_wgc_raw) <- rnaseq_wgc_raw$gene_id
    rnaseq_wgc_raw <- rnaseq_wgc_raw[!rownames(rnaseq_wgc_raw) %in% zero_all_genes, , drop = FALSE]
    rnaseq_wgc <- rnaseq_wgc_raw[, !names(rnaseq_wgc_raw) %in% "pos", drop = FALSE]

    set.seed(random_seed)
    train_indices <- createDataPartition(rnaseq_wgc$sqrt_V, p = 0.8, list = FALSE)
    rnaseq_wgc_train <- rnaseq_wgc[train_indices, , drop = FALSE]
    rnaseq_wgc_test <- rnaseq_wgc[-train_indices, , drop = FALSE]

    rnaseq_wgc_train_X <- rnaseq_wgc_train[, !names(rnaseq_wgc_train) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_train_y <- rnaseq_wgc_train$sqrt_V
    rnaseq_wgc_test_X <- rnaseq_wgc_test[, !names(rnaseq_wgc_test) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_test_y <- rnaseq_wgc_test$sqrt_V
    rnaseq_wgc_all_X <- rnaseq_wgc[, !names(rnaseq_wgc) %in% "sqrt_V", drop = FALSE]
    rnaseq_wgc_all_y <- rnaseq_wgc$sqrt_V

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

    # CV Save Dir
    save_cv_dir <- paste0(save_dir, "/cv_plots")
    dir.create(save_cv_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Remove gene_id for modeling
    rnaseq_wgc_all_X_numeric <- rnaseq_wgc_all_X[, !names(rnaseq_wgc_all_X) %in% "gene_id", drop = FALSE]
    
    # Question: despite mixed, any other variable, keep this par?
    for (frag_type in c("mixed")) {
        print(frag_type)
        
        models <- list(
            "random_forest" = all_models[[gene_select_name]][["rf"]],
            "lr" = all_models[[gene_select_name]][["lr"]]
        )
        
        set.seed(random_seed)
        folds <- createFolds(rnaseq_wgc_all_y, k = 5, list = TRUE, returnTrain = FALSE)
        
        cv_scores <- list()

        # Clean column names once before k-fold
        names(rnaseq_wgc_all_X_numeric) <- make.names(names(rnaseq_wgc_all_X_numeric), unique = TRUE)

        # k-fold
        for (model_name in names(models)) {
            n_folds <- length(folds)
            rmse_scores <- numeric(n_folds)
            cors <- numeric(n_folds)
            spearman_cors <- numeric(n_folds)
            
            for (i in seq_along(folds)) {
                fold_name <- names(folds)[i]
                test_index <- folds[[fold_name]]
                train_index <- setdiff(1:nrow(rnaseq_wgc_all_X_numeric), test_index)
                
                X_train <- rnaseq_wgc_all_X_numeric[train_index, , drop = FALSE]
                X_test <- rnaseq_wgc_all_X_numeric[test_index, , drop = FALSE]
                y_train <- rnaseq_wgc_all_y[train_index]
                y_test <- rnaseq_wgc_all_y[test_index]
                
                # train model
                if (model_name == "random_forest") {
                    rf_model <- ranger(y ~ ., data = cbind(y = y_train, X_train), num.trees = 100, num.threads = n_cores)
                    y_pred <- predict(rf_model, X_test)$predictions
                } else {
                    lm_model <- lm(y ~ ., data = cbind(y = y_train, X_train))
                    y_pred <- predict(lm_model, X_test)
                }
                
                # Scatter Plot With Density Color
                # Compute 2D Density
                if (length(y_pred) > 1 && length(y_test) > 1) {
                    tryCatch({
                        density_est <- kde2d(y_pred, y_test, n = 50)
                        
                        density_values <- sapply(1:length(y_pred), function(j) {
                            x_idx <- which.min(abs(density_est$x - y_pred[j]))
                            y_idx <- which.min(abs(density_est$y - y_test[j]))
                            density_est$z[x_idx, y_idx]
                        })
                        
                        plot_data <- data.frame(y_pred = y_pred, y_test = y_test, density = density_values)
                        
                        p <- ggplot(plot_data, aes(x = y_pred, y = y_test, color = density)) +
                            geom_point(size = 2) +
                            scale_color_viridis_c() +
                            labs(
                                x = "y_pred", 
                                y = "y_test", 
                                title = paste(model_name, fold_name)) +
                            theme_bw()
                        
                    }, error = function(e) {
                        # If Density Approx Fail, Create Scatter Plot
                        plot_data <- data.frame(y_pred = y_pred, y_test = y_test)
                        p <- ggplot(plot_data, aes(x = y_pred, y = y_test)) +
                            geom_point(size = 2) +
                            labs(
                                x = "y_pred", 
                                y = "y_test", 
                                title = paste(model_name, fold_name)) +
                            theme_bw()
                    })
                    
                    plot_filename <- paste0(save_cv_dir, "/", model_name, "_", fold_name, ".png")
                    ggsave(plot_filename, p, width = 8, height = 6, dpi = 300)
                    print(p)
                }
                
                # Compute Criteria
                rmse_scores[i] <- sqrt(mean((y_test - y_pred)^2))
                cors[i] <- cor(y_test, y_pred, method = "pearson")
                spearman_cors[i] <- cor(y_test, y_pred, method = "spearman")
            }
            
            # Calculate Mean
            average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
            average_rmse <- mean(rmse_scores, na.rm = TRUE)
            average_cor <- mean(cors, na.rm = TRUE)
            
            cat(model_name, "\n")
            cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
            cat("Average RMSE:", average_rmse, "\n")
            cat("Pearson for each fold:", cors, "\n")
            cat("Average Pearson:", average_cor, "\n")
            cat("Average Spearman:", average_spearman_cor, "\n")
        }
    }

    # Get Models
    models <- list(
        "random_forest" = all_models[[gene_select_name]][["rf"]],
        "lr" = all_models[[gene_select_name]][["lr"]]
    )

    set.seed(random_seed)
    # Create index for 5-fold validation
    n_samples <- nrow(rnaseq_wgc_all_X_numeric)
    indices <- 1:n_samples

    # Disrupt Index
    set.seed(random_seed)
    shuffled_indices <- sample(indices)

    # 5-fold
    k <- 5
    fold_size <- floor(n_samples / k)
    folds <- list()

    for (i in 1:k) {
        start_idx <- (i - 1) * fold_size + 1
        if (i == k) {
            end_idx <- length(shuffled_indices)
        } else {
            end_idx <- i * fold_size
        }
        folds[[paste0("fold_", i)]] <- shuffled_indices[start_idx:end_idx]
    }

    cv_scores <- list()

    # Clean column names once before k-fold
    names(rnaseq_wgc_all_X_numeric) <- make.names(names(rnaseq_wgc_all_X_numeric), unique = TRUE)

    # k-fold validation
    for (model_name in names(models)) {
        cat("Processing model:", model_name, "\n")
        
        n_folds <- length(folds)
        
        if (n_cores > 1 && n_folds > 1) {
            # Parallel Version
            cat("Computing CV in parallel using", n_cores, "cores...\n")
            
            cl <- makeCluster(min(n_cores, n_folds))
            registerDoParallel(cl)
            
            clusterExport(cl, c("rnaseq_wgc_all_X_numeric", "rnaseq_wgc_all_y", "folds", "model_name"), envir = environment())
            
            fold_results <- foreach(
                i = seq_along(folds),
                .packages = c("randomForest"),
                .combine = 'rbind'
            ) %dopar% {
                fold_name <- names(folds)[i]
                test_index <- folds[[fold_name]]
                train_index <- setdiff(1:nrow(rnaseq_wgc_all_X_numeric), test_index)
                
                X_train <- rnaseq_wgc_all_X_numeric[train_index, , drop = FALSE]
                X_test <- rnaseq_wgc_all_X_numeric[test_index, , drop = FALSE]
                y_train <- rnaseq_wgc_all_y[train_index]
                y_test <- rnaseq_wgc_all_y[test_index]
                
                tryCatch({
                    if (model_name == "random_forest") {
                        train_data <- cbind(y = y_train, X_train)
                        model <- randomForest(y ~ ., data = train_data, ntree = 100)
                        y_pred <- predict(model, X_test)
                    } else if (model_name == "lr") {
                        train_data <- cbind(y = y_train, X_train)
                        model <- lm(y ~ ., data = train_data)
                        y_pred <- predict(model, X_test)
                    }
                    
                    if (is.matrix(y_pred)) y_pred <- as.vector(y_pred)
                    
                    data.frame(
                        rmse = sqrt(mean((y_test - y_pred)^2, na.rm = TRUE)),
                        cor = cor(y_test, y_pred, method = "pearson", use = "complete.obs"),
                        spearman = cor(y_test, y_pred, method = "spearman", use = "complete.obs")
                    )
                }, error = function(e) {
                    data.frame(rmse = NA, cor = NA, spearman = NA)
                })
            }
            
            stopCluster(cl)
            
            rmse_scores <- fold_results$rmse
            cors <- fold_results$cor
            spearman_cors <- fold_results$spearman
            
        } else {
            # Sequential loop version
            rmse_scores <- numeric(n_folds)
            cors <- numeric(n_folds)
            spearman_cors <- numeric(n_folds)
            
            for (i in seq_along(folds)) {
                fold_name <- names(folds)[i]
                cat("  Fold", i, "/", n_folds, "\n")
                
                test_index <- folds[[fold_name]]
                train_index <- setdiff(1:nrow(rnaseq_wgc_all_X_numeric), test_index)
                
                X_train <- rnaseq_wgc_all_X_numeric[train_index, , drop = FALSE]
                X_test <- rnaseq_wgc_all_X_numeric[test_index, , drop = FALSE]
                y_train <- rnaseq_wgc_all_y[train_index]
                y_test <- rnaseq_wgc_all_y[test_index]
                
                tryCatch({
                    if (model_name == "random_forest") {
                        train_data <- cbind(y = y_train, X_train)
                        model <- randomForest(y ~ ., data = train_data, ntree = 100)
                        y_pred <- predict(model, X_test)
                    } else if (model_name == "lr") {
                        train_data <- cbind(y = y_train, X_train)
                        model <- lm(y ~ ., data = train_data)
                        y_pred <- predict(model, X_test)
                    }
                    
                    if (is.matrix(y_pred)) y_pred <- as.vector(y_pred)
                    
                    rmse_scores[i] <- sqrt(mean((y_test - y_pred)^2, na.rm = TRUE))
                    cors[i] <- cor(y_test, y_pred, method = "pearson", use = "complete.obs")
                    spearman_cors[i] <- cor(y_test, y_pred, method = "spearman", use = "complete.obs")
                    
                    cat("    RMSE:", rmse_scores[i], "Pearson:", cors[i], "\n")
                    
                }, error = function(e) {
                    cat("    Error:", e$message, "\n")
                    rmse_scores[i] <<- NA
                    cors[i] <<- NA
                    spearman_cors[i] <<- NA
                })
            }
        }
        
        # Calculate means
        average_rmse <- mean(rmse_scores, na.rm = TRUE)
        average_cor <- mean(cors, na.rm = TRUE)
        average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
        
        # Print results
        cat("\n")
        cat(model_name, "\n")
        cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
        cat("Average RMSE:", average_rmse, "\n")
        cat("Pearson for each fold:", cors, "\n")
        cat("Average Pearson:", average_cor, "\n")
        cat("Average Spearman:", average_spearman_cor, "\n\n")
        
        # Store results
        cv_scores[[model_name]] <- list(
            rmse_scores = rmse_scores,
            average_rmse = average_rmse,
            cors = cors,
            average_cor = average_cor,
            spearman_cors = spearman_cors,
            average_spearman_cor = average_spearman_cor
        )
    }

    cat("=== Cross-Validation Results Summary ===\n")
    for (model_name in names(cv_scores)) {
        results <- cv_scores[[model_name]]
        cat(sprintf("%s: RMSE=%.4f, Pearson=%.4f, Spearman=%.4f\n", 
                    model_name, results$average_rmse, results$average_cor, results$average_spearman_cor))
    }

    results_summary <- data.frame(
        model = names(cv_scores),
        average_rmse = sapply(cv_scores, function(x) x$average_rmse),
        average_pearson = sapply(cv_scores, function(x) x$average_cor),
        average_spearman = sapply(cv_scores, function(x) x$average_spearman_cor),
        stringsAsFactors = FALSE
    )

    write.csv(results_summary, file = paste0(save_cv_dir, "/cv_summary_", gene_select_name, "_", frag_type, ".csv"), row.names = FALSE)


    # Kmeans Clusters
    save_cluster_dir <- paste0(save_dir, "/clusters")
    generate_kmeans_cluster(wgc_raw = rnaseq_wgc_all_X, save_cluster_dir = save_cluster_dir, rnaseq_dir_filename = rnaseq_dir_filename, kmeans_random_seed = kmeans_random_seed, max_clusters = max_clusters)

    # Calculate Feature Importance
    results <- calculate_feature_importance(save_cluster_dir = save_cluster_dir, result_dir = apply_par_and_pdp_output_dir, target_pair_mapping_df = target_pair_mapping_df, top_num_features = top_num_features, random_seeds = all_seeds)

    wgc_raw <- results$wgc_raw
    cluster_ids <- results$cluster_ids
    top_features <- results$top_features
    feature_importance <- results$feature_importance   

    # Read Feature Colors File
    feature_colors_df <- read.csv(feature_colors_path, row.names = 1)

    # Output 1 - Bar Plot
    # Figure 1
    feature_importance_bar_plot(feature_colors_df = feature_colors_df, feature_importance = feature_importance, save_bar_plot_dir = save_dir)

    # Output 2 - PDP Grid Plot
    # Figure 2
    results <- pdp_grid_plot(wgc_raw = wgc_raw, top_features = top_features, cluster_ids = cluster_ids, apply_par_and_pdp_output_dir = apply_par_and_pdp_output_dir, gene_select_dict = gene_select_dict, gene_select_name = gene_select_name, save_pdp_grid_dir = save_dir, target_pair_mapping_df = target_pair_mapping_df, random_seeds = all_seeds)

    pdp_results <- results$pdp_results
    cluster_results_stacked <- results$cluster_results_stacked

    # Output 3 - Draw PDP Plot For Chosen Epitope Pairs In Each Cluster
    # Figure 3
    save_single_cluster_dir <- paste0(save_dir, "/sqrt_coding_all_centered_same_range")

    single_cluster_pdp(save_single_cluster_dir = save_single_cluster_dir, cluster_ids = cluster_ids, color_list = color_list, target_pairs_list = target_pairs_list, pdp_results = pdp_results, cluster_results_stacked = cluster_results_stacked)

    # Output 4 - Summary PDP Plot
    # Figure 4
    summary_pdp_plot(top_features = top_features, pdp_results = pdp_results, cluster_results_stacked = cluster_results_stacked, feature_colors_df = feature_colors_df, save_summary_pdp_dir = save_dir, width = summary_pdp_width, height = summary_pdp_height)

    # Output 5 - effect size calculation
    # Data 1: csv file
    effect_size_calculation(top_features = top_features,  pdp_results = pdp_results, cluster_results_stacked = cluster_results_stacked, target_pair_mapping_df = target_pair_mapping_df, cluster_ids = cluster_ids, gene_select_name = gene_select_name, save_csv_dir = save_dir, effect_size_method = effect_size_method, random_seeds = all_seeds)

}
