


# helper functions
map_target_names_pdp <- function(target_pair_list, target_pair_mapping_df, from_col = "targets", to_col = "shorthand") {
  cur_names <- target_pair_mapping_df[[from_col]]
  new_names <- target_pair_mapping_df[[to_col]]
  
  cat("\n")
  result <- target_pair_list
  for (i in seq_along(cur_names)) {
    cur_name <- cur_names[i]
    new_name <- new_names[i]
    result <- sapply(result, function(target_pair) {
      gsub(cur_name, new_name, target_pair, fixed = TRUE)
    }, USE.NAMES = FALSE)
  }
  
  return(result)
}

column_to_rownames <- function(wgc, var = "pos") {
  pos_list_full <- wgc[[var]]
  rownames(wgc) <- wgc[[var]]
  wgc <- wgc[, !names(wgc) %in% var, drop = FALSE]
  return(wgc)
}

remove_model_in_param_grid <- function(best_params) {
  results <- list()
  for (param in names(best_params)) {
    new_param <- gsub("model__", "", param, fixed = TRUE)
    results[[new_param]] <- best_params[[param]]
  }
  return(results)
}

safe_minmax <- function(df) {
  df_scaled <- df
  for (col in colnames(df)) {
    col_min <- min(df[[col]], na.rm = TRUE)
    col_max <- max(df[[col]], na.rm = TRUE)
    if (col_max != col_min) {
      df_scaled[[col]] <- (df[[col]] - col_min) / (col_max - col_min)
    } else {
      df_scaled[[col]] <- 0.0
    }
  }
  return(df_scaled)
}

gene_select_name <- "coding_all"
model_design <- "rnaseq_vs_hiplex_rm_outlier_log"
random_seed <- 42


library(arrow)
library(caret)

gene_select_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island/gene_categories.json"
library(jsonlite)
gene_select_dict <- fromJSON(gene_select_dir_filename)
gene_select_dict[["coding_all"]] <- c(gene_select_dict[["coding_cpg"]], gene_select_dict[["coding_non_cpg"]])

print(gene_select_name)
gene_select_list <- gene_select_dict[[gene_select_name]]

rnaseq_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/RNA_seq_TPM_all.csv"
rnaseq <- read.csv(rnaseq_dir_filename, header = TRUE, row.names = 1)
rnaseq$sqrt_V <- log10(rnaseq$V1V2 + 1)

wgc_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/promoter_-1000-1000/V_mixed_promoter_-1000-1000_colQC-all_libnorm.feather"
wgc_raw <- read_feather(wgc_dir_filename)

features_path <- "/dcs05/hongkai/data/next_cutntag/script/utils/filtered_target_pairs.json"
features <- fromJSON(features_path)

wgc_vals <- wgc_raw[, features]
zero_hiplex_genes <- wgc_raw$pos[rowSums(wgc_vals) == 0]
wgc_raw[, features] <- log10(wgc_raw[, features] + 1)
wgc_raw[, features] <- apply(wgc_raw[, features], 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
})

features <- c("pos", features)
wgc_raw <- wgc_raw[, features, drop = FALSE]

wgc_raw_list <- wgc_raw$pos
overlap_gene_list <- intersect(wgc_raw_list, gene_select_dict[["coding_all"]])
rnaseq_avail <- rnaseq[overlap_gene_list, c("gene_id", "sqrt_V"), drop = FALSE]

rnaseq_cutoffs_path <- "/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_hiplex_cutoff.json"
rnaseq_cutoffs <- fromJSON(rnaseq_cutoffs_path)
q99 <- rnaseq_cutoffs[[model_design]]
rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, , drop = FALSE]
zero_rnaseq_genes <- rnaseq_avail$gene_id[rnaseq_avail$sqrt_V == 0]
zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)

rnaseq_wgc_raw <- merge(rnaseq_avail, wgc_raw, by.x = "gene_id", by.y = "pos", all = FALSE)  # inner join
rownames(rnaseq_wgc_raw) <- rnaseq_wgc_raw$gene_id
rnaseq_wgc_raw <- rnaseq_wgc_raw[!rownames(rnaseq_wgc_raw) %in% zero_all_genes, , drop = FALSE]
rnaseq_wgc <- rnaseq_wgc_raw[, !names(rnaseq_wgc_raw) %in% c("pos", "gene_id"), drop = FALSE]

set.seed(42)
n <- nrow(rnaseq_wgc)
test_size <- 0.2
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




library(jsonlite)
library(randomForest)
library(caret)
library(ggplot2)
library(MASS)  # for kde2d
library(viridis)
library(Metrics)  # for rmse


frag_type <- "mixed"
all_models <- list()
all_models[["coding_all"]] <- list()

model_params_file <- "/dcs05/hongkai/data/next_cutntag/script/explainability/rnaseq_vs_hiplex_rm_outlier_log/coding_all.json"
random_seed <- 42
model_params <- fromJSON(model_params_file, simplifyDataFrame = FALSE)

library(R6)
library(randomForest)

library(R6)
library(ranger)

# replace RandomForestRegressor class - use x/y avoid formula
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
                                     
                                     # 使用 ranger 的 x/y 接口，避免 formula 问题
                                     ranger_params <- list(
                                       x = as.data.frame(X),
                                       y = y,
                                       num.threads = parallel::detectCores() - 1,
                                       seed = self$random_state,
                                       importance = "impurity"  # 添加这行！对应 sklearn 的 feature_importances_
                                     )
                                     
                                     # 转换 sklearn 参数名到 ranger
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

# Define LinearRegression Class
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
                                # 用反引号包裹列名，处理特殊字符
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

print(gene_select_name)

dim(rnaseq_wgc_all_X)
length(rnaseq_wgc_all_y)

save_dir <- paste0("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step2_output/", random_seed)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()

library(parallel)
library(doParallel)

# parralel core
n_cores <- detectCores() - 1 
cat("Using", n_cores, "cores\n")

results <- list()

for (frag_type in c("mixed")) {
  print(frag_type)
  
  models <- list(
    rf = all_models[[gene_select_name]][["rf"]],
    lr = all_models[[gene_select_name]][["lr"]]
  )
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # RF 使用并行
    if (model_name == "rf") {
      # 方法1: 直接在 randomForest 中设置并行（推荐）
      model$params$num.threads <- n_cores
      # 或者如果用的是原生 randomForest，可以用 ranger 替代（更快）
    }
    
    model$fit(rnaseq_wgc_all_X, rnaseq_wgc_all_y)
    rnaseq_wgc_pred_y <- model$predict(rnaseq_wgc_all_X)
    
    rmse <- sqrt(mean((rnaseq_wgc_all_y - rnaseq_wgc_pred_y)^2))
    cor_val <- cor(rnaseq_wgc_all_y, rnaseq_wgc_pred_y)^2
    
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


library(jsonlite)

top_features_dict <- list()
top_features_dict[[gene_select_name]] <- list()

print(gene_select_name)

model <- all_models[[gene_select_name]][["rf"]]

# Get feature importance from the fitted model
feature_importance <- data.frame(
  feature = names(model$model$variable.importance),
  importance = model$model$variable.importance
)
## guiyi
feature_importance$importance <- feature_importance$importance / sum(feature_importance$importance)
feature_importance <- feature_importance[order(-feature_importance$importance), ]

write.csv(feature_importance, paste0(save_dir, "/top_features.csv"), row.names = FALSE)

top_features_dict[[gene_select_name]][[model_design]] <- feature_importance$feature

print(feature_importance)

write_json(top_features_dict, paste0(save_dir, "/top_feature_importance.json"), pretty = TRUE, auto_unbox = TRUE)

pdp_lines_result <- list()

pdp_lines_result[[gene_select_name]] <- list()




library(pdp)  

write.csv(rnaseq_wgc_all_X, paste0(save_dir, "/rnaseq_wgc_all_X.csv"), row.names = FALSE)

model <- all_models[[gene_select_name]][["rf"]]
pdp_lines_dict <- list()
target_pairs <- top_features_dict[[gene_select_name]][[model_design]]

remove_TV <- FALSE
if (remove_TV) {
  target_pairs <- target_pairs[!grepl("^(T_|V_)", target_pairs)]
}

print(target_pairs)

# replace PDP loop
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
  
  # ice_result 
  n_grid <- length(grid_values)
  n_samples <- nrow(rnaseq_wgc_all_X)
  
  individual_matrix <- matrix(ice_result$yhat, nrow = n_samples, ncol = n_grid, byrow = TRUE)
  
  pdp_lines_result[[gene_select_name]][[target_pair]] <- list(
    average = list(average_values),  
    individual = list(list(individual_matrix)), 
    grid_values = list(grid_values)
  )
}


saveRDS(pdp_lines_result, paste0(save_dir, "/top_features_pdp.rds"))