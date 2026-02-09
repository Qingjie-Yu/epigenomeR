


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

calculate_effect_size <- function(x, y, method = "slope") {  
    if (method == "slope") {
        result <- (y[length(y)] - y[1]) / (x[length(x)] - x[1])
        
    } else if (method == "auc") {
        sorted_indices <- order(x)
        x_sorted <- x[sorted_indices]
        y_sorted <- y[sorted_indices]
        result <- pracma::trapz(x_sorted, y_sorted)
    } else if (method == "end") {
        result <- y[length(y)]
    }
    return(result)
}


get_average_lines <- function(pdp_results, target_pair) {
    all_arrays <- list()
    
    # grid length
    first_seed <- names(pdp_results)[1]
    grid_length <- length(pdp_results[[first_seed]][["coding_all"]][[target_pair]][["grid_values"]][[1]])
    
    for (random_seed in names(pdp_results)) {
        individual_lines <- pdp_results[[random_seed]][["coding_all"]][[target_pair]][["individual"]]
        
        n_samples <- length(individual_lines) / grid_length
        
        if (length(individual_lines) %% grid_length == 0) {
            reshaped_matrix <- matrix(individual_lines, 
                                    nrow = n_samples, 
                                    ncol = grid_length, 
                                    byrow = TRUE)
            all_arrays <- append(all_arrays, list(reshaped_matrix))
        } else {
            warning(paste("Cannot reshape data for seed", random_seed))
        }
    }
    
    if (length(all_arrays) > 0) {
        stacked <- do.call(rbind, all_arrays)
        return(stacked)
    } else {
        stop("Could not process individual_lines data")
    }
}


get_average_lines <- function(pdp_results, target_pair) {
    all_arrays <- list()
    
    seed_names <- names(pdp_results)
    seed_names_sorted <- seed_names[order(as.numeric(seed_names))]
    
    cat("Processing:", seed_names_sorted, "\n")  
    for (random_seed in seed_names_sorted) {
        individual_lines <- pdp_results[[random_seed]][["coding_all"]][[target_pair]][["individual"]]
        grid_length <- length(pdp_results[[random_seed]][["coding_all"]][[target_pair]][["grid_values"]][[1]])
        
        n_samples <- length(individual_lines) / grid_length
        
        if (length(individual_lines) %% grid_length == 0) {
            reshaped_matrix <- matrix(individual_lines, 
                                    nrow = n_samples, 
                                    ncol = grid_length, 
                                    byrow = FALSE)
            
            all_arrays <- append(all_arrays, list(reshaped_matrix))
        }
    }
    
    stacked <- do.call(rbind, all_arrays)
    return(stacked)
}



gene_select_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island/gene_categories.json"

library(jsonlite)
gene_select_dict <- fromJSON(gene_select_dir_filename)

gene_select_dict[["coding_all"]] <- c(gene_select_dict[["coding_cpg"]], gene_select_dict[["coding_non_cpg"]])

target_pair_mapping_file_36x <- "/dcs05/hongkai/data/next_cutntag/script/utils/target_pair_short_hand.csv"
target_pair_mapping_df <- read.csv(target_pair_mapping_file_36x, sep = " ")

library(readr)
target_pair_mapping_df <- read_delim(target_pair_mapping_file_36x, delim = " ")

library(jsonlite)
target_pairs_selected_file <- "/dcs05/hongkai/data/next_cutntag/bulk/explainability/pdp_explain_results/target_pairs_selected.json"
target_pairs <- fromJSON(target_pairs_selected_file)

cluster_result <- read.csv("/dcs05/hongkai/data/next_cutntag/bulk/explainability/leave_one_out/V/all/coding_cpg_clusters=5.csv", row.names = 1)

rownames(cluster_result) <- cluster_result$gene_id
cluster_result <- cluster_result[, !names(cluster_result) %in% "gene_id", drop = FALSE]


library(arrow)
library(jsonlite)
library(dplyr)
library(caret)

data <- list()
model_designs <- c("rnaseq_vs_hiplex_rm_outlier_log")
model_design <- "rnaseq_vs_hiplex_rm_outlier_log"
frag_type <- "mixed"
gene_select_names <- c("coding_all")
gene_select_name <- "coding_all"


rnaseq_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/RNA_seq_TPM_all.csv"
rnaseq <- read.csv(rnaseq_dir_filename, sep = ",", header = TRUE, row.names = 1)

rnaseq$sqrt_V <- log10(rnaseq$V1V2 + 1)

wgc_dir_filename <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/promoter_-1000-1000/V_mixed_promoter_-1000-1000_colQC-all_libnorm.feather"
wgc_raw <- read_feather(wgc_dir_filename)

features <- fromJSON('/dcs05/hongkai/data/next_cutntag/script/utils/filtered_target_pairs.json')
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

rnaseq_cutoffs <- fromJSON("/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_hiplex_cutoff.json")
q99 <- rnaseq_cutoffs[[model_design]]
rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, , drop = FALSE]
zero_rnaseq_genes <- rnaseq_avail$gene_id[rnaseq_avail$sqrt_V == 0]
zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)

rnaseq_wgc_raw <- merge(rnaseq_avail, wgc_raw, by.x = "gene_id", by.y = "pos", all = FALSE)  # inner join
rownames(rnaseq_wgc_raw) <- rnaseq_wgc_raw$gene_id
rnaseq_wgc_raw <- rnaseq_wgc_raw[!rownames(rnaseq_wgc_raw) %in% zero_all_genes, , drop = FALSE]
rnaseq_wgc <- rnaseq_wgc_raw[, !names(rnaseq_wgc_raw) %in% "pos", drop = FALSE]

set.seed(42)
library(caret)
train_indices <- createDataPartition(rnaseq_wgc$sqrt_V, p = 0.8, list = FALSE)
rnaseq_wgc_train <- rnaseq_wgc[train_indices, , drop = FALSE]
rnaseq_wgc_test <- rnaseq_wgc[-train_indices, , drop = FALSE]

rnaseq_wgc_train_X <- rnaseq_wgc_train[, !names(rnaseq_wgc_train) %in% "sqrt_V", drop = FALSE]
rnaseq_wgc_train_y <- rnaseq_wgc_train$sqrt_V
rnaseq_wgc_test_X <- rnaseq_wgc_test[, !names(rnaseq_wgc_test) %in% "sqrt_V", drop = FALSE]
rnaseq_wgc_test_y <- rnaseq_wgc_test$sqrt_V
rnaseq_wgc_all_X <- rnaseq_wgc[, !names(rnaseq_wgc) %in% "sqrt_V", drop = FALSE]
rnaseq_wgc_all_y <- rnaseq_wgc$sqrt_V


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

# Define RandomForestRegressor Class
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
      data <- data.frame(y = y, X)
      formula <- as.formula(paste("y ~", paste(names(X), collapse = " + ")))
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


# for (frag_type in c('mixed')) {
#     print(frag_type)
    
#     models <- list(
#         "random_forest" = all_models[[gene_select_name]][["rf"]],
#         "lr" = all_models[[gene_select_name]][["lr"]]
#     )
    
#     set.seed(42)
#     folds <- createFolds(rnaseq_wgc_all_y, k = 5, list = TRUE, returnTrain = FALSE)
    
#     cv_scores <- list()
    
#     # 执行K折交叉验证
#     for (model_name in names(models)) {
#         rmse_scores <- c()
#         cors <- c()
#         spearman_cors <- c()
        
#         for (fold_name in names(folds)) {
#             test_index <- folds[[fold_name]]
#             # 在交叉验证之前，清理列名
#             clean_column_names <- function(df) {
#                 names(df) <- make.names(names(df), unique = TRUE)
#                 return(df)
#             }

#             # 清理数据的列名
#             rnaseq_wgc_all_X <- clean_column_names(rnaseq_wgc_all_X)

#             # 或者更简单的方式
#             names(rnaseq_wgc_all_X) <- make.names(names(rnaseq_wgc_all_X), unique = TRUE)

#             train_index <- setdiff(1:nrow(rnaseq_wgc_all_X), test_index)
            
#             X_train <- rnaseq_wgc_all_X[train_index, , drop = FALSE]
#             X_test <- rnaseq_wgc_all_X[test_index, , drop = FALSE]
#             y_train <- rnaseq_wgc_all_y[train_index]
#             y_test <- rnaseq_wgc_all_y[test_index]
            
#             # 训练模型
#             if (model_name == "random_forest") {
#                 # 创建训练数据
#                 train_data <- cbind(y = y_train, X_train)
#                 model <- randomForest(y ~ ., data = train_data, ntree = 100)
#                 y_pred <- predict(model, X_test)
#             } else {  # 线性回归
#                 train_data <- cbind(y = y_train, X_train)
#                 model <- lm(y ~ ., data = train_data)
#                 y_pred <- predict(model, X_test)
#             }
            
#             # 创建散点图with密度颜色
#             # 计算2D核密度
#             if (length(y_pred) > 1 && length(y_test) > 1) {
#                 tryCatch({
#                     # 创建密度估计
#                     density_est <- kde2d(y_pred, y_test, n = 50)
                    
#                     # 为每个点分配密度值
#                     density_values <- sapply(1:length(y_pred), function(i) {
#                         x_idx <- which.min(abs(density_est$x - y_pred[i]))
#                         y_idx <- which.min(abs(density_est$y - y_test[i]))
#                         density_est$z[x_idx, y_idx]
#                     })
                    
#                     # 创建数据框
#                     plot_data <- data.frame(
#                         y_pred = y_pred,
#                         y_test = y_test,
#                         density = density_values
#                     )
                    
#                     # 创建ggplot
#                     p <- ggplot(plot_data, aes(x = y_pred, y = y_test, color = density)) +
#                         geom_point(size = 2) +
#                         scale_color_viridis_c() +
#                         labs(
#                             x = "y_pred",
#                             y = "y_test",
#                             title = paste(model_name, fold_name)
#                         ) +
#                         theme_bw()
                    
#                     print(p)
                    
#                 }, error = function(e) {
#                     # 如果密度估计失败，创建简单散点图
#                     plot_data <- data.frame(y_pred = y_pred, y_test = y_test)
#                     p <- ggplot(plot_data, aes(x = y_pred, y = y_test)) +
#                         geom_point(size = 2) +
#                         labs(
#                             x = "y_pred",
#                             y = "y_test",
#                             title = paste(model_name, fold_name)
#                         ) +
#                         theme_bw()
#                     print(p)
#                 })
#             }
            
#             # 计算评估指标
#             rmse <- sqrt(mean((y_test - y_pred)^2))
#             cor_val <- cor(y_test, y_pred, method = "pearson")
#             spearman_cor <- cor(y_test, y_pred, method = "spearman")
            
#             spearman_cors <- c(spearman_cors, spearman_cor)
#             rmse_scores <- c(rmse_scores, rmse)
#             cors <- c(cors, cor_val)
#         }
        
#         # 计算平均值
#         average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
#         average_rmse <- mean(rmse_scores, na.rm = TRUE)
#         average_cor <- mean(cors, na.rm = TRUE)
        
#         cat(model_name, "\n")
#         cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
#         cat("Average Mean Squared Error:", average_rmse, "\n")
#         cat("Pearson for each fold:", cors, "\n")
#         cat("Average Mean Pearson:", average_cor, "\n")
#         cat("Average Mean Spearman:", average_spearman_cor, "\n\n")
#     }
# }

save_dir <- "/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/mine"

for (frag_type in c("mixed")) {
    print(frag_type)
    
    models <- list(
        "random_forest" = all_models[[gene_select_name]][["rf"]],
        "lr" = all_models[[gene_select_name]][["lr"]]
    )
    
    set.seed(42)
    folds <- createFolds(rnaseq_wgc_all_y, k = 5, list = TRUE, returnTrain = FALSE)
    
    cv_scores <- list()
    
    # save dir
    plot_dir <- paste0(save_dir, "/cv_plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # k-fold
    for (model_name in names(models)) {
        rmse_scores <- c()
        cors <- c()
        spearman_cors <- c()
        
        for (fold_name in names(folds)) {
            test_index <- folds[[fold_name]]
            # clear col names before k-fold
            clean_column_names <- function(df) {
                names(df) <- make.names(names(df), unique = TRUE)
                return(df)
            }

            # clear col names
            rnaseq_wgc_all_X <- clean_column_names(rnaseq_wgc_all_X)

            train_index <- setdiff(1:nrow(rnaseq_wgc_all_X), test_index)
            
            X_train <- rnaseq_wgc_all_X[train_index, , drop = FALSE]
            X_test <- rnaseq_wgc_all_X[test_index, , drop = FALSE]
            y_train <- rnaseq_wgc_all_y[train_index]
            y_test <- rnaseq_wgc_all_y[test_index]
            
            # train model
            if (model_name == "random_forest") {
                # train data
                train_data <- cbind(y = y_train, X_train)
                model <- randomForest(y ~ ., data = train_data, ntree = 100)
                y_pred <- predict(model, X_test)
            } else {  # linear reg
                train_data <- cbind(y = y_train, X_train)
                model <- lm(y ~ ., data = train_data)
                y_pred <- predict(model, X_test)
            }
            
            # calc 2d density
            if (length(y_pred) > 1 && length(y_test) > 1) {
                tryCatch({
                    density_est <- kde2d(y_pred, y_test, n = 50)
                    
                    density_values <- sapply(1:length(y_pred), function(i) {
                        x_idx <- which.min(abs(density_est$x - y_pred[i]))
                        y_idx <- which.min(abs(density_est$y - y_test[i]))
                        density_est$z[x_idx, y_idx]
                    })
                    
                    plot_data <- data.frame(
                        y_pred = y_pred,
                        y_test = y_test,
                        density = density_values
                    )
                    
                    p <- ggplot(plot_data, aes(x = y_pred, y = y_test, color = density)) +
                        geom_point(size = 2) +
                        scale_color_viridis_c() +
                        labs(
                            x = "y_pred",
                            y = "y_test",
                            title = paste(model_name, fold_name)
                        ) +
                        theme_bw()
                    
                    plot_filename <- paste0(plot_dir, "/", model_name, "_", fold_name, ".png")
                    ggsave(plot_filename, p, width = 8, height = 6, dpi = 300)
                    
                    print(p)
                    
                }, error = function(e) {
                    plot_data <- data.frame(y_pred = y_pred, y_test = y_test)
                    p <- ggplot(plot_data, aes(x = y_pred, y = y_test)) +
                        geom_point(size = 2) +
                        labs(
                            x = "y_pred",
                            y = "y_test",
                            title = paste(model_name, fold_name)
                        ) +
                        theme_bw()
                    
                    plot_filename <- paste0(plot_dir, "/", model_name, "_", fold_name, ".png")
                    ggsave(plot_filename, p, width = 8, height = 6, dpi = 300)
                    
                    print(p)
                })
            }
            
            rmse <- sqrt(mean((y_test - y_pred)^2))
            cor_val <- cor(y_test, y_pred, method = "pearson")
            spearman_cor <- cor(y_test, y_pred, method = "spearman")
            
            spearman_cors <- c(spearman_cors, spearman_cor)
            rmse_scores <- c(rmse_scores, rmse)
            cors <- c(cors, cor_val)
        }
        
        # cal mean
        average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
        average_rmse <- mean(rmse_scores, na.rm = TRUE)
        average_cor <- mean(cors, na.rm = TRUE)
        
        cat(model_name, "\n")
        cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
        cat("Average Mean Squared Error:", average_rmse, "\n")
        cat("Pearson for each fold:", cors, "\n")
        cat("Average Mean Pearson:", average_cor, "\n")
        cat("Average Mean Spearman:", average_spearman_cor, "\n")
    }
}


# library
library(jsonlite)
library(ggplot2)
library(viridis)
library(MASS)
library(caret)

print(frag_type)

models <- list(
    "random_forest" = all_models[[gene_select_name]][["rf"]],
    "lr" = all_models[[gene_select_name]][["lr"]]
)

set.seed(42)

# 5-fold index
n_samples <- nrow(rnaseq_wgc_all_X)
indices <- 1:n_samples

set.seed(42)
shuffled_indices <- sample(indices)

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

# # k-fold val
# for (model_name in names(models)) {
#     cat("Processing model:", model_name, "\n")
    
#     rmse_scores <- c()
#     cors <- c()
#     spearman_cors <- c()
    
#     model_obj <- models[[model_name]]
#     fold_counter <- 1
    
#     for (fold_name in names(folds)) {
#         cat("  Processing", fold_name, "\n")
        
#         test_index <- folds[[fold_name]]
#         train_index <- setdiff(1:nrow(rnaseq_wgc_all_X), test_index)
        
#         # 准备数据
#         X_train <- rnaseq_wgc_all_X[train_index, , drop = FALSE]
#         X_test <- rnaseq_wgc_all_X[test_index, , drop = FALSE]
#         y_train <- rnaseq_wgc_all_y[train_index]
#         y_test <- rnaseq_wgc_all_y[test_index]
        
#         # 清理列名
#         names(X_train) <- make.names(names(X_train), unique = TRUE)
#         names(X_test) <- make.names(names(X_test), unique = TRUE)
        
#         # 训练和预测
#         tryCatch({
#             # 训练模型
#             model_obj$fit(X_train, y_train)
            
#             # 预测
#             y_pred <- model_obj$predict(X_test)
            
#             # 确保是向量
#             if (is.matrix(y_pred)) {
#                 y_pred <- as.vector(y_pred)
#             }
            
#             # 计算评估指标
#             rmse_val <- sqrt(mean((y_test - y_pred)^2, na.rm = TRUE))
#             cor_val <- cor(y_test, y_pred, method = "pearson", use = "complete.obs")
#             spearman_val <- cor(y_test, y_pred, method = "spearman", use = "complete.obs")
            
#             # 存储结果 - 修复这里的语法
#             rmse_scores <- c(rmse_scores, rmse_val)
#             cors <- c(cors, cor_val)
#             spearman_cors <- c(spearman_cors, spearman_val)
            
#             cat("    RMSE:", rmse_val, "Pearson:", cor_val, "Spearman:", spearman_val, "\n")
            
#         }, error = function(e) {
#             cat("    Error in fold", fold_name, ":", e$message, "\n")
#             # 添加NA值
#             rmse_scores <<- c(rmse_scores, NA)
#             cors <<- c(cors, NA)
#             spearman_cors <<- c(spearman_cors, NA)
#         })
        
#         fold_counter <- fold_counter + 1
#     }
    
#     # 计算平均值
#     average_rmse <- mean(rmse_scores, na.rm = TRUE)
#     average_cor <- mean(cors, na.rm = TRUE)
#     average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
    
#     # 打印结果
#     cat("\n")
#     cat(model_name, "\n")
#     cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
#     cat("Average Mean Squared Error:", average_rmse, "\n")
#     cat("Pearson for each fold:", cors, "\n")
#     cat("Average Mean Pearson:", average_cor, "\n")
#     cat("Average Mean Spearman:", average_spearman_cor, "\n")
#     cat("\n")
    
#     # 存储结果
#     cv_scores[[model_name]] <- list(
#         rmse_scores = rmse_scores,
#         average_rmse = average_rmse,
#         cors = cors,
#         average_cor = average_cor,
#         spearman_cors = spearman_cors,
#         average_spearman_cor = average_spearman_cor
#     )
# }

library(parallel)
library(foreach)
library(doParallel)

n_cores <- detectCores() - 1  
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# only once
names(rnaseq_wgc_all_X) <- make.names(names(rnaseq_wgc_all_X), unique = TRUE)

# k-fold validation 
for (model_name in names(models)) {
    cat("Processing model:", model_name, "\n")
    
    # 预分配向量
    n_folds <- length(folds)
    rmse_scores <- numeric(n_folds)
    cors <- numeric(n_folds)
    spearman_cors <- numeric(n_folds)
    
    for (i in seq_along(folds)) {
        fold_name <- names(folds)[i]
        cat("  Fold", i, "/", n_folds, "\n")
        
        test_index <- folds[[fold_name]]
        train_index <- setdiff(1:nrow(rnaseq_wgc_all_X), test_index)
        
        # 准备数据
        X_train <- rnaseq_wgc_all_X[train_index, , drop = FALSE]
        X_test <- rnaseq_wgc_all_X[test_index, , drop = FALSE]
        y_train <- rnaseq_wgc_all_y[train_index]
        y_test <- rnaseq_wgc_all_y[test_index]
        
        # 训练和预测
        tryCatch({
            # 根据模型类型训练
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
            
            # 存储结果
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
    
    # 计算平均值
    average_rmse <- mean(rmse_scores, na.rm = TRUE)
    average_cor <- mean(cors, na.rm = TRUE)
    average_spearman_cor <- mean(spearman_cors, na.rm = TRUE)
    
    # 打印结果
    cat("\n")
    cat(model_name, "\n")
    cat("Root Mean Squared Error for each fold:", rmse_scores, "\n")
    cat("Average RMSE:", average_rmse, "\n")
    cat("Pearson for each fold:", cors, "\n")
    cat("Average Pearson:", average_cor, "\n")
    cat("Average Spearman:", average_spearman_cor, "\n\n")
    
    # 存储结果
    cv_scores[[model_name]] <- list(
        rmse_scores = rmse_scores,
        average_rmse = average_rmse,
        cors = cors,
        average_cor = average_cor,
        spearman_cors = spearman_cors,
        average_spearman_cor = average_spearman_cor
    )
}

stopCluster(cl)

# print
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

write.csv(results_summary, file = paste0("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/mine/cv_summary_", gene_select_name, "_", frag_type, ".csv"), row.names = FALSE)


# kmeans cluster
model_design <- "rnaseq_vs_hiplex_rm_outlier_log"
#wgc_raw <- data[[gene_select_name]][[model_design]][[frag_type]][["rnaseq_wgc_all_X"]]
wgc_raw <- rnaseq_wgc_all_X
rownames(wgc_raw) <- wgc_raw$gene_id
wgc_raw$gene_id <- NULL
wcss <- c()
set.seed(42)
for (n_clusters in 2:10) {
  kmeans_result <- kmeans(wgc_raw, centers = n_clusters, iter.max = 300, nstart = 10)
  wcss <- c(wcss, kmeans_result$tot.withinss)
}
wcss

find_knee <- function(x, y) {
    x_norm <- (x - min(x)) / (max(x) - min(x))
    y_norm <- (y - min(y)) / (max(y) - min(y))
    
    n <- length(x)
    x1 <- x_norm[1]; y1 <- y_norm[1]
    x2 <- x_norm[n]; y2 <- y_norm[n]
    
    distances <- c()
    for (i in 1:n) {
        d <- abs((y2 - y1) * x_norm[i] - (x2 - x1) * y_norm[i] + x2 * y1 - y2 * x1) /
            sqrt((y2 - y1)^2 + (x2 - x1)^2)
        distances <- c(distances, d)
    }
    
    return(x[which.max(distances)])
}

knee_point <- find_knee(2:10, wcss)
print(knee_point)
best_n_clusters <- knee_point

library(ggplot2)

# data
elbow_data <- data.frame(
  n_clusters = 2:10,
  wcss = wcss
)

# ggplot
p <- ggplot(elbow_data, aes(x = n_clusters, y = wcss)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = knee_point, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = knee_point, y = max(wcss) * 0.9, 
           label = paste("Optimal k =", knee_point), 
           color = "red", hjust = -0.1) +
  labs(
    title = "Elbow Method",
    x = "Number of Clusters",
    y = "Within-Cluster Sum of Squares (WCSS)"
  ) +
  scale_x_continuous(breaks = 2:10) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

save_dir <- "/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/mine"

# save picture
ggsave(
  filename = paste0(save_dir, "/kmeans_elbow_method.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


cat("Elbow plot saved to:", paste0(save_dir, "/kmeans_elbow_method.png"), "\n")

set.seed(42)
cluster_method <- kmeans(wgc_raw, centers = best_n_clusters, iter.max = 300, nstart = 10)
wgc_raw$cluster <- cluster_method$cluster



# reorder cluster label based on average gene expression
rnaseq_dir_filename = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/RNA_seq_TPM_all.csv"
rnaseq <- read.csv(rnaseq_dir_filename, sep = ",", header = TRUE, row.names = 1)
rnaseq$sqrt_V <- rnaseq$sqrt_V1V2
wgc_raw$sqrt_V <- rnaseq[rownames(wgc_raw), "sqrt_V"]

cluster_means <- aggregate(sqrt_V ~ cluster, data = wgc_raw, FUN = mean)
cluster_means <- cluster_means[order(cluster_means$sqrt_V), ]
gene_exp_order <- cluster_means$cluster
swap_dict <- setNames(0:(length(gene_exp_order) - 1), gene_exp_order)
print(cluster_means)
print(swap_dict)
wgc_raw$cluster <- swap_dict[as.character(wgc_raw$cluster)]

dir.create(paste0(save_dir, "/", model_design, "/fig/coding_all/"), recursive = TRUE, showWarnings = FALSE)
write.csv(wgc_raw, paste0(save_dir, "/", model_design, "/fig/coding_all/cluster.csv"), row.names = FALSE)

wgc_raw <- read.csv(paste0(save_dir, "/", model_design, "/fig/coding_all/cluster.csv"))

# unique cluster iD and arrange
cluster_ids <- unique(wgc_raw$cluster)
cluster_ids <- sort(cluster_ids)

# read feature colors file
feature_colors_df <- read.csv("/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_vs_hiplex_rm_outlier_log/fig/feature_colors.csv", row.names = 1)
print(feature_colors_df)

result_dir <- paste0("/dcs05/hongkai/data/next_cutntag/bulk/explainability/", model_design, "/")

library(dplyr)
library(purrr)

random_seeds <- c(0, 1, 7, 42, 123, 999, 1234, 1337, 2021, 31415)
feature_importance_list <- list()

for (i in seq_along(random_seeds)) {
  random_seed <- random_seeds[i]
  feature_importance <- read.csv(paste0(result_dir, "/", random_seed, "/top_features.csv"))
  colnames(feature_importance)[colnames(feature_importance) == "importance"] <- paste0(random_seed, "_importance")
  feature_importance_list[[i]] <- feature_importance
}

# Merge all dataframes by "feature"
feature_importance_all <- reduce(feature_importance_list, function(left, right) merge(left, right, by = "feature"))

# Set feature as rownames
rownames(feature_importance_all) <- feature_importance_all$feature
feature_importance_all$feature <- NULL

# Calculate mean importance across all seeds
feature_importance_all$importance <- rowMeans(feature_importance_all)

# Sort by importance (descending) and get top 20
feature_importance <- feature_importance_all[order(-feature_importance_all$importance), ]
feature_importance <- data.frame(
  feature = rownames(feature_importance),
  importance = feature_importance$importance
)
feature_importance <- feature_importance[1:20, ]

top_features <- feature_importance$feature
short_target_pairs <- map_target_names_pdp(feature_importance$feature, target_pair_mapping_df)
feature_importance$feature <- short_target_pairs

print(feature_importance)


feature_importance <- read.csv("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step2_output/42/top_features.csv")
feature_importance <- head(feature_importance, 20)


# output 1

# feture importance plotting
library(ggplot2)
colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#6a3d9a",
  "#b15928", "#01665e", "#8c510a", "#2c3e50", "#4b0082",
  "#003f5c", "#2f4f4f", "#800000", "#191970", "#3c1053"
)

# color 
feature_importance$colors <- colors[1:nrow(feature_importance)]
feature_importance$norm_feature_importance <- sqrt(feature_importance$importance)

# Create the plot
p <- ggplot(feature_importance, aes(x = norm_feature_importance, y = reorder(feature, norm_feature_importance), fill = feature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = setNames(feature_importance$colors, feature_importance$feature)) +
  labs(
    title = "",
    x = "Feature Importance",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x^2))
# Save the plot

out_dir <- "/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/mine"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "feature_importance_scaled.pdf"), plot = p, width = 6, height = 6, device = cairo_pdf)



# output 2 pdp

# /dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_vs_hiplex/0/top_features_pdp.pkl
library(reticulate)

# Create a Python function to load pickle files
py_run_string("
import pickle

def load_pickle_file(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)
")

result_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_vs_hiplex_rm_outlier_log"
random_seeds <- c(0, 1, 7, 42, 123, 999, 1234, 1337, 2021, 31415)
pdp_results <- list()

for (random_seed in random_seeds) {
    pdp_result_dir <- paste0(result_dir, "/", random_seed, "/top_features_pdp.pkl")
    
    # Use the Python function
    data <- py$load_pickle_file(pdp_result_dir)
    pdp_lines_result <- py_to_r(data)
    pdp_results[[as.character(random_seed)]] <- pdp_lines_result
}

cluster_results <- list()
for (random_seed in random_seeds) {
    rnaseq_wgc_all_X <- read.csv("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step2_output/42/rnaseq_wgc_all_X.csv", row.names = 1)
    rnaseq_wgc_all_X$cluster <- wgc_raw[rownames(rnaseq_wgc_all_X), "cluster"]
    rnaseq_wgc_all_X$category <- "coding_cpg"
    rnaseq_wgc_all_X[rownames(rnaseq_wgc_all_X) %in% gene_select_dict[["coding_non_cpg"]], "category"] <- "coding_non_cpg"
    cluster_results[[length(cluster_results) + 1]] <- rnaseq_wgc_all_X
}
cluster_results_stacked <- do.call(rbind, cluster_results) # same as kuai

rownames(cluster_results_stacked) <- NULL  # equivalent to ignore_index=True


## draw pdp plot for each epitope pair in each cluster

effect_size_dict <- list()
gene_select_name <- "coding_all"

# Create save directory and open PDF device FIRST
save_dir <- paste0(out_dir, "/fig/pdp_grid/")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)


# Prepare the figure with subplots
n_rows <- length(top_features)  # one row per rm_target_pair
n_cols <- 2  # two columns: centered and not_centered

# Open PDF device BEFORE plotting
pdf(file.path(save_dir, paste0(gene_select_name, "_pdp_grid.pdf")), width = 14, height = 5 * n_rows)

par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

for (col_idx in 1:2) {
    plt_type <- c("centered", "not_centered")[col_idx]
    for (row_idx in 1:length(top_features)) {
        rm_target_pair <- top_features[row_idx]
        
        print(paste0("Plotting: gene_select_name=", gene_select_name, ", plt_type=", plt_type, ", rm_target_pair=", rm_target_pair))
        
        grid_values <- pdp_results[[as.character(random_seed)]][["coding_all"]][[rm_target_pair]][["grid_values"]][[1]]
        x <- grid_values
        
        colors <- c("#D76532", "#F6C546", "#8574A4", "#5583C2", "#8CB463")
        
        # 预计算所有cluster的y值来确定y轴范围
        all_y_values <- c()
        cluster_data <- list()
        
        for (i in 1:length(cluster_ids)) {
            cluster_id <- cluster_ids[i]
            
            individual_lines <- get_average_lines(pdp_results, rm_target_pair) # same as kuai
            
            cluster_mask <- cluster_results_stacked$cluster == cluster_id
            cluster_lines <- individual_lines[cluster_mask, , drop = FALSE]
            
            if (nrow(cluster_lines) > 0) {
                y <- colMeans(cluster_lines)
                
                if (plt_type == "centered") {
                    if (min(grid_values) < 1) {
                        y_at_zero <- y[which.min(abs(x))]
                        y <- y - y_at_zero
                    } else {
                        y <- y - y[1]
                    }
                }
                
                y <- y[1:length(x)]
                cluster_data[[i]] <- y
                all_y_values <- c(all_y_values, y)
            }
        }
        
        # 计算y轴范围
        if (length(all_y_values) > 0) {
            y_range <- range(all_y_values, na.rm = TRUE)
            y_margin <- diff(y_range) * 0.1
            ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)
        } else {
            ylim <- c(-2, 2)
        }
        
        # 创建图形
        plot(NULL, xlim = c(-0.1, 1.5), ylim = ylim,
             xlab = "Hi-Plex Signal", ylab = "Gene Expression",
             cex.lab = 1.2, cex.axis = 1.1)
        
        # 绘制所有cluster的线
        for (i in 1:length(cluster_ids)) {
            if (!is.null(cluster_data[[i]])) {
                cluster_id <- cluster_ids[i]
                line_color <- colors[i]
                y <- cluster_data[[i]]
                
                lines(x, y, lwd = 1.8, col = line_color)
                text(x[length(x)] + 0.0125, y[length(y)], 
                     paste0("cluster ", cluster_id + 1),
                     col = line_color, cex = 1.3)
            }
        }
        
        rm_target_pair_short <- map_target_names_pdp(rm_target_pair, target_pair_mapping_df)[1]
        title(paste0(rm_target_pair_short, "\n(", plt_type, ")"), cex.main = 1.5)
        grid(FALSE)
    }
}

dev.off()

# output 3: ## draw pdp plot for each epitope pair in each cluster with selected epitope pairs


















# pdp function


plot_pdp_core <- function(pdp_results, rm_target_pair, cluster_results_stacked, 
                         cluster_id = NULL, plt_type = "centered", random_seed = 42) {
    # 1. 获取grid_values
    grid_values <- pdp_results[[as.character(random_seed)]][["coding_all"]][[rm_target_pair]][["grid_values"]][[1]]
    x <- grid_values
    
    # 2. 获取individual_lines
    individual_lines <- get_average_lines(pdp_results, rm_target_pair)
    
    # 3. 根据cluster_id筛选数据
    if (!is.null(cluster_id)) {
        cluster_mask <- cluster_results_stacked$cluster == cluster_id
        individual_lines <- individual_lines[cluster_mask, , drop = FALSE]
    }
    
    # 4. 计算平均值
    y <- colMeans(individual_lines)
    
    # 5. 处理centered逻辑
    if (plt_type == "centered") {
        if (min(grid_values) < 1) {
            y_at_zero <- y[which.min(abs(x))]
            y <- y - y_at_zero
        } else {
            y <- y - y[1]
        }
    }
    
    # 6. 截取到x的长度
    y <- y[1:length(x)]
    
    return(list(x = x, y = y))
}



# output 2: pair
# Figure 1
n_rows <- length(top_features)
n_cols <- 2

pdf(file.path(save_dir, paste0(gene_select_name, "_pdp_grid.pdf")), width = 14, height = 5 * n_rows)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

for (col_idx in 1:2) {
    plt_type <- c("centered", "not_centered")[col_idx]
    for (row_idx in 1:length(top_features)) {
        rm_target_pair <- top_features[row_idx]
        
        print(paste0("Plotting: gene_select_name=", gene_select_name, ", plt_type=", plt_type, ", rm_target_pair=", rm_target_pair))
        
        colors <- c("#D76532", "#F6C546", "#8574A4", "#5583C2", "#8CB463")
        
        all_y_values <- c()
        cluster_data <- list()
        
        for (i in 1:length(cluster_ids)) {
            cluster_id <- cluster_ids[i]
            
            # 调用核心函数
            result <- plot_pdp_core(pdp_results, rm_target_pair, cluster_results_stacked, 
                                   cluster_id = cluster_id, plt_type = plt_type, random_seed = random_seed)
            
            if (length(result$y) > 0) {
                cluster_data[[i]] <- result$y
                all_y_values <- c(all_y_values, result$y)
            }
        }
        
        # 计算y轴范围
        if (length(all_y_values) > 0) {
            y_range <- range(all_y_values, na.rm = TRUE)
            y_margin <- diff(y_range) * 0.1
            ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)
        } else {
            ylim <- c(-2, 2)
        }
        
        # 创建图形
        plot(NULL, xlim = c(-0.1, 1.5), ylim = ylim,
             xlab = "Hi-Plex Signal", ylab = "Gene Expression",
             cex.lab = 1.2, cex.axis = 1.1)
        
        for (i in 1:length(cluster_ids)) {
            if (!is.null(cluster_data[[i]])) {
                cluster_id <- cluster_ids[i]
                line_color <- colors[i]
                
                result <- plot_pdp_core(pdp_results, rm_target_pair, cluster_results_stacked, 
                                       cluster_id = cluster_id, plt_type = plt_type, random_seed = random_seed)
                x <- result$x
                y <- result$y
                
                lines(x, y, lwd = 1.8, col = line_color)
                text(x[length(x)] + 0.0125, y[length(y)], 
                     paste0("cluster ", cluster_id + 1),
                     col = line_color, cex = 1.3)
            }
        }
        
        rm_target_pair_short <- map_target_names_pdp(rm_target_pair, target_pair_mapping_df)[1]
        title(paste0(rm_target_pair_short, "\n(", plt_type, ")"), cex.main = 1.5)
        grid(FALSE)
    }
}

dev.off()


# output 3: all
# Figure 2
## chosen epitope pairs，each cluster a plot

save_dir <- paste0("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/", model_design, "/fig/coding_all/sqrt_coding_all_centered_same_range")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)


for (i in 1:length(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    
    file_path <- file.path(save_dir, paste0(cluster_id, ".pdf"))
    cat("Saving to:", file_path, "\n")
    
    pdf(file_path, width = 6, height = 6)
    
    colors <- c("#D76532", "#F6C546", "#8574A4", "#4C9F70")
    target_pairs <- c("H3K27me3-H3K4me3", "H3K4me1-H3K4me3", "H3K27ac-H3K4me3", "H3K27me3-H3K9me3")
    
    plot(NULL, xlim = c(-0.05, 1.5), ylim = c(-0.5, 0.75),
         xlab = "Hi-Plex Signal", ylab = "Gene Expression")
    
    for (row_idx in 1:length(target_pairs)) {
        rm_target_pair <- target_pairs[row_idx]
        
        result <- plot_pdp_core(pdp_results, rm_target_pair, cluster_results_stacked, 
                               cluster_id = cluster_id, plt_type = "centered", random_seed = random_seed)
        
        x <- result$x
        y <- result$y
        y <- sign(y) * sqrt(abs(y))  
        
        line_color <- colors[row_idx]
        lines(x, y, lwd = 1.8, col = line_color)
        text(x[length(x)] - 0.5, y[length(y)] + 0.025, rm_target_pair,
             col = line_color, cex = 1.3)
    }
    
    title(paste0("Cluster ", cluster_id + 1), cex.main = 1.5)
    dev.off()

}

print(list.files(save_dir))




# output 4
# Figure 3
## change label
# Figure 3  - 
save_dir <- paste0("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/", model_design, "/fig/", gene_select_name, "/")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

width <- 10
height <- 10
pdf(file.path(save_dir, paste0("summary_", width, "x", height, ".pdf")), width = width, height = height)

colors <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#6a3d9a",
    "#b15928", "#01665e", "#8c510a", "#2c3e50", "#4b0082",
    "#003f5c", "#2f4f4f", "#800000", "#191970", "#3c1053"
)


all_x_values <- c()
all_y_values <- c()
plot_data <- list()
text_positions <- list()

for (i in 1:length(top_features)) {
    target_pair <- top_features[i]
    
    result <- plot_pdp_core(pdp_results, target_pair, cluster_results_stacked, 
                           cluster_id = NULL, plt_type = "centered", random_seed = random_seed)
    
    x <- result$x
    y <- result$y
    
    all_x_values <- c(all_x_values, x)
    all_y_values <- c(all_y_values, y)
    
    plot_data[[i]] <- list(x = x, y = y, target_pair = target_pair, color = colors[i])
    
    text_positions[[i]] <- list(
        x = x[length(x)],
        y = y[length(y)],
        label = target_pair,
        color = colors[i]
    )
}

x_range <- range(all_x_values, na.rm = TRUE)
y_range <- range(all_y_values, na.rm = TRUE)
x_margin <- diff(x_range) * 0.25  
y_margin <- diff(y_range) * 0.2
xlim <- c(x_range[1] - x_margin * 0.2, x_range[2] + x_margin)
ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)

plot(NULL, xlim = xlim, ylim = ylim,
     xlab = "Hi-Plex Signal", ylab = "Gene Expression",
     cex.lab = 1.2, cex.axis = 1.1)

for (i in 1:length(plot_data)) {
    data <- plot_data[[i]]
    lines(data$x, data$y, lwd = 2, col = data$color, lty = "solid")
}

text_y_positions <- sapply(text_positions, function(tp) tp$y)
text_x_positions <- sapply(text_positions, function(tp) tp$x)

sorted_indices <- order(text_y_positions)

n_labels <- length(text_positions)
y_spread <- diff(ylim) * 0.8  
y_center <- mean(ylim)
y_start <- y_center - y_spread/2
y_step <- y_spread / (n_labels - 1)

for (i in 1:length(text_positions)) {
    idx <- sorted_indices[i]
    tp <- text_positions[[idx]]
    
    new_y <- y_start + (i - 1) * y_step
    new_x <- tp$x + 0.05  
    
    segments(tp$x, tp$y, new_x - 0.02, new_y, 
             col = tp$color, lty = 2, lwd = 0.8, alpha = 0.7)
    
    text(new_x, new_y, tp$label,
         col = tp$color, cex = 0.9, adj = 0)
}

title("Coding Genes", cex.main = 1.8)
par(bty = "l")
dev.off()


## original
# # all epitope 
# save_dir <- paste0("/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/Step3_output/", model_design, "/fig/", gene_select_name, "/")
# dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# width <- 10
# height <- 10

# pdf(file.path(save_dir, paste0("summary_", width, "x", height, ".pdf")), width = width, height = height)

# colors <- c(
#     "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
#     "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#6a3d9a",
#     "#b15928", "#01665e", "#8c510a", "#2c3e50", "#4b0082",
#     "#003f5c", "#2f4f4f", "#800000", "#191970", "#3c1053"
# )

# all_x_values <- c()
# all_y_values <- c()
# plot_data <- list()

# for (i in 1:length(top_features)) {
#     target_pair <- top_features[i]
    
#     result <- plot_pdp_core(pdp_results, target_pair, cluster_results_stacked, 
#                            cluster_id = NULL, plt_type = "centered", random_seed = random_seed)
    
#     x <- result$x
#     y <- result$y
    
#     all_x_values <- c(all_x_values, x)
#     all_y_values <- c(all_y_values, y)
    
#     plot_data[[i]] <- list(x = x, y = y, target_pair = target_pair, color = colors[i])
# }

# x_range <- range(all_x_values, na.rm = TRUE)
# y_range <- range(all_y_values, na.rm = TRUE)

# x_margin <- diff(x_range) * 0.1  
# y_margin <- diff(y_range) * 0.1  

# xlim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
# ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)


# plot(NULL, xlim = xlim, ylim = ylim,
#      xlab = "Hi-Plex Signal", ylab = "Gene Expression",
#      cex.lab = 1.2, cex.axis = 1.1)

# for (i in 1:length(plot_data)) {
#     data <- plot_data[[i]]
    
#     lines(data$x, data$y, lwd = 2, col = data$color, lty = "solid")
    
#     text(data$x[length(data$x)] + 0.0125, data$y[length(data$y)], data$target_pair,
#          col = data$color, cex = 1.3)
# }

# title("Coding Genes", cex.main = 1.8)

# par(bty = "l")

# dev.off()






# output 5: csv file
# effect size calculation
top_features <- as.list(top_features)
gene_select_name <- 'coding_all'

for (effect_size_method in c("auc")) {
    # Iterate over plot types and features (rm_target_pair)
    for (col_idx in 1:1) {  
        plt_type <- "centered"
        
        # create effect_size_df grid
        effect_size_df <- matrix(0, nrow = length(cluster_ids), ncol = length(top_features))
        rownames(effect_size_df) <- cluster_ids
        colnames(effect_size_df) <- top_features
        effect_size_df <- as.data.frame(effect_size_df)
        
        for (row_idx in 1:length(top_features)) {
            rm_target_pair <- top_features[[row_idx]]
            
            grid_values <- pdp_results[[as.character(random_seed)]][["coding_all"]][[rm_target_pair]][["grid_values"]][[1]]
            x <- grid_values
            
            for (i in 1:length(cluster_ids)) {
                cluster_id <- cluster_ids[i]
                
                individual_lines <- get_average_lines(pdp_results, rm_target_pair)
                
                if (gene_select_name == "coding_all") {
                    cluster_mask <- cluster_results_stacked$cluster == cluster_id
                    individual_lines <- individual_lines[cluster_mask, , drop = FALSE]
                } else {
                    cluster_mask <- (cluster_results_stacked$cluster == cluster_id) & 
                                   (cluster_results_stacked$category == gene_select_name)
                    individual_lines <- individual_lines[cluster_mask, , drop = FALSE]
                }
                
                y <- colMeans(individual_lines)
                
                if (plt_type == "centered") {
                    if (min(grid_values) < 1) {
                        y_at_zero <- y[which.min(abs(x))]
                        y <- y - y_at_zero
                    } else {
                        y <- y - y[1]
                    }
                }
                
                effect_size <- calculate_effect_size(x, y, method = effect_size_method)
                effect_size_df[as.character(cluster_id), rm_target_pair] <- effect_size
            }
        }
        
        effect_size_df <- effect_size_df / apply(effect_size_df, 2, sd)
        
        # sort index
        effect_size_df <- effect_size_df[order(as.numeric(rownames(effect_size_df))), ]
        
        # col name
        short_target_pairs <- map_target_names_pdp(colnames(effect_size_df), target_pair_mapping_df)
        colnames(effect_size_df) <- short_target_pairs
        
        # index +1
        rownames(effect_size_df) <- as.numeric(rownames(effect_size_df)) + 1
        
        # save CSV
        write.csv(effect_size_df, file.path(save_dir, paste0('effect_size_', effect_size_method, '.csv')), row.names = TRUE)
    }
}