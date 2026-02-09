#!/usr/bin/env Rscript

library(arrow)
library(jsonlite)
library(dplyr)
library(caret)
library(randomForest)

# Function to perform model grid search
# perform_model_gridsearch <- function(gene_select_path, gene_select_name, rna_path, model_design, hiplex_data_path, filtered_target_pair_path, cutoff_path) {
#   suppressPackageStartupMessages({
#     library(arrow)
#     library(jsonlite)
#     library(dplyr)
#     library(caret)
#     library(randomForest)
#   }) 

#   # Load gene selection dictionary
#   gene_select_dict <- fromJSON(gene_select_path)
#   gene_select_dict$coding_all <- c(gene_select_dict$coding_cpg, gene_select_dict$coding_non_cpg)
#   gene_select_list <- gene_select_dict[[gene_select_name]]
  
#   # Load RNA-seq data
#   rnaseq <- read.csv(rna_path, row.names = 1)
#   rnaseq$sqrt_V <- log10(rnaseq$V1V2 + 1)
  
#   # Load HiPlex data
#   wgc_raw <- read_feather(hiplex_data_path)
  
#   features <- fromJSON(filtered_target_pair_path)
#   wgc_vals <- wgc_raw[, features]
#   zero_hiplex_genes <- wgc_raw$pos[rowSums(wgc_vals) == 0]
  
#   # Log transform and normalize HiPlex features
#   wgc_raw[, features] <- log10(wgc_raw[, features] + 1)

#   global_min <- min(wgc_raw[, features])
#   global_max <- max(wgc_raw[, features])
#   wgc_raw[, features] <- (wgc_raw[, features] - global_min) / (global_max - global_min)  

#   # Filter and merge data
#   rnaseq_cutoffs <- fromJSON(cutoff_path)
#   q99 <- rnaseq_cutoffs[[model_design]]
  
#   overlap_genes <- intersect(wgc_raw$pos, gene_select_list)
#   rnaseq_avail <- rnaseq[overlap_genes, c("gene_id", "sqrt_V")]
#   rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, ]
  
#   zero_rnaseq_genes <- rownames(rnaseq_avail)[rnaseq_avail$sqrt_V == 0]
#   zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)
  
#   # Merge datasets
#   rnaseq_wgc <- merge(rnaseq_avail, wgc_raw, by.x = "row.names", by.y = "pos")
#   rnaseq_wgc <- rnaseq_wgc[!rnaseq_wgc$Row.names %in% zero_all_genes, ]
#   rnaseq_wgc <- rnaseq_wgc %>% select(-Row.names, -gene_id)
  
#   # Train-test split
#   set.seed(42)
#   train_idx <- createDataPartition(rnaseq_wgc$sqrt_V, p = 0.8, list = FALSE)
#   train_data <- rnaseq_wgc[train_idx, ]
#   test_data <- rnaseq_wgc[-train_idx, ]
  
#   # Prepare for grid search
#   n_features <- ncol(train_data) - 1
#   train_control <- trainControl(method = "cv", number = 5, verboseIter = TRUE)
  
#   # Create comprehensive parameter grid for Random Forest
#   rf_param_grid <- expand.grid(
#     mtry = c(log2(n_features), sqrt(n_features), n_features),  # max_features: 'log2', 'sqrt', 1.0
#     ntree = c(50, 100, 150),  # n_estimators
#     nodesize = c(1, 2, 4),  # min_samples_leaf
#     maxnodes = c(10, 20, 30, -1),  # max_depth (use -1 for None/unlimited)
#     replace = c(TRUE, FALSE)  # bootstrap
#   )
  
#   # Additional parameter for min_samples_split (handled via sampsize in RF)
#   min_samples_split <- c(2, 5, 10)
  
#   result <- list()
  
#   # ============ Random Forest Grid Search ============
#   print("start training model rf")
#   flush.console()
  
#   rf_results <- list()
#   best_rf_rmse <- Inf
#   best_rf_params <- NULL
#   best_rf_model <- NULL
  
#   for (i in 1:nrow(rf_param_grid)) {
#     for (split_val in min_samples_split) {
#       params <- rf_param_grid[i, ]
      
#       # Set maxnodes (NULL if -1, otherwise use the value)
#       maxnodes_val <- if(params$maxnodes == -1) NULL else params$maxnodes
      
#       tryCatch({
#         set.seed(42)
#         rf_temp <- train(
#           sqrt_V ~ ., 
#           data = train_data, 
#           method = "rf",
#           trControl = train_control,
#           tuneGrid = data.frame(mtry = round(params$mtry)),
#           ntree = params$ntree,
#           nodesize = params$nodesize,
#           maxnodes = maxnodes_val,
#           replace = params$replace,
#           metric = "RMSE"
#         )
        
#         current_rmse <- min(rf_temp$results$RMSE)
        
#         if (current_rmse < best_rf_rmse) {
#           best_rf_rmse <- current_rmse
#           best_rf_params <- list(
#             mtry = round(params$mtry),
#             ntree = params$ntree,
#             nodesize = params$nodesize,
#             maxnodes = if(params$maxnodes == -1) "None" else params$maxnodes,
#             replace = params$replace,
#             min_samples_split = split_val
#           )
#           best_rf_model <- rf_temp
#         }
        
#         cat(sprintf("RF - Params: mtry=%d, ntree=%d, nodesize=%d, maxnodes=%s, replace=%s, split=%d | RMSE=%.4f\n",
#                     round(params$mtry), params$ntree, params$nodesize, 
#                     ifelse(params$maxnodes == -1, "None", as.character(params$maxnodes)),
#                     params$replace, split_val, current_rmse))
#         flush.console()
        
#       }, error = function(e) {
#         cat(sprintf("Error with params: %s\n", paste(params, collapse=", ")))
#       })
#     }
#   }
  
#   result[[1]] <- list(
#     model_name = "rf",
#     best_score = -best_rf_rmse,  # Negative to match Python's neg_root_mean_squared_error
#     best_params = best_rf_params,
#     model = best_rf_model
#   )
  
#   # ============ Linear Regression ============
#   print("start training model lr")
#   flush.console()
  
#   set.seed(42)
#   lr_model <- train(
#     sqrt_V ~ ., 
#     data = train_data, 
#     method = "lm",
#     trControl = train_control,
#     metric = "RMSE"
#   )
  
#   lr_rmse <- min(lr_model$results$RMSE)
  
#   result[[2]] <- list(
#     model_name = "lr",
#     best_score = -lr_rmse,  # Negative to match Python
#     best_params = list(),
#     model = lr_model
#   )
  
#   # Sort results by best score (descending, since scores are negative)
#   result <- result[order(sapply(result, function(x) x$best_score), decreasing = TRUE)]
  
#   # Print summary
#   cat("\n========== Grid Search Results ==========\n")
#   for(i in 1:length(result)) {
#     cat(sprintf("\nRank %d - Model: %s\n", i, result[[i]]$model_name))
#     cat(sprintf("Best Score (negative RMSE): %.6f\n", result[[i]]$best_score))
#     cat("Best Params:\n")
#     print(result[[i]]$best_params)
#   }
  
#   # Save results (without the model objects to keep file size small)
#   results_to_save <- lapply(result, function(x) {
#     list(
#       model_name = x$model_name,
#       best_score = x$best_score,
#       best_params = x$best_params
#     )
#   })
  
#   dir.create(paste0('/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/', model_design), 
#              showWarnings = FALSE, recursive = TRUE)
#   write_json(results_to_save, 
#              paste0('/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/', 
#                     model_design, '/', gene_select_name, '.json'),
#              pretty = TRUE, auto_unbox = TRUE)
  
#   return(result)
# }

perform_model_gridsearch <- function(gene_select_path, gene_select_name, rna_path, model_design, hiplex_data_path, filtered_target_pair_path, cutoff_path, GE_col_name, out_dir) {
  suppressPackageStartupMessages({
    library(arrow)
    library(jsonlite)
    library(dplyr)
    library(caret)
    library(randomForest)
    library(BiocParallel)
  }) 

  # Load gene selection dictionary
  gene_select_dict <- fromJSON(gene_select_path) # list with 5
  gene_select_dict$coding_all <- c(gene_select_dict$coding_cpg, gene_select_dict$coding_non_cpg)
  gene_select_list <- gene_select_dict[[gene_select_name]]
  
  # Load RNA-seq data
  rnaseq <- read.csv(rna_path, row.names = 1)
  rnaseq$sqrt_V <- log10(rnaseq[[GE_col_name]] + 1)
  
  # Load HiPlex data
  wgc_raw <- read_feather(hiplex_data_path)
  features <- fromJSON(filtered_target_pair_path)
  wgc_vals <- wgc_raw[, features]
  zero_hiplex_genes <- wgc_raw$pos[rowSums(wgc_vals) == 0]
  
  # Log transform and normalize HiPlex features
  wgc_raw[, features] <- log10(wgc_raw[, features] + 1)
  global_min <- min(wgc_raw[, features])
  global_max <- max(wgc_raw[, features])
  wgc_raw[, features] <- (wgc_raw[, features] - global_min) / (global_max - global_min)  

  # Filter and merge data
  rnaseq_cutoffs <- fromJSON(cutoff_path)
  q99 <- rnaseq_cutoffs[[model_design]]
  
  overlap_genes <- intersect(wgc_raw$pos, gene_select_list)
  rnaseq_avail <- rnaseq[overlap_genes, c("gene_id", "sqrt_V")]
  rnaseq_avail <- rnaseq_avail[rnaseq_avail$sqrt_V <= q99, ]
  
  zero_rnaseq_genes <- rownames(rnaseq_avail)[rnaseq_avail$sqrt_V == 0]
  zero_all_genes <- intersect(zero_hiplex_genes, zero_rnaseq_genes)
  
  # Merge datasets
  rnaseq_wgc <- merge(rnaseq_avail, wgc_raw, by.x = "row.names", by.y = "pos")
  rnaseq_wgc <- rnaseq_wgc[!rnaseq_wgc$Row.names %in% zero_all_genes, ]
  rnaseq_wgc <- rnaseq_wgc %>% select(-Row.names, -gene_id)
  
  # Train-test split
  set.seed(42)
  train_idx <- createDataPartition(rnaseq_wgc$sqrt_V, p = 0.8, list = FALSE)
  train_data <- rnaseq_wgc[train_idx, ]
  test_data <- rnaseq_wgc[-train_idx, ]
  
  # Prepare for grid search
  n_features <- ncol(train_data) - 1
  train_control <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
  
  # Create comprehensive parameter grid for Random Forest
  rf_param_grid <- expand.grid(
    mtry = c(log2(n_features), sqrt(n_features), n_features),
    ntree = c(50, 100, 150),
    nodesize = c(1, 2, 4),
    maxnodes = c(10, 20, 30, -1),
    replace = c(TRUE, FALSE),
    min_samples_split = c(2, 5, 10)
  )
  
  result <- list()
  
  # ============ Random Forest Grid Search with BiocParallel ============
  # Setup BiocParallel backend
  num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
  bp_param <- MulticoreParam(workers = num_cores, progressbar = TRUE)
  register(bp_param)
  print("start training model rf")
  cat(sprintf("Using %d cores for parallel processing\n", num_cores))
  flush.console()

  # Function to train one parameter set
  train_one_rf <- function(i, train_data, train_control, rf_param_grid) {
    library(caret)
    library(randomForest)
    
    params <- rf_param_grid[i, ]
    maxnodes_val <- if(params$maxnodes == -1) NULL else params$maxnodes
    
    tryCatch({
      set.seed(42)
      rf_temp <- train(
        sqrt_V ~ ., 
        data = train_data, 
        method = "rf",
        trControl = train_control,
        tuneGrid = data.frame(mtry = round(params$mtry)),
        ntree = params$ntree,
        nodesize = params$nodesize,
        maxnodes = maxnodes_val,
        replace = params$replace,
        metric = "RMSE"
      )
      
      current_rmse <- min(rf_temp$results$RMSE)
      
      list(
        mtry = round(params$mtry),
        ntree = params$ntree,
        nodesize = params$nodesize,
        maxnodes = ifelse(params$maxnodes == -1, -1, params$maxnodes),
        replace = params$replace,
        min_samples_split = params$min_samples_split,
        rmse = current_rmse
      )
    }, error = function(e) {
      list(
        mtry = round(params$mtry),
        ntree = params$ntree,
        nodesize = params$nodesize,
        maxnodes = ifelse(params$maxnodes == -1, -1, params$maxnodes),
        replace = params$replace,
        min_samples_split = params$min_samples_split,
        rmse = Inf
      )
    })
  }

  
  # Run parallel grid search using BiocParallel
  rf_results_list <- bplapply(1:nrow(rf_param_grid), train_one_rf, 
                               train_data = train_data,
                               train_control = train_control,
                               rf_param_grid = rf_param_grid,
                               BPPARAM = bp_param)
  
  # Convert list to dataframe
  rf_results <- do.call(rbind, lapply(rf_results_list, function(x) as.data.frame(x)))
  
  # Find best parameters
  best_idx <- which.min(rf_results$rmse)
  best_rf_params <- as.list(rf_results[best_idx, ])
  best_rf_rmse <- rf_results$rmse[best_idx]
  
  cat(sprintf("\nBest RF params found: mtry=%d, ntree=%d, nodesize=%d, maxnodes=%s, replace=%s, split=%d | RMSE=%.4f\n",
              best_rf_params$mtry, best_rf_params$ntree, best_rf_params$nodesize,
              ifelse(best_rf_params$maxnodes == -1, "None", as.character(best_rf_params$maxnodes)),
              best_rf_params$replace, best_rf_params$min_samples_split, best_rf_rmse))
  
  # Retrain best model
  cat("Retraining best RF model...\n")
  maxnodes_val <- if(best_rf_params$maxnodes == -1) NULL else best_rf_params$maxnodes
  set.seed(42)
  best_rf_model <- train(
    sqrt_V ~ ., 
    data = train_data, 
    method = "rf",
    trControl = train_control,
    tuneGrid = data.frame(mtry = best_rf_params$mtry),
    ntree = best_rf_params$ntree,
    nodesize = best_rf_params$nodesize,
    maxnodes = maxnodes_val,
    replace = best_rf_params$replace,
    metric = "RMSE"
  )
  
  result[[1]] <- list(
    model_name = "rf",
    best_score = -best_rf_rmse,
    best_params = lapply(best_rf_params[1:6], function(x) if(x == -1) "None" else x),
    model = best_rf_model
  )
  
  # ============ Linear Regression ============
  print("start training model lr")
  flush.console()
  
  set.seed(42)
  lr_model <- train(
    sqrt_V ~ ., 
    data = train_data, 
    method = "lm",
    trControl = train_control,
    metric = "RMSE"
  )
  
  lr_rmse <- min(lr_model$results$RMSE)
  
  result[[2]] <- list(
    model_name = "lr",
    best_score = -lr_rmse,
    best_params = list(),
    model = lr_model
  )
  
  # Sort results
  result <- result[order(sapply(result, function(x) x$best_score), decreasing = TRUE)]
  
  # Print summary
  cat("\n========== Grid Search Results ==========\n")
  for(i in 1:length(result)) {
    cat(sprintf("\nRank %d - Model: %s\n", i, result[[i]]$model_name))
    cat(sprintf("Best Score (negative RMSE): %.6f\n", result[[i]]$best_score))
    cat("Best Params:\n")
    print(result[[i]]$best_params)
  }
  
  # Save results
  results_to_save <- lapply(result, function(x) {
    list(
      model_name = x$model_name,
      best_score = x$best_score,
      best_params = x$best_params
    )
  })
  
  dir.create(paste0(out_dir, model_design), showWarnings = FALSE, recursive = TRUE)
  write_json(results_to_save, 
             paste0(out_dir, 
                    model_design, '/', gene_select_name, '.json'),
             pretty = TRUE, auto_unbox = TRUE)
  
  return(result)
}

# Set default paths
gene_select_path <- "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island/gene_categories.json"
rna_path <- "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/RNA_seq_TPM_all.csv"
hiplex_data_path <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/promoter_-1000-1000/V_mixed_promoter_-1000-1000_colQC-all_libnorm.feather"
filtered_target_pair_path <- "/dcs05/hongkai/data/next_cutntag/script/utils/filtered_target_pairs.json"
cutoff_path <- "/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_hiplex_cutoff.json"
gene_select_name <- "coding_all"
model_design <- "rnaseq_vs_hiplex_rm_outlier_log"
GE_col_name <- "V1V2"

# Run function
results <- perform_model_gridsearch(
  gene_select_path = gene_select_path,
  gene_select_name = gene_select_name,
  rna_path = rna_path,
  model_design = model_design,
  hiplex_data_path = hiplex_data_path,
  filtered_target_pair_path = filtered_target_pair_path,
  cutoff_path = cutoff_path,
  GE_col_name = GE_col_name,
  out_dir = "/dcs05/hongkai/data/yhu1/next_cutntag/AAA_last_part_from_kuai/"
)

print(results)