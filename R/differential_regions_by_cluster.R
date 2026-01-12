

differential_regions_by_cluster <- function(col_cluster_file_path, cm_path, conditions, sample_names = NULL, out_dir = "./", cpm_scale = 1E6, lfc_threshold = 0.5, fdr_threshold = 0.05, mean_quantile = 0.25, log_pseudocount = 1) {
  # Load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(latex2exp)
    library(edgeR)
    library(matrixStats)
    library(limma)
  })

  if (length(cm_path) != length(conditions)) {
    stop("cm_path and conditions must have the same length.")
  }

  cond_table <- table(conditions)
  if (any(cond_table < 2)) {
    stop(paste0("Each condition must have at least 2 replicates. Failed for: ", paste(names(cond_table)[cond_table < 2], collapse = ", ")))
  }

  # Reorder by conditions
  ord <- order(conditions)
  cm_path <- cm_path[ord]
  conditions <- conditions[ord]

  if (is.null(sample_names)) {
    sample_names <- paste0(conditions, ave(seq_along(conditions), conditions, FUN = seq_along))
  } else {
    if (length(sample_names) != length(conditions)) {
      stop("Length of sample_names must equal length of conditions.")
    }
    if (anyDuplicated(sample_names)) {
      stop("sample_names must be unique.")
    }
    sample_names <- sample_names[ord]
  }


  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Read .feather file
  cm_list <- lapply(cm_path, function(f) column_to_rownames(read_feather(f), var = "pos"))
  names(cm_list) <- sample_names

  # Check if all count matrices have the same row/col names
  all_regions <- lapply(cm_list, rownames)
  if (length(unique(sapply(all_regions, paste, collapse = "|"))) > 1) {
    warning("Row names (pos) do not match across all count matrices. Taking intersection.")
    common_regions <- Reduce(intersect, all_regions)
    if (length(common_regions) == 0) {
      stop("No common row names (pos) found across all count matrices.")
    }
    message(glue("Using {length(common_regions)} common rows out of {paste(sapply(all_regions, length), collapse = ', ')} rows in each file."))
    cm_list <- lapply(cm_list, function(mat) mat[common_regions, , drop = FALSE])
  }

  all_pairs <- lapply(cm_list, colnames)
  if (length(unique(sapply(all_pairs, paste, collapse = "|"))) > 1) {
    warning("Column names (pairs) do not match across all count matrices. Taking intersection.")
    common_pairs <- Reduce(intersect, all_pairs)
    if (length(common_pairs) == 0) {
      stop("No common column names found across all count matrices.")
    }
    message(glue("Using {length(common_pairs)} common columns out of {paste(sapply(all_pairs, length), collapse = ', ')} columns in each file."))
    cm_list <- lapply(cm_list, function(mat) mat[, common_pairs, drop = FALSE])
  }

  cm_regions <- rownames(cm_list[[1]])
  cm_pairs <- colnames(cm_list[[1]])

  # Read column cluster file
  col_cluster_table <- read.table(col_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  cluster_pairs <- col_cluster_table$pair
  col_cluster_table <- col_cluster_table[col_cluster_table$pair %in% cm_pairs, ]

  cm_not_in_cluster <- setdiff(cm_pairs, cluster_pairs)
  cluster_not_in_cm <- setdiff(cluster_pairs, cm_pairs)
  if (length(cm_not_in_cluster) > 0) {
    warning("Count matrices contain ", length(cm_not_in_cluster), " column(s) not found in col_cluster_file")
  }
  if (length(cluster_not_in_cm) > 0) {
    warning("col_cluster_file contains ", length(cluster_not_in_cm), " pair(s) not found in count matrices")
  }

  # Setup design matrix
  conditions <- factor(conditions)
  design <- model.matrix(~ 0 + conditions)
  colnames(design) <- levels(conditions)
  
  # Generate all pairwise comparisons
  cond_levels <- levels(conditions)
  n_cond <- length(cond_levels)
  comparisons <- list()

  if (n_cond == 2) {
    # For 2 conditions, just one comparison
    comparisons[[1]] <- c(cond_levels[2], cond_levels[1])
  } else {
    # For >2 conditions, all pairwise comparisons
    idx <- 1
    for (i in 1:(n_cond - 1)) {
      for (j in (i + 1):n_cond) {
        comparisons[[idx]] <- c(cond_levels[j], cond_levels[i])
        idx <- idx + 1
      }
    }
  }

}

# Differential Analysis
# differential -1
# two condition, 4 columns
# Post: Perform differential analysis between two conditions using limma-voom pipeline, analyzing each column cluster separately to identify significantly different genomic regions.
# Parameter: sample_names: Vector of sample names in order (first case_num samples for condition1, then control_num samples for condition2)
#            case_num: Number of replicates in condition 1 (minimum 2)
#            control_num: Number of replicates in condition 2 (minimum 2)
#            col_cluster_file: Path to column cluster assignment TSV file
#            wgc_file_path: Vector of paths to count matrix feather files for each sample
#            sig_result_dir: Directory to save differential analysis results
#            pseudocount: Pseudocount for normalization (default: 0.5)
#            normalization_factor: Scaling factor for CPM calculation (default: 1E6)
#            lowess_span: Span parameter for voom lowess fitting (default: 0.5)
#            l2fc_thres: Log2 fold change threshold for significance (default: 0.5)
#            mean_per_thres_list: Vector of mean expression percentile thresholds for filtering (default: c(0.25))
#            fdr_thres_list: Vector of FDR thresholds for significance calling (default: c(0.25))
# Output: Saves differential analysis results, filtered count matrices, and significant regions for each cluster and threshold combination

differential_regions_by_cluster <- function(sample_names, case_num, control_num, wgc_file_path, sig_result_dir, col_cluster_file = NULL, normalization_factor = 1E6, lowess_span = 0.5, l2fc_thres = 0.5, mean_per_thres_list = c(0.25), fdr_thres_list = c(0.25), pseudocount_for_log = 1) {
  
  group1 <- sample_names[1:case_num]
  group2 <- sample_names[(case_num + 1):(case_num + control_num)]
  conditions <- c(rep("condition1", case_num), rep("condition2", control_num))
  coldata <- data.frame("condition" = conditions, row.names = sample_names)
  group <- factor(coldata$condition)
  condition_levels <- levels(group)
  mm <- model.matrix(~ 0 + group)

  wgc_list <- lapply(wgc_file_path, function(f) column_to_rownames(read_feather(f), var = "pos"))
  names(wgc_list) <- sample_names

  # Check if all WGC files have matching column names
  all_colnames <- lapply(wgc_list, colnames)
  all_colnames_sorted <- lapply(all_colnames, sort)
  colnames_strings <- sapply(all_colnames_sorted, paste, collapse = "|")

  if (length(unique(colnames_strings)) > 1) {
    mismatch_info <- c()
    for (i in 1:(length(all_colnames) - 1)) {
      for (j in (i + 1):length(all_colnames)) {
        if (!identical(all_colnames_sorted[[i]], all_colnames_sorted[[j]])) {
          mismatch_info <- c(
            mismatch_info,
            paste0("  File ", i, " (", sample_names[i], ") vs File ", j, " (", sample_names[j], ")")
          )
        }
      }
    }
    stop("Column names do not match across WGC files:\n", paste(unique(mismatch_info), collapse = "\n"))
  }

  # Handle column cluster file
  if (is.null(col_cluster_file)) {
    # Get all unique column names from all WGC files
    all_features <- unique(unlist(lapply(wgc_list, colnames)))
    # Create a data frame where each feature gets its own label
    col_cluster <- data.frame(
      feature = all_features,
      label = seq_along(all_features),
      stringsAsFactors = FALSE
    )
  } else {
    # Load column cluster file
    col_cluster <- read.table(col_cluster_file, header = TRUE, sep = "\t", row.names = NULL)
  }

  # Check intersection between WGC columns and col_cluster features
  wgc_cols <- colnames(wgc_list[[1]])
  cluster_features <- col_cluster$feature

  wgc_not_in_cluster <- setdiff(wgc_cols, cluster_features)
  cluster_not_in_wgc <- setdiff(cluster_features, wgc_cols)

  if (length(wgc_not_in_cluster) > 0) {
    warning("WGC files contain ", length(wgc_not_in_cluster), " column(s) not found in col_cluster file")
  }

  if (length(cluster_not_in_wgc) > 0) {
    warning("col_cluster file contains ", length(cluster_not_in_wgc), " feature(s) not found in WGC files")
  }

  # Loop over column clusters
  col_cluster <- col_cluster[col_cluster$feature %in% wgc_cols, ]
  col_label_list <- unique(col_cluster$label)
  for (col_label in col_label_list) {
    target_pair_select_list <- col_cluster[col_cluster$label == col_label, "feature"]

    # Combine counts for this cluster
    tmp_combine_orig <- do.call(
      cbind,
      lapply(sample_names, function(s) rowSums(wgc_list[[s]][, target_pair_select_list, drop = FALSE]))
    )
    colnames(tmp_combine_orig) <- sample_names

    # Keep rows with at least one group meeting both conditions:
    # Count non-zero replicates per group for each region
    group1_nonzero_count <- rowSums(tmp_combine_orig[, group1] != 0)
    group2_nonzero_count <- rowSums(tmp_combine_orig[, group2] != 0)

    # Condition a: more than 50% non-zero in at least one group
    group1_has_majority <- group1_nonzero_count > (case_num * 0.5)
    group2_has_majority <- group2_nonzero_count > (control_num * 0.5)

    # Condition b: at least 2 non-zero in at least one group
    group1_has_min2 <- group1_nonzero_count >= 2
    group2_has_min2 <- group2_nonzero_count >= 2

    keep_rows <- (group1_has_majority & group1_has_min2) | (group2_has_majority & group2_has_min2)
    tmp_combine <- tmp_combine_orig[keep_rows, ]


    # EdgeR object
    d0 <- DGEList(tmp_combine)
    d <- calcNormFactors(d0)
    norm_info <- d$samples
    lib.size <- norm_info$lib.size
    effect_libsize <- norm_info$lib.size * norm_info$norm.factors

    # voom (mean-variance plot)
    voom_plot_filename <- glue("{sig_result_dir}/voom_plot_lowess_span-{lowess_span}_post-filter-one_condition_nonzero-2_column_cluster-{col_label}.pdf")
    pdf(voom_plot_filename)
    y <- voom(d, mm, span = lowess_span, plot = TRUE)
    dev.off()

    mean_log2_cpm <- rowMeans(y$E) + log2(exp(mean(log(lib.size + 1)))) - log2(normalization_factor)

    # # Fit full model
    # fit <- lmFit(y, mm)
    # fitted_values <- fit$coefficients %*% t(mm)
    # residuals <- y$E - fitted_values
    # sqrt_res_std <- sqrt(apply(residuals, 1, sd))

    # # Save "all regions" with p-value/log2FC
    # contr <- makeContrasts(grouptreated - groupuntreated, levels = colnames(coef(fit)))
    # tmp <- contrasts.fit(fit, contr)
    # tmp <- eBayes(tmp)
    # top.table <- topTable(tmp, sort.by = "P", n = Inf)
    # top.table_all <- rownames_to_column(top.table, var = "pos")
    # write_feather(top.table_all, glue("{sig_result_dir}/all_regions_col_cluster-{col_label}.feather"))

    # Loop over mean threshold filters
    for (mean_per_thres in mean_per_thres_list) {
      tmp_combine_filtered <- tmp_combine[mean_log2_cpm >= quantile(mean_log2_cpm, mean_per_thres), ]

      # # Save filtered count matrix
      # # double filter
      # filtered_matrix_filename <- glue("{sig_result_dir}/filtered_counts_col_cluster-{col_label}_rowmean-{mean_per_thres}.feather")
      # write_feather(rownames_to_column(as.data.frame(tmp_combine_filtered), var = "pos"), filtered_matrix_filename)

      # Re-fit with filtered
      d0_filtered <- DGEList(tmp_combine_filtered)
      d_filtered <- calcNormFactors(d0_filtered)

      filtered_voom_plot_filename <- glue("{sig_result_dir}/voom_plot_lowess_span-{lowess_span}_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}.pdf")
      pdf(filtered_voom_plot_filename)
      y_filtered <- voom(d_filtered, mm, span = lowess_span, plot = TRUE)
      dev.off()
      E_value_filtered <- y_filtered$E

      fit <- lmFit(y_filtered, mm)
      contr <- makeContrasts(groupcondition2 - groupcondition1, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contr)
      tmp <- eBayes(tmp)
      top.table <- topTable(tmp, sort.by = "P", n = Inf)

      # save the variance histogram
      E_value_filtered_var <- rowVars(E_value_filtered)
      var_hist_dir_filename <- glue("{sig_result_dir}/hist_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_for_limma.pdf")
      pdf(var_hist_dir_filename)
      hist(E_value_filtered_var, breaks = 100, main = glue("filter: {mean_per_thres}, col cluster: {col_label}\nregion #: {nrow(E_value_filtered)}, dof: {unique(tmp$df.total)}"), xlab = "Variance", col = "blue", border = FALSE)
      dev.off()

      # Plot two histograms
      top.table$P.Value[top.table$P.Value == 0] <- 1E-300
      top.table$adj.P.Val[top.table$adj.P.Val == 0] <- 1E-300
      p_list <- -log10(top.table$P.Value)
      fdr_list <- -log10(top.table$adj.P.Val)

      # Save full table for this threshold
      # significant count matrix
      top.table_for_save <- rownames_to_column(top.table, var = "pos")
      top.table_for_save_dir_filename <- glue("{sig_result_dir}/result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_limma.feather")
      write_feather(top.table_for_save, top.table_for_save_dir_filename)

      # Save significant regions for each FDR threshold
      for (fdr_thres in fdr_thres_list) {
        top.table_sig <- top.table[(top.table$adj.P.Val < fdr_thres) & (abs(top.table$logFC) > l2fc_thres), ]
        sig_region_num <- nrow(top.table_sig)
        if (sig_region_num == 0) {
          print(glue("Will not save for cluster {col_label}!"))
        } else if (sig_region_num > 0) {
          top.table_sig_filename <- glue("{sig_result_dir}/result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_limma_FDR.tsv")
          write.table(top.table_sig, top.table_sig_filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
          sig_region_list <- rownames(top.table_sig)

          for (sample_name in sample_names) {
            wgc_orig_select_pre <- wgc_list[[sample_name]][, target_pair_select_list]
            wgc_orig_select <- wgc_orig_select_pre[sig_region_list, target_pair_select_list]
            wgc_log2_select <- log2(wgc_orig_select + pseudocount_for_log)
            wgc_log2_select <- rownames_to_column(wgc_log2_select, var = "pos")
            wgc_log2_select_filename <- glue("{sample_name}_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
            wgc_log2_select_dir_filename <- file.path(sig_result_dir, wgc_log2_select_filename)
            write_feather(wgc_log2_select, wgc_log2_select_dir_filename)
          }
        }
      }
    }
  }
}
