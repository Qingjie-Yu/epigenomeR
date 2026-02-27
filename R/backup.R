
# Backup
limma_column_cluster_differential_regions <- function(sample_names, case_num, control_num, wgc_file_path, sig_result_dir, col_cluster_file = NULL, normalization_factor = 1E6, lowess_span = 0.5, l2fc_thres = 0.5, mean_per_thres_list = c(0.25), fdr_thres_list = c(0.25), pseudocount_for_log = 1) {
  
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



# differential -2
# Post: Generate differential expression heatmaps for significant regions by reading pre-computed log2-transformed expression matrices, visualizing expression patterns across conditions with customizable formatting and clustering options.
# Parameter: sample_names: Vector of sample names to be displayed as column titles and used for data organization
#            load_dir: Directory containing log2-transformed expression matrix feather files from limma_column_cluster_differential_regions
#            sig_result_dir: Output directory path where heatmap PDF files will be saved
#            col_cluster_file: Path to TSV file containing column cluster assignments with 'feature' and 'label' columns (default: NULL, treat each col as a cluster)
#            cluster_idx_list: Vector of cluster indices to process (default: NULL, uses all unique labels from col_cluster_file)
#            mean_per_thres: Mean expression percentile threshold used in differential analysis (default: "0.25")
#            fdr_thres: FDR threshold used in differential analysis (default: "0.25")
#            l2fc_thres: Log2 fold change threshold used in differential analysis (default: 0.5)
#            show_heatmap_legend: Whether to display heatmap legend ("on"/"off", default: "off")
#            show_colnames: Whether to display column names ("on"/"off", default: "off")
#            col_size_coef: Coefficient for adjusting column width (default: 20)
#            colnames_fontsize: Font size for column names (default: 10)
#            width_base: Base width in mm for heatmap sizing (default: 8)
#            random_seed: Random seed for reproducible raster rendering (default: 42)
#            font_size: Font size for various text elements (default: 50)
#            target_pair_mapping_df_path: Path to target name mapping file, or NULL for no mapping (default: NULL)
# Output: Saves PDF heatmap files for each cluster showing log2 expression values with blue-white-red color scheme, organized by sample groups
plot_diff_heatmaps <- function(sample_names, load_dir, sig_result_dir, col_cluster_file = NULL, cluster_idx_list = NULL, show_heatmap_legend = "off", show_colnames = "off", col_size_coef = 20, colnames_fontsize = 10, width_base = 8, random_seed = 42, font_size = 50, target_pair_mapping_df_path = NULL, mean_per_thres = "0.25", fdr_thres = "0.25", l2fc_thres = 0.5) {
  # Load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(svglite)
    library(ComplexHeatmap)
    library(circlize)
    library(tidyr)
    library(dplyr)
    library(latex2exp)
  })

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  # Handle column cluster file
  if (is.null(col_cluster_file)) {
    # Get all unique column names from all WGC files
    all_features <- unique(unlist(lapply(wgc_list, colnames)))
    # Create a data frame where each feature gets its own label
    col_cluster_full <- data.frame(
      feature = all_features,
      label = seq_along(all_features),
      stringsAsFactors = FALSE
    )
  } else {
    # Load column cluster file
    col_cluster_full <- read.table(col_cluster_file, header = TRUE, sep = "\t", row.names = NULL)
  }

  if (is.null(cluster_idx_list)) {
    cluster_idx_list <- sort(unique(col_cluster_full$label))
    message(glue("cluster_idx_list not specified. Using all unique labels from column cluster file: {paste(cluster_idx_list, collapse=', ')}"))
  }

  for (cluster_idx in cluster_idx_list) {
    wgc_file_path <- glue("{load_dir}/{sample_names}_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
    wgc_list <- lapply(wgc_file_path, function(f) {
      df <- column_to_rownames(read_feather(f), var = "pos")
      colnames(df) <- map_target_names(colnames(df), target_pair_mapping_df_path = target_pair_mapping_df_path)
      df
    })
    names(wgc_list) <- sample_names

    df_list <- lapply(seq_along(sample_names), function(i) {
      get_cluster_df(wgc_list[[i]], sample_names[i])
    })
    col_cluster_df <- bind_rows(df_list) %>%
      mutate(label = factor(label, levels = sample_names)) %>%
      mutate(order = match(nonprefix, col_cluster_full$feature)) %>%
      arrange(label, order)

    col_order <- col_cluster_df$feature
    col_split <- col_cluster_df$label

    for (s in sample_names) {
      colnames(wgc_list[[s]]) <- paste0(s, ":", colnames(wgc_list[[s]]))
    }

    wgc_log2_cbind <- do.call(cbind, lapply(wgc_list, as.matrix))
    wgc_log2_cbind <- wgc_log2_cbind[, col_order]

    col_fun <- colorRamp2(c(min(wgc_log2_cbind), 0.9, 1.8), c("#3155C3", "white", "#AF0525"))

    show_colnames_bool <- (show_colnames == "on")

    heatmap_prefix <- glue("diff_heatmap_col_cluster-{cluster_idx}_colname-{show_colnames}_col-reorder_size-{col_size_coef}")
    col_num <- ncol(wgc_log2_cbind)

    ht <- Heatmap(as.matrix(wgc_log2_cbind),
      name = "log2",
      col = col_fun,
      show_row_names = FALSE,
      show_column_names = show_colnames_bool,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_order = col_order,
      column_split = col_split,
      width = col_size_coef * unit(width_base, "mm"),
      height = (3750 / 89) * unit(10, "mm"),
      column_title = sample_names,
      column_title_gp = gpar(col = c("#3155C3", "#3155C3", "#AF0525", "#AF0525"), fontsize = 90),
      column_gap = unit(8, "mm"),
      column_names_gp = gpar(fontsize = colnames_fontsize),
      show_row_dend = FALSE,
      show_heatmap_legend = FALSE,
      heatmap_legend_param = list(
        title = "log2",
        grid_width = 5 * unit(5, "mm"),
        legend_height = 77 * 1.3 * unit(5, "mm"),
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 40),
        legend_direction = "vertical"
      ),
      use_raster = TRUE
    )
    set.seed(random_seed)
    size <- calc_ht_size1(ht, unit = "inch", show_annotation_legend = FALSE)

    # save
    pdf_h_heatmap_filename <- glue("{heatmap_prefix}.pdf")
    pdf_h_heatmap_dir_filename <- file.path(sig_result_dir, pdf_h_heatmap_filename)
    pdf(pdf_h_heatmap_dir_filename, width = 1.01 * size[1], height = 1.01 * size[2])
    set.seed(random_seed)
    draw(ht, background = "transparent", show_annotation_legend = FALSE)
    dev.off()

    message(glue("cluster_idx: {cluster_idx}; {nrow(wgc_log2_cbind)} rows"))
  }
}
