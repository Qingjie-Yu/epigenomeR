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
