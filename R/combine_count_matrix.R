# Combine cluster dataframe to count matrix
# Post: Apply clustering order to a count matrix and plot heatmap.
# Parameters: count_matrix: a df matrix with rownames = pos
#             row_cluster_file_path: path to row cluster assignment .tsv
#             col_cluster_file_path: path to column cluster assignment .tsv
#             basename: base name for output files
#             output_dir_path: directory to save plot
#             seed: random seed (default is 123).
#             show_dend_boolean: whether to show column names (default to FALSE).
# Output: saves a single heatmap as PDF
combine_count_matrix <- function(count_matrix, row_cluster_file_path, col_cluster_file_path, basename, target_pair_mapping_df_path = NULL, output_dir_path = NULL, seed = 123, show_dend_boolean = FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = 40, col_title_fontsize = 22, legend_title_fontsize = 40, legend_label_fontsize = 30) {
  # load library
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("The latex2exp package is required but not installed.")
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(tibble)
    library(arrow)
    library(latex2exp)
    library(svglite)
    library(viridis)
    library(glue)
    library(dplyr)
    library(tools)
  })
  # source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
  calc_ht_size <- function(ht, heatmap_legend_side = NULL, unit = "inch") {
    pdf(NULL)
    ht <- draw(ht, heatmap_legend_side = heatmap_legend_side, background = "transparent")
    w <- ComplexHeatmap:::width(ht)
    w <- convertX(w, unit, valueOnly = TRUE)
    h <- ComplexHeatmap:::height(ht)
    h <- convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
  }
  # load row cluster info
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  row_order <- row_cluster$feature
  row_split <- row_cluster$label

  # load col cluster info
  col_cluster <- read.table(col_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  col_order <- col_cluster$feature
  col_order_shorten <- map_target_names(col_cluster$feature, target_pair_mapping_df_path = target_pair_mapping_df_path)

  # inner_join
  col_order_valid <- intersect(col_order, colnames(count_matrix))
  row_order_valid <- intersect(row_order, rownames(count_matrix))
  combined_df <- count_matrix[row_order_valid, col_order_valid]

  filtered_col_cluster <- col_cluster[col_cluster$feature %in% colnames(combined_df), ]
  col_order_shorten <- filtered_col_cluster$feature
  col_split <- as.numeric(filtered_col_cluster$label)

  input_name <- basename

  if (is.null(lower_range) || lower_range == "") {
    lower_range <- min(combined_df, na.rm = TRUE)
  } else {
    lower_range <- lower_range
  }
  if (is.null(upper_range) || upper_range == "") {
    upper_range <- max(combined_df, na.rm = TRUE)
  } else {
    upper_range <- upper_range
  }
  avg <- (lower_range + upper_range) / 2
  col_fun <- colorRamp2(c(lower_range, avg, upper_range), c("#3155C3", "white", "#AF0525"))

  set.seed(seed)
  V1V2_ht <- Heatmap(as.matrix(combined_df),
    name = input_name,
    col = col_fun,
    cluster_columns = FALSE, cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_dend_reorder = FALSE,
    row_order = row_order,
    row_split = row_split,
    show_row_dend = FALSE,
    column_order = col_order_shorten,
    column_split = col_split,
    border = TRUE,
    rect_gp = gpar(col = NA, lwd = 0),
    border_gp = gpar(col = "white", lwd = 0),
    show_row_names = FALSE,
    show_column_names = show_dend_boolean,
    width = ncol(combined_df) * unit(5, "mm"),
    height = ncol(combined_df) * unit(5, "mm"),
    row_gap = unit(3, "mm"),
    column_gap = unit(3, "mm"),
    # row_gap=unit(0, "mm"),
    # column_gap=unit(0, "mm"),
    row_title_gp = gpar(fontsize = row_title_fontsize), # module size
    column_title_gp = gpar(fontsize = col_title_fontsize),
    row_title_rot = 0,
    show_heatmap_legend = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    heatmap_legend_param = list(
      title = "Normalized Read Counts",
      legend_width = ncol(combined_df) / 2 * unit(5, "mm"),
      grid_height = 2 * unit(5, "mm"),
      title_position = "topcenter",
      title_gp = gpar(fontsize = legend_title_fontsize),
      labels_gp = gpar(fontsize = legend_label_fontsize),
      legend_direction = "horizontal"
    ),
    column_labels = TeX(colnames(combined_df)),
    use_raster = TRUE
  )

  set.seed(seed)
  V1V2_size <- calc_ht_size(V1V2_ht, heatmap_legend_side = "bottom")

  # save file
  fig_dir <- file.path(output_dir_path, "figures")

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  cluster_heatmap_dir_filename <- file.path(fig_dir, paste0(input_name, "_figS3A.pdf"))
  pdf(cluster_heatmap_dir_filename, width = V1V2_size[1], height = V1V2_size[2])
  set.seed(seed)
  V1V2_ht_draw <- draw(V1V2_ht, heatmap_legend_side = "bottom", background = "transparent")
  dev.off()
}
