# Biclustering Heatmap Visualization
# Post: Create a heatmap with predefined row and column clustering assignments.
#       Automatically calculates optimal cell sizes based on label dimensions.
# Parameters: mat: a numeric matrix with rownames and colnames
#             row_cluster_file_path: path to row cluster assignment .tsv file
#                                   (required columns: feature, label)
#             col_cluster_file_path: path to column cluster assignment .tsv file
#                                   (required columns: feature, label)
#             out_dir: directory to save output PDF (default: "./")
#             show_column_names: whether to show column names at bottom
#             lower_range: minimum value for color scale (default: NULL, auto-calculated)
#             upper_range: maximum value for color scale (default: NULL, auto-calculated)
#             row_title_fontsize: font size for row cluster titles (default: 40)
#             col_title_fontsize: font size for column cluster titles (default: 22)
#             legend_title_fontsize: font size for legend title (default: 40)
#             legend_label_fontsize: font size for legend labels (default: 30)
# Output: saves a PDF file named "biclustering_heatmap.pdf" in out_dir
#         returns the Heatmap object for further customization

biclustering_heatmap <- function(mat, row_cluster_file_path, col_cluster_file_path, out_dir = "./", show_column_names = FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  # Load Library
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

  calculate_ht_size <- function(ht, heatmap_legend_side = NULL, unit = "inch") {
    pdf(NULL)
    ht <- draw(ht, heatmap_legend_side = heatmap_legend_side, background = "transparent")
    w <- ComplexHeatmap:::width(ht)
    w <- convertX(w, unit, valueOnly = TRUE)
    h <- ComplexHeatmap:::height(ht)
    h <- convertY(h, unit, valueOnly = TRUE)
    dev.off()
    c(w, h)
  }

  calculate_cell_size <- function(row_labels, col_labels, row_fontsize, col_fontsize, default_cell_width = 5, safety_factor = 1.5, lower_quantile = 0.8, upper_quantile = 0.5) {
    default_cell_height <- default_cell_width * length(col_labels) / length(row_labels)

    row_counts <- table(row_labels)
    unique_row_labels <- names(row_counts)
    row_text_height <- row_fontsize * safety_factor
    lower_count_row <- quantile(row_counts, lower_quantile, type = 1)
    min_cell_height_lower <- row_text_height / lower_count_row
    upper_count_row <- quantile(row_counts, upper_quantile, type = 1)
    max_cell_height_upper <- row_text_height * 4 / upper_count_row

    col_counts <- table(col_labels)
    unique_col_labels <- names(col_counts)
    col_text_widths <- sapply(unique_col_labels, function(label) {
      nchar(as.character(label)) * col_fontsize * 0.6 * safety_factor
    })
    lower_count_col <- quantile(col_counts, lower_quantile, type = 1)
    min_cell_width_lower <- max(col_text_widths / lower_count_col)
    upper_count_col <- quantile(col_counts, upper_quantile, type = 1)
    max_cell_width_upper <- min(col_text_widths * 4 / upper_count_col)

    if (min_cell_height_lower > max_cell_height_upper) {
      final_cell_height <- min_cell_height_lower
    } else {
      final_cell_height <- max(default_cell_height, min_cell_height_lower)
      final_cell_height <- min(final_cell_height, max_cell_height_upper)
    }

    if (min_cell_width_lower > max_cell_width_upper) {
      final_cell_width <- min_cell_width_lower
    } else {
      final_cell_width <- max(default_cell_width, min_cell_width_lower)
      final_cell_width <- min(final_cell_width, max_cell_width_upper)
    }

    list(
      cell_width_mm = final_cell_width,
      cell_height_mm = final_cell_height
    )
  }

  calculate_legend_fontsize <- function(col_labels, cell_width_mm, legend_title = "Normalized Read Counts", safety_factor = 1.2) {
    available_width_mm <- length(col_labels) * cell_width_mm * safety_factor
    title_chars <- nchar(legend_title)
    max_title_fontsize <- available_width_mm / (title_chars * 0.2)
    max_title_fontsize
  }

  # load row cluster info
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  row_cluster <- row_cluster[row_cluster$feature %in% rownames(mat), ]
  row_order <- row_cluster$feature
  row_split <- row_cluster$label

  # load col cluster info
  col_cluster <- read.table(col_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  col_cluster <- col_cluster[col_cluster$feature %in% colnames(mat), ]
  col_order <- col_cluster$feature
  col_split <- as.numeric(factor(col_cluster$label))

  mat <- mat[row_order, col_order]
  if (is.null(lower_range) || lower_range == "") {
    lower_range <- min(mat, na.rm = TRUE)
  }
  if (is.null(upper_range) || upper_range == "") {
    upper_range <- max(mat, na.rm = TRUE)
  }
  avg <- (lower_range + upper_range) / 2
  col_fun <- colorRamp2(c(lower_range, avg, upper_range), c("#3155C3", "white", "#AF0525"))

  if (is.null(row_title_fontsize)) {
    row_title_fontsize <- 20
  }
  if (is.null(col_title_fontsize)) {
    col_title_fontsize <- 20
  }
  if (is.null(legend_title_fontsize)) {
    legend_title_fontsize <- 15
  }
  if (is.null(legend_label_fontsize)) {
    legend_label_fontsize <- 15
  }

  cell_size <- calculate_cell_size(row_labels = row_split, col_labels = col_split, row_fontsize = row_title_fontsize, col_fontsize = col_title_fontsize)
  cell_width <- cell_size$cell_width_mm
  cell_height <- cell_size$cell_height_mm

  max_legend_title_fontsize <- calculate_legend_fontsize(col_labels = col_split, cell_width_mm = cell_width)
  legend_title_fontsize <- min(legend_label_fontsize, max_legend_title_fontsize)

  ht <- Heatmap(
    mat,
    col = col_fun,
    # Clustering setting
    cluster_columns = FALSE, cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_order = row_order, column_order = col_order,
    row_dend_reorder = FALSE,
    row_split = row_split, column_split = col_split,
    show_row_dend = FALSE, show_column_dend = FALSE,

    # General setting
    width = ncol(mat) * unit(cell_width, "mm"), height = nrow(mat) * unit(cell_height, "mm"),
    row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
    row_title_gp = gpar(fontsize = row_title_fontsize), column_title_gp = gpar(fontsize = col_title_fontsize),
    show_row_names = FALSE, show_column_names = show_column_names,
    column_labels = TeX(colnames(mat)),
    row_title_rot = 0,
    use_raster = TRUE,

    # Cell setting
    border = TRUE,
    rect_gp = gpar(col = NA, lwd = 0),
    border_gp = gpar(col = "white", lwd = 0),

    # Legend setting
    heatmap_legend_param = list(
      title = "Normalized Read Counts",
      legend_width = ncol(mat) / 2 * unit(cell_width, "mm"),
      grid_height = 2 * unit(cell_width, "mm"),
      title_position = "topcenter",
      title_gp = gpar(fontsize = legend_title_fontsize),
      labels_gp = gpar(fontsize = legend_label_fontsize),
      legend_direction = "horizontal"
    )
  )
  ht_size <- calculate_ht_size(ht, heatmap_legend_side = "bottom")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  out_path <- file.path(out_dir, "biclustering_heatmap.pdf")
  pdf(out_path, width = ht_size[1], height = ht_size[2])
  draw(ht, heatmap_legend_side = "bottom", background = "transparent")
  dev.off()
  message("Heatmap saved to: ", out_path)
}
