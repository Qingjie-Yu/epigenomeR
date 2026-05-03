# Clustering Heatmap Visualization
# Post: Create a heatmap with predefined row clustering assignments.
#       Automatically calculates optimal cell sizes based on label dimensions.
# Parameters: mat: a numeric matrix with rownames and colnames
#             row_cluster_file_path: path to row cluster assignment .tsv file
#                                   (required columns: region, cluster)
#             out_dir: directory to save output PDF (default: "./")
#             show_column_names: whether to show column names (default: FALSE)
#             column_names_rot: rotation angle for column names (default: 45)
#             row_title_fontsize: font size for row cluster titles (default: 8)
#             col_title_fontsize: font size for column titles (default: 8)
#             legend_title_fontsize: font size for legend title (default: auto-calculated)
#             legend_label_fontsize: font size for legend labels (default: 6)
#             cell_width: cell width in inches (default: NULL, auto-calculated)
#             cell_height: cell height in inches (default: NULL, auto-calculated)
#             fig_width: output PDF width in inches (default: NULL, auto-calculated)
#             fig_height: output PDF height in inches (default: NULL, auto-calculated)
# Output: saves a PDF file named "biclustering_heatmap.pdf" in out_dir
#         returns the Heatmap object invisibly


clustering_heatmap <- function(mat, row_cluster_file_path, out_dir = "./", pdf_name = "clustering_heatmap.pdf", show_column_names = FALSE, fig_width = NULL, fig_height = NULL, cell_width = 0.5 / 2.54, cell_height = 0.003 / 2.54, row_title_fontsize = 8, col_title_fontsize = 8, legend_title_fontsize = 6, legend_label_fontsize = 6, column_names_rot = 45) {
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

  # Load row cluster info
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t", row.names = NULL)
  row_cluster <- row_cluster[row_cluster$region %in% rownames(mat), ]
  row_order   <- row_cluster$region
  row_split   <- row_cluster$cluster
  mat         <- mat[row_order, , drop = FALSE]

  # Color scale: blue → white → red
  lo      <- min(mat, na.rm = TRUE)
  hi      <- max(mat, na.rm = TRUE)
  avg     <- (lo + hi) / 2
  col_fun <- colorRamp2(c(lo, avg, hi), c("#3155C3", "white", "#AF0525"))

  # Build heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,

    cluster_columns    = FALSE,
    cluster_rows       = FALSE,
    cluster_row_slices = FALSE,
    row_order          = row_order,
    column_order       = colnames(mat),
    row_dend_reorder   = FALSE,
    row_split          = row_split,
    show_row_dend      = FALSE,
    show_column_dend   = FALSE,

    width             = ncol(mat) * unit(cell_width,  "inches"),
    height            = nrow(mat) * unit(cell_height, "inches"),
    row_gap           = unit(1, "mm"),
    column_gap        = unit(1, "mm"),
    row_title_gp      = gpar(fontsize = row_title_fontsize),
    column_title_gp   = gpar(fontsize = col_title_fontsize),
    show_row_names    = FALSE,
    show_column_names = show_column_names,
    column_names_rot  = column_names_rot,
    column_labels     = TeX(colnames(mat)),
    row_title_rot     = 0,
    use_raster        = TRUE,

    border    = TRUE,
    rect_gp   = gpar(col = NA, lwd = 0),
    border_gp = gpar(col = "white", lwd = 0),

    heatmap_legend_param = list(
      title            = "Z-Score",
      title_position   = "topcenter",
      legend_direction = "horizontal",
      legend_width     = ncol(mat) / 2 * unit(cell_width, "inches"),
      grid_height      = 2 * unit(cell_width, "inches"),
      title_gp         = gpar(fontsize = legend_title_fontsize),
      labels_gp        = gpar(fontsize = legend_label_fontsize)
    )
  )

  # Save PDF
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_path <- file.path(out_dir, pdf_name)

  pdf(NULL)
  drawn <- draw(ht, heatmap_legend_side = "right", background = "transparent")
  w <- convertX(ComplexHeatmap:::width(drawn),  "inches", valueOnly = TRUE)
  h <- convertY(ComplexHeatmap:::height(drawn), "inches", valueOnly = TRUE)
  dev.off()

  if (!is.null(fig_width))  w <- fig_width
  if (!is.null(fig_height)) h <- fig_height

  pdf(out_path, width = w, height = h)
  draw(ht, heatmap_legend_side = "right", background = "transparent")
  dev.off()
  message("Heatmap saved to: ", out_path)
  invisible(ht)
}

clustering_multi_heatmap <- function(
  crf_mat_list,             # named list of matrices, each peaks × samples, names = CRF names
  row_cluster_file_path,    # path to row cluster .tsv (columns: region, cluster)
  out_dir = "./",
  pdf_name = "clustering_multi_heatmap",
  show_column_names = TRUE,
  column_names_rot = 45     # rotation angle for sample names
) {

  # ============================================================
  # Visual parameters (edit here, no need to pass as arguments)
  # ============================================================

  # Font sizes (pt)
  row_title_fontsize      <- 8   # pt, row cluster label (e.g. "Cluster 1")
  col_title_fontsize      <- 8   # pt, CRF title on top (e.g. "SMARCA4")
  legend_title_fontsize   <- 6    # pt, legend title "Z-Score"
  legend_label_fontsize   <- 5    # pt, legend tick labels
  column_names_fontsize   <- 6    # pt, sample names at bottom

  # Figure sizes
  fig_width <- 17 / 2.54
  fig_height <- NULL

  # Cell dimensions (inches → divide cm value by 2.54)
  cell_width   <- 0.5 / 2.54   # width per cell
  cell_height  <- 0.001 / 2.54 # height per cell (peaks are many, keep small)

  # Gap between CRF heatmaps (inches)
  heatmap_gap_inch  <- 0.3 / 2.54

  # Legend dimensions (inches)
  legend_width_inch  <- 0.3 / 2.54    # total legend bar width
  legend_height_inch <- 5 / 2.54  # legend bar height

  # Color scale range (NULL = auto per CRF)
  lower_range <- -2
  upper_range <-  2

  # Row gap between clusters (mm)
  row_gap_mm <- 1

  # Per-CRF palettes: white → deep color (single hue)
  crf_palettes <- list(
    c("white", "#AF0525"),   # CRF1: white-red
    c("white", "#c46e00"),   # CRF2: white-orange
    c("white", "#1a6e2e"),   # CRF3: white-green
    c("white", "#0277bd"),   # CRF4: white-green blue
    c("white", "#3155C3"),   # CRF5: white-blue
    c("white", "#7b1fa2")    # CRF6: white-puper
  )

  # ============================================================
  # Load libraries
  # ============================================================
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })

  # ============================================================
  # Load row cluster info
  # ============================================================
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
  # filter to peaks present in matrices
  all_peaks <- rownames(crf_mat_list[[1]])
  row_cluster <- row_cluster[row_cluster$region %in% all_peaks, ]
  row_order  <- row_cluster$region
  row_split  <- row_cluster$cluster   # character vector e.g. "A","B",...

  crf_names <- names(crf_mat_list)
  n_crf     <- length(crf_names)

  # ============================================================
  # Build one Heatmap object per CRF
  # ============================================================
  ht_list <- NULL

  for (k in seq_len(n_crf)) {
    crf  <- crf_names[k]
    mat  <- crf_mat_list[[crf]]

    # reorder rows
    mat  <- mat[row_order, , drop = FALSE]

    # color function for this CRF
    pal_idx <- ((k - 1) %% length(crf_palettes)) + 1
    pal     <- crf_palettes[[pal_idx]]
    lo      <- if (is.null(lower_range)) min(mat, na.rm = TRUE) else lower_range
    hi      <- if (is.null(upper_range)) max(mat, na.rm = TRUE) else upper_range
    col_fun <- colorRamp2(c(lo, hi), pal)

    ht <- Heatmap(
      mat,
      col  = col_fun,
      name = crf,   # used as legend title

      # pre-computed order, no re-clustering
      cluster_rows         = FALSE,
      cluster_columns      = FALSE,
      cluster_row_slices   = FALSE,
      row_order            = row_order,
      column_order         = colnames(mat),
      row_split            = row_split,

      # titles
      column_title     = crf,
      column_title_gp  = gpar(fontsize = col_title_fontsize, fontface = "bold"),
      row_title_gp     = gpar(fontsize = row_title_fontsize),
      row_title_rot    = 0,
      row_title_side   = "left",

      # show row cluster labels only on leftmost heatmap
      show_row_names   = FALSE,
      show_column_names = show_column_names,
      column_names_gp  = gpar(fontsize = column_names_fontsize),
      column_names_rot = column_names_rot,

      # dendrograms
      show_row_dend    = FALSE,
      show_column_dend = FALSE,

      # cell dimensions
      width  = ncol(mat) * unit(cell_width,  "inches"),
      height = nrow(mat) * unit(cell_height, "inches"),

      # gaps
      row_gap    = unit(row_gap_mm, "mm"),
      column_gap = unit(0, "mm"),

      # borders
      border    = TRUE,
      rect_gp   = gpar(col = NA),
      border_gp = gpar(col = "white", lwd = 0.5),

      use_raster = TRUE,

      # per-CRF legend
      heatmap_legend_param = list(
        title            = "Z-Score",
        title_position   = "topcenter",
        legend_direction = "vertical",
        legend_width     = unit(legend_width_inch, "inches"),
        grid_height      = unit(legend_height_inch, "inches"),
        title_gp         = gpar(fontsize = legend_title_fontsize, fontface = "bold"),
        labels_gp        = gpar(fontsize = legend_label_fontsize)
      )
    )

    if (is.null(ht_list)) {
      ht_list <- ht
    } else {
      ht_list <- ht_list + ht
    }
  }

  # ============================================================
  # Draw and save
  # ============================================================
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_path <- file.path(out_dir, pdf_name)

  # calculate size
  pdf(NULL)
  drawn <- draw(
    ht_list,
    heatmap_legend_side  = "right",
    ht_gap               = unit(heatmap_gap_inch, "inches"),
    background           = "transparent",
    merge_legend         = FALSE
  )
  w <- ComplexHeatmap:::width(drawn)
  h <- ComplexHeatmap:::height(drawn)
  w <- convertX(w, "inches", valueOnly = TRUE)
  h <- convertY(h, "inches", valueOnly = TRUE)
  dev.off()

  if (!is.null(fig_width))  w <- fig_width
  if (!is.null(fig_height)) h <- fig_height
  pdf(out_path, width = w, height = h)
  draw(
    ht_list,
    heatmap_legend_side  = "right",
    ht_gap               = unit(heatmap_gap_inch, "inches"),
    background           = "transparent",
    merge_legend         = FALSE
  )
  dev.off()

  message("Multi-CRF heatmap saved to: ", out_path)
  invisible(ht_list)
}