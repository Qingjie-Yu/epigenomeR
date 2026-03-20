differential_heatmap_single <- function(counts_path, cmp_tag, grp_name, out_dir = "./", show_colnames = FALSE, col_width_mm = 10, row_height_mm = 0.4, random_seed = 42) {
  suppressPackageStartupMessages({
    library(glue)
    library(ComplexHeatmap)
    library(circlize)
    library(data.table)
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dt        <- data.table::fread(counts_path)
  pos       <- as.character(dt$pos)
  dt$pos    <- NULL
  all_cols  <- colnames(dt)
  count_mat <- as.matrix(dt)
  mode(count_mat) <- "numeric"
  rownames(count_mat) <- pos

  if (nrow(count_mat) == 0) {
    message(glue("| {cmp_tag} | {grp_name} | empty counts file, skipping heatmap"))
    return(invisible(NULL))
  }

  log2_mat <- log2(count_mat + 1)

  # auto-detect cluster mode
  use_cluster  <- any(grepl(":", all_cols, fixed = TRUE))
  sample_names <- if (use_cluster) unique(sub(":.*$", "", all_cols)) else all_cols
  col_split    <- factor(sub(":.*$", "", all_cols), levels = sample_names)

  # column title fontsize: must fit within the narrowest cluster (80% of its width)
  min_cluster_width_pt <- min(table(col_split)) * col_width_mm * 2.835
  fontsize <- min(
    min_cluster_width_pt * 0.8,                                     # cluster width constraint
    min_cluster_width_pt * 0.8 / (max(nchar(sample_names)) * 0.6), # name length constraint
    72                                                               # hard ceiling
  )
  fontsize <- max(fontsize, 6)

  col_fun <- colorRamp2(
    c(min(log2_mat), mean(log2_mat), max(log2_mat)),
    c("#3155C3", "white", "#AF0525")
  )

  ht <- Heatmap(log2_mat,
    name                = "log2",
    col                 = col_fun,
    show_row_names      = FALSE,
    show_column_names   = show_colnames,
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    column_split        = col_split,
    width               = ncol(log2_mat) * unit(col_width_mm,  "mm"),
    height              = nrow(log2_mat) * unit(row_height_mm, "mm"),
    column_title        = sample_names,
    column_title_gp     = gpar(fontsize = fontsize),
    column_gap          = unit(8, "mm"),
    show_row_dend       = FALSE,
    show_heatmap_legend = FALSE,
    use_raster          = TRUE
  )

  pdf_w <- (ncol(log2_mat) * col_width_mm + (length(sample_names) - 1) * 8) / 25.4
  pdf_h <- nrow(log2_mat) * row_height_mm / 25.4 + fontsize / 72 * 1.5

  pdf_filename <- file.path(out_dir, glue("{cmp_tag}_{grp_name}.pdf"))
  pdf(pdf_filename, width = pdf_w, height = pdf_h)
  set.seed(random_seed)
  draw(ht, background = "transparent", show_annotation_legend = FALSE)
  dev.off()

  message(glue("| {cmp_tag} | {grp_name} | {nrow(log2_mat)} regions | fontsize: {round(fontsize, 1)}"))
  invisible(pdf_filename)
}


differential_heatmap <- function(dr_result, out_dir = "./", show_colnames = FALSE, col_width_mm  = 10, row_height_mm = 0.4, random_seed = 42) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (grp_name in names(dr_result)) {
    for (cmp_tag in names(dr_result[[grp_name]])) {
      counts_path <- dr_result[[grp_name]][[cmp_tag]]
      if (is.null(counts_path) || !file.exists(counts_path)) next
      differential_heatmap_single(counts_path = counts_path, cmp_tag = cmp_tag, grp_name = grp_name, out_dir = out_dir, show_colnames = show_colnames, col_width_mm  = col_width_mm, row_height_mm = row_height_mm, random_seed = random_seed)
    }
  }
}