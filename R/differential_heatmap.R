calc_ht_size <- function(ht, unit = "inch", show_annotation_legend = FALSE, column_title = NULL, column_title_fontsize = 14) {
  pdf(NULL)
  if (show_annotation_legend) {
    ht <- draw(ht, background = "transparent", column_title = column_title, column_title_gp = gpar(fontsize = column_title_fontsize), merge_legend = TRUE, annotation_legend_side = "top")
  } else {
    ht <- draw(ht, background = "transparent", column_title = column_title, column_title_gp = gpar(fontsize = column_title_fontsize), show_annotation_legend = FALSE)
  }
  ht <- draw(ht, background = "transparent", column_title = column_title, column_title_gp = gpar(fontsize = column_title_fontsize))
  w <- ComplexHeatmap:::width(ht)
  w <- convertX(w, unit, valueOnly = TRUE)
  h <- ComplexHeatmap:::height(ht)
  h <- convertY(h, unit, valueOnly = TRUE)
  dev.off()

  c(w, h)
}


differential_heatmap <- function(da_cm_path, out_dir, show_colnames = FALSE, col_size_coef = 20, width_base = 8, random_seed = 42) {
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fnames <- basename(da_cm_path)
  sample_names <- sub(".*_cluster\\d+_(.+)_log2\\.feather$", "\\1", fnames)
  comparison <- sub("^(.+)_cluster\\d+_.+_log2\\.feather$", "\\1", fnames[1])
  cluster_idx  <- sub(".*_cluster(\\d+)_.+_log2\\.feather$", "\\1", fnames[1])

  wgc_list <- lapply(seq_along(da_cm_path), function(i) {
    df <- read_feather(da_cm_path[i])
    rownames(df) <- df$pos
    df$pos <- NULL
    as.data.frame(df)
  })
  names(wgc_list) <- sample_names

  wgc_prefixed <- lapply(seq_along(sample_names), function(i) {
    m <- wgc_list[[i]]
    colnames(m) <- paste0(sample_names[i], ":", colnames(m))
    m
  })
  wgc_log2_cbind <- do.call(cbind, lapply(wgc_prefixed, as.matrix))

  col_order <- colnames(wgc_log2_cbind)
  col_split <- factor(
    sub("^(.+):.+$", "\\1", col_order),
    levels = sample_names
  )

  col_fun <- colorRamp2(
    c(min(wgc_log2_cbind), mean(wgc_log2_cbind), max(wgc_log2_cbind)),
    c("#3155C3", "white", "#AF0525")
  )

  ht <- Heatmap(as.matrix(wgc_log2_cbind),
    name                = "log2",
    col                 = col_fun,
    show_row_names      = FALSE,
    show_column_names   = show_colnames,
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    column_split        = col_split,
    width               = col_size_coef * unit(width_base, "mm"),
    height              = (3750 / 89) * unit(10, "mm"),
    column_title        = sample_names,
    column_title_gp     = gpar(fontsize = 90),
    column_gap          = unit(8, "mm"),
    show_row_dend       = FALSE,
    show_heatmap_legend = FALSE,
    use_raster          = TRUE
  )

  set.seed(random_seed)
  size <- calc_ht_size(ht, unit = "inch", show_annotation_legend = FALSE)

  pdf_filename <- file.path(out_dir,
    glue("{comparison}_cluster{cluster_idx}.pdf"))

  pdf(pdf_filename, width = 1.01 * size[1], height = 1.01 * size[2])
  set.seed(random_seed)
  draw(ht, background = "transparent", show_annotation_legend = FALSE)
  dev.off()

  message(glue("cluster {cluster_idx} | {comparison} | {nrow(wgc_log2_cbind)} rows"))
}