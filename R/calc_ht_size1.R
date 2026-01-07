calc_ht_size1 <- function(ht, unit = "inch", show_annotation_legend = FALSE, column_title = NULL, column_title_fontsize = 14) {
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
