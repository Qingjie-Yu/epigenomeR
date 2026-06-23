# Top-Percentile Region Filter
#
# Alternative to detect_hvr(). For each CRF pair (column) in the transformed count
# matrix, computes the value at the top `top_pct` threshold (e.g. top_pct = 0.01 ->
# the 99th percentile of that column), then takes the UNION across all pairs of bins
# whose value meets or exceeds their column's threshold. The transformed matrix is
# subsetted to that union of bins (all pair columns retained, transformed values
# unchanged) and written out in the same feather format as filter_hvr()'s output,
# so it is a drop-in alternative wherever a filtered_cm_path is expected downstream.
#
# Parameters:
#   transformed_cm_path: Path to the transformed count matrix `.feather` file
#                         (output of apply_transformations()).
#   out_dir : Directory to save the filtered matrix (and diagnostic plot, if requested).
#   top_pct : Numeric in (0, 1). Fraction defining the "top" cutoff per pair.
#             e.g. 0.01 keeps, for each pair, the bins in that pair's top 1% by
#             transformed value. Default: 0.01
#   plot    : Logical. Whether to generate a diagnostic plot showing, per pair, the
#             top_pct threshold value and number of bins selected. Default: FALSE
#
# Output: Writes a Feather file (same layout as filter_hvr(): first column 'pos',
#         remaining columns = transformed values per pair, restricted to the union
#         of top-pct bins across all pairs). Returns the full output file path (character).

filter_top_pct <- function(transformed_cm_path, out_dir = "./", top_pct = 0.001, plot = FALSE) {
  suppressPackageStartupMessages({
    library(arrow)
  })

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  if (!is.numeric(top_pct) || length(top_pct) != 1 || top_pct <= 0 || top_pct >= 1) {
    stop("Error: 'top_pct' must be a single numeric value strictly between 0 and 1.")
  }

  cm <- read_feather(transformed_cm_path)
  if (!"pos" %in% colnames(cm)) {
    stop("Error: Expected a 'pos' column in the transformed count matrix, not found.")
  }

  pair_cols <- setdiff(colnames(cm), "pos")
  if (length(pair_cols) == 0) {
    stop("Error: No pair columns found besides 'pos'.")
  }

  # Per-pair threshold = value at the (1 - top_pct) quantile of that pair's column
  thresholds <- vapply(pair_cols, function(p) {
    quantile(cm[[p]], probs = 1 - top_pct, na.rm = TRUE, names = FALSE)
  }, numeric(1))
  names(thresholds) <- pair_cols

  # Union across pairs: a bin is kept if it clears ITS column's threshold in ANY pair
  keep_mat <- vapply(pair_cols, function(p) {
    cm[[p]] >= thresholds[p]
  }, logical(nrow(cm)))
  keep_idx <- which(rowSums(keep_mat) > 0)

  if (length(keep_idx) == 0) {
    stop("Error: No bins passed the top_pct threshold for any pair. Check 'top_pct' and input data.")
  }

  filtered_cm <- cm[keep_idx, , drop = FALSE]

  if (isTRUE(plot)) {
    n_selected_per_pair <- colSums(keep_mat)
    plot_df <- data.frame(
      pair = pair_cols,
      threshold = thresholds,
      n_selected = n_selected_per_pair
    )
    plot_path <- file.path(out_dir, "filter_top_pct_diagnostic.pdf")
    pdf(plot_path, width = 7, height = 5)
    op <- par(mar = c(10, 4, 4, 2))
    bar_x <- barplot(
      plot_df$n_selected,
      names.arg = NA,  # draw labels manually below, rotated at 45 degrees
      main = sprintf("Bins selected per pair (top %.2g%%)", top_pct * 100),
      ylab = "Number of bins above threshold"
    )
    text(
      x = bar_x, y = par("usr")[3] - 0.03 * diff(par("usr")[3:4]),
      labels = plot_df$pair, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 0.65
    )
    par(op)
    dev.off()
    cat("Saved diagnostic plot to: ", plot_path, "\n")
  }

  output_path <- file.path(out_dir, paste0("Filtered_TopPct_", format(top_pct, scientific = FALSE), ".feather"))
  write_feather(filtered_cm, output_path)
  cat("Successfully saved to: ", output_path, "\n")
  cat(sprintf("Selected %d / %d bins (union across %d pairs, top %.2g%% each)\n",
              length(keep_idx), nrow(cm), length(pair_cols), top_pct * 100))

  return(output_path)
}