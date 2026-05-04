# Automatically performs clustering, generates cluster files internally, and draws heatmap
# Description: This function reads a count matrix from a `.feather` file, performs
#              bidirectional clustering on both rows (genomic features) and columns (samples),
#              saves the cluster assignments, and optionally generates a publication-ready
#              heatmap in a single workflow. Supports two clustering methods: consensus k-means
#              with euclidean distance, or hierarchical clustering with 1 - Pearson correlation
#              distance.
# Parameters:  cm_path: Path to count matrix `.feather` file
#                        (must contain a `pos` column for feature IDs and sample columns).
#              row_km: Number of clusters for rows (genomic features).
#              col_km: Number of clusters for columns (samples).
#              out_dir: Directory to save cluster tables and heatmap output.
#              seed: Random seed for reproducible clustering (default: 42).
#              cluster_method: Clustering method for both rows and columns (default: "kmeans").
#                              "kmeans": consensus k-means with euclidean distance.
#                              "correlation": hierarchical clustering with 1 - Pearson
#                                             correlation distance, using average linkage.
#              order_clusters: Whether to hierarchically order clusters by centroid similarity
#                              (default: TRUE). If FALSE, clusters are ordered by mean expression.
#              cluster_linkage: Linkage method for hierarchical ordering of clusters (default: "complete").
#              order_within_clusters: Whether to reorder features within each cluster using
#                                     hierarchical clustering (default: TRUE). If FALSE,
#                                     features maintain their original order within clusters.
#              feature_distance: Distance metric for within-cluster feature reordering
#                                (default: "euclidean"). Any method accepted by dist() is valid.
#              feature_linkage: Linkage method for within-cluster feature reordering
#                               (default: "complete"). Any method accepted by hclust() is valid.
#              plot: Whether to generate a heatmap plot (default: TRUE).
#              show_column_names: Whether to show column names at the bottom of the
#                                 heatmap (default: FALSE).
#              lower_range: Lower bound for the heatmap color scale
#                           (default: NULL, automatically determined).
#              upper_range: Upper bound for the heatmap color scale
#                           (default: NULL, automatically determined).
#              row_title_fontsize: Font size for row cluster titles (default: NULL).
#              col_title_fontsize: Font size for column cluster titles (default: NULL).
#              legend_title_fontsize: Font size for the legend title (default: NULL).
#              legend_label_fontsize: Font size for the legend labels (default: NULL).
# Output:     A named list with:
#              - "row_table": Path to the saved row cluster assignment table (`row_table.tsv`).
#              - "col_table": Path to the saved column cluster assignment table (`col_table.tsv`).
#             (Cluster tables and heatmap files are written to `out_dir`.)

biclustering <- function(cm_path, row_km, col_km, out_dir, seed = 42, order_clusters = TRUE, cluster_linkage = "complete", order_within_clusters = TRUE, feature_distance = "euclidean", feature_linkage = "complete", plot = TRUE, show_column_names = FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  # Load Library
  suppressPackageStartupMessages({
    library(tibble)
    library(arrow)
    library(dplyr)
  })

  cluster_method <- match.arg(cluster_method)

  if (is.null(cm_path) || length(cm_path) == 0) {
    stop("`count_matrix_file_path` is required", call. = FALSE)
  }
  df <- arrow::read_feather(cm_path)
  pos <- df$pos
  df$pos <- NULL
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  rownames(mat) <- as.character(pos)

  mat <- mat[, sort(colnames(mat))]
  col_km <- min(col_km, ncol(mat))

  message("Performing bidirectional k-means clustering...")
  result <- bidirectional_kmeans_clustering(
    mat = mat, 
    row_k = row_km, col_k = col_km, 
    seed = seed, 
    order_clusters = order_clusters, cluster_linkage = cluster_linkage, 
    order_within_clusters = order_within_clusters, feature_distance = feature_distance, feature_linkage = feature_linkage
  )
  row_letter <- result$row_letter
  col_num <- result$col_num

  df_row <- data.frame(
    region = names(row_letter),
    cluster = unname(row_letter),
    stringsAsFactors = FALSE
  )

  df_col <- data.frame(
    pair = names(col_num),
    cluster = unname(col_num),
    stringsAsFactors = FALSE
  )

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path1 <- file.path(out_dir, "row_table.tsv")
  path2 <- file.path(out_dir, "col_table.tsv")
  write.table(df_row, path1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(df_col, path2, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved cluster assignments to:")
  message("  - row clusters: ", path1)
  message("  - col clusters: ", path2)

  if (plot) {
    message("Generating heatmap...")
    mat_sorted <- mat[names(result$row_letter), names(result$col_num)]
    biclustering_heatmap(mat = mat_sorted, row_cluster_file_path = path1, col_cluster_file_path = path2, out_dir = out_dir, show_column_names = show_column_names, lower_range = lower_range, upper_range = upper_range, row_title_fontsize = row_title_fontsize, col_title_fontsize = col_title_fontsize, legend_title_fontsize = legend_title_fontsize, legend_label_fontsize = legend_label_fontsize)
  }

  return(list("row_table" = path1, "col_table" = path2))
}