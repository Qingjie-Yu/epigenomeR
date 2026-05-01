# Automatically performs k-means clustering, generates cluster files internally, and draw heatmap
# Description: This function reads a count matrix from a `.feather` file, performs
#              k-means clustering on both rows (genomic features) and columns (samples),
#              saves the cluster assignments, and optionally generates a publication-ready
#              heatmap in a single workflow.
# Parameters:  cm_path: Path to count matrix `.feather` file
#                        (must contain a `pos` column for feature IDs and sample columns).
#              row_km: Number of k-means clusters for rows (genomic features).
#              col_km: Number of k-means clusters for columns (samples).
#              out_dir: Directory to save cluster tables and heatmap output.
#              seed: Random seed for reproducible clustering (default: 42).
#              plot: Whether to generate a heatmap plot (default: TRUE).
#              show_column_names: Whether to show column names at the bottom of the
#                                 heatmap (default: FALSE).
#              lower_range: Lower bound for the heatmap color scale
#                          (default: NULL, automatically determined).
#              upper_range: Upper bound for the heatmap color scale
#                          (default: NULL, automatically determined).
#              row_title_fontsize: Font size for row cluster titles (default: NULL).
#              col_title_fontsize: Font size for column cluster titles (default: NULL).
#              legend_title_fontsize: Font size for the legend title (default: NULL).
#              legend_label_fontsize: Font size for the legend labels (default: NULL).
# Output:     A named list with:
#              - "row_table": Path to the saved row cluster assignment table (`row_table.tsv`).
#              - "col_table": Path to the saved column cluster assignment table (`col_table.tsv`).
#             (Cluster tables and heatmap files are written to `out_dir`.)

biclustering <- function(cm_path, row_km, col_km, out_dir, seed = 42, plot = TRUE, show_column_names = FALSE, lower_range = NULL, upper_range = NULL, row_title_fontsize = NULL, col_title_fontsize = NULL, legend_title_fontsize = NULL, legend_label_fontsize = NULL) {
  library(arrow)
  library(tibble)
  library(dplyr)

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
  result <- bidirectional_kmeans_clustering(mat = mat, row_k = row_km, col_k = col_km, seed = seed)
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

biclustering_by_crf <- function(cm_paths, sample_names, row_km, out_dir, crf_names = NULL, seed = 42, plot = TRUE, show_column_names = TRUE) {
  # peaks × CRFs
  sample_list <- lapply(seq_along(cm_paths), function(j) {
    df <- arrow::read_feather(cm_paths[j])
    pos <- df$pos
    df$pos <- NULL
    m <- as.matrix(df)
    mode(m) <- "numeric"
    rownames(m) <- as.character(pos)
    m
  })

  if (is.null(crf_names)) {
    common_crfs <- Reduce(intersect, lapply(sample_list, colnames))
    crf_names <- sort(common_crfs)
  }
  sample_list <- lapply(sample_list, function(m) m[, crf_names, drop = FALSE])

  # peaks × samples
  crf_mats <- lapply(crf_names, function(crf) {
    m <- sapply(sample_list, function(s) s[, crf])
    colnames(m) <- sample_names
    m
  })
  names(crf_mats) <- crf_names

  # per-(sample,CRF) Z-score
  crf_mats_z_orig <- lapply(crf_mats, function(m) {
    scale(m)
  })
  names(crf_mats_z_orig) <- crf_names
  crf_mats_z <- lapply(seq_along(crf_mats_z_orig), function(k) {
    m <- crf_mats_z_orig[[k]]
    colnames(m) <- paste0(crf_names[k], "_", colnames(m))
    m
  })
  mat_concat <- do.call(cbind, crf_mats_z)

  # clustering
  result <- bidirectional_kmeans_clustering(
    mat = mat_concat, row_k = row_km, col_k = NULL, seed = seed
  )
  row_letter <- result$row_letter

  df_row <- data.frame(
    region = names(row_letter),
    cluster = unname(row_letter),
    stringsAsFactors = FALSE
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  row_path <- file.path(out_dir, "row_table.tsv")
  write.table(df_row, row_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved cluster assignments to:", row_path)

  # multi-heatmap
  if (plot) {
    message("Generating heatmap...")
    # restore original (non-prefixed) crf_mats_z for plotting
    crf_mats_plot <- lapply(crf_mats_z_orig, function(m) {
      m[names(row_letter), , drop = FALSE]
    })
    biclustering_multi_heatmap(
      crf_mat_list          = crf_mats_plot,
      row_cluster_file_path = row_path,
      out_dir               = out_dir,
      show_column_names     = show_column_names
    )
  }

  return(row_path)
}