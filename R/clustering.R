# Parameters: cm_paths: character vector of feather file paths, one per sample
#             sample_names: character vector of sample names (same order as cm_paths)
#             row_km: number of row clusters for k-means
#             out_dir: root output directory; each CRF gets a subdirectory
#             crf_names: character vector of CRF names to analyze (default: NULL = all common CRFs)
#             seed: random seed for k-means reproducibility (default: 42)
#             apply_filter: whether to apply HVR filtering before clustering (default: TRUE)
#             order_clusters: whether to hierarchically order clusters by centroid similarity (default: TRUE)
#             cluster_linkage: linkage method for hierarchical ordering of clusters (default: "complete")
#             order_within_clusters: whether to reorder features within each cluster using hierarchical clustering (default: TRUE)
#             feature_distance: distance metric for within-cluster feature reordering (default: "euclidean")
#             feature_linkage: linkage method for within-cluster feature reordering (default: "complete")
#             plot: whether to generate heatmaps (default: TRUE)
#             apply_zscore: whether to z-score rows before plotting heatmap (default: TRUE).
#                           When TRUE, heatmap colors reflect relative differences across samples per region.
#                           When FALSE, heatmap colors reflect absolute normalized (libnorm + log2p1) values,
#                           which better reveals the magnitude of signal changes across samples.
#             show_column_names: whether to show sample names on heatmap (default: TRUE)
#            lower_range: lower bound for heatmap color scale (default: NULL, auto-calculated)
#            upper_range: upper bound for heatmap color scale (default: NULL, auto-calculated)

clustering <- function(cm_paths, sample_names, row_km, out_dir, crf_names = NULL, seed = 42, apply_filter = TRUE, order_clusters = TRUE, cluster_linkage = "complete", order_within_clusters = TRUE, feature_distance = "euclidean", feature_linkage = "complete", plot = TRUE, apply_zscore = TRUE, show_column_names = TRUE, lower_range = NULL, upper_range = NULL) {
  # peaks × CRFs
  sample_list <- lapply(seq_along(cm_paths), function(j) {
    df <- arrow::read_feather(cm_paths[j])
    pos <- df$pos
    df$pos <- NULL
    m <- as.matrix(df)
    mode(m) <- "numeric"
    rownames(m) <- as.character(pos)
    transform_mat(m, transformations = c("libnorm", "log2p1"))
  })

  if (is.null(crf_names)) {
    crf_names <- Reduce(intersect, lapply(sample_list, colnames))
  } else {
    # validate specified crf_names all exist
    common_crfs <- Reduce(intersect, lapply(sample_list, colnames))
    missing <- setdiff(crf_names, common_crfs)
    if (length(missing) > 0) {
      stop("The following CRFs are not found in all samples: ", paste(missing, collapse = ", "))
    }
  }
  sample_list <- lapply(sample_list, function(m) m[, crf_names, drop = FALSE])

  # peaks × samples
  crf_mats <- lapply(crf_names, function(crf) {
    m <- sapply(sample_list, function(s) s[, crf])
    colnames(m) <- sample_names
    m
  })
  names(crf_mats) <- crf_names

  # per-CRF clustering + heatmap
  row_paths <- list()

  for (crf in crf_names) {
    message("Processing CRF: ", crf)
    mat <- crf_mats[[crf]]
    crf_out_dir <- file.path(out_dir, crf)
    if (!dir.exists(crf_out_dir)) dir.create(crf_out_dir, recursive = TRUE, showWarnings = FALSE)

    # filtering
    if (apply_filter) {
      tmp_path <- file.path(crf_out_dir, paste0(crf, "_transformed.feather"))
      tmp_df <- as.data.frame(mat)
      tmp_df$pos <- rownames(mat)
      tmp_df <- tmp_df[, c("pos", setdiff(colnames(tmp_df), "pos"))]
      arrow::write_feather(tmp_df, tmp_path)

      filtered_path <- detect_hvr(transformed_cm_path = tmp_path, out_dir = crf_out_dir)
      filtered_df  <- arrow::read_feather(filtered_path)
      pos_filtered <- filtered_df$pos
      filtered_df$pos <- NULL
      mat <- as.matrix(filtered_df)
      mode(mat) <- "numeric"
      rownames(mat) <- pos_filtered
    } else {
      # save cleaned matrix
      tmp_df <- data.frame(pos = rownames(mat), as.data.frame(mat), check.names = FALSE)
      arrow::write_feather(tmp_df, file.path(crf_out_dir, paste0(crf, "_transformed.feather")))
    }

    # clustering
    result <- bidirectional_correlation_clustering(
      mat = mat, row_k = row_km, col_k = NULL, seed = seed, order_clusters = order_clusters, cluster_linkage = cluster_linkage, order_within_clusters = order_within_clusters, feature_distance = feature_distance, feature_linkage = feature_linkage
    )
    row_letter <- result$row_letter

    # optional zscore for heatmap
    if (apply_zscore) {
      mat_plot <- transform_mat(mat[names(row_letter), , drop = FALSE], transformations = "zscore")
      bad_z <- apply(mat_plot, 1, function(x) any(!is.finite(x)))
      if (any(bad_z)) {
        warning("Removing ", sum(bad_z), " zero-variance rows after z-scoring")
        mat_plot <- mat_plot[!bad_z, , drop = FALSE]
        row_letter <- row_letter[rownames(mat_plot)]
      }
      legend_title <- "Z-Score"
    } else {
      mat_plot <- mat[names(row_letter), , drop = FALSE]
      legend_title <- "log2(cpm)"
    }

    # save cluster assignments
    df_row <- data.frame(
      region  = names(row_letter),
      cluster = unname(row_letter),
      stringsAsFactors = FALSE
    )
    row_path <- file.path(crf_out_dir, paste0(crf, ".tsv"))
    write.table(df_row, row_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved cluster assignments to: ", row_path)
    row_paths[[crf]] <- row_path

    if (plot) {
      message("Generating heatmap for CRF: ", crf)
      clustering_heatmap(
        mat = mat_plot,
        row_cluster_file_path = row_path,
        out_dir = crf_out_dir,
        pdf_name = paste0(crf, ".pdf"),
        show_column_names = show_column_names,
        lower_range = lower_range,
        upper_range = upper_range,
        legend_title = legend_title
      )
    }
  }

  return(row_paths)
}