transform_mat <- function(mat, transformations = c("log2p1", "zscore")
) {
  for (t in transformations) {
    if (t == "libnorm") {
      mat <- t(t(mat) / colSums(mat) * 1e6)
    } else if (t == "log2p1") {
      mat <- log2(mat + 1)
    } else if (t == "qnorm") {
      rn  <- rownames(mat)
      cn  <- colnames(mat)
      mat <- preprocessCore::normalize.quantiles(mat)
      rownames(mat) <- rn
      colnames(mat) <- cn
    } else if (t == "zscore") {
      mat <- t(scale(t(mat)))
    } else {
      warning("Unrecognized transformation: '", t, "'. Skipping!")
    }
  }
  mat
}

# Parameters: cm_paths: character vector of feather file paths, one per sample
#             sample_names: character vector of sample names (same order as cm_paths)
#             row_km: number of row clusters for k-means
#             out_dir: root output directory; each CRF gets a subdirectory
#             crf_names: character vector of CRF names to analyze (default: NULL = all common CRFs)
#             seed: random seed for k-means reproducibility (default: 42)
#             plot: whether to generate heatmaps (default: TRUE)
#             show_column_names: whether to show sample names on heatmap (default: TRUE)
#             apply_filter: whether to apply HVR filtering before clustering (default: TRUE)
#             apply_qnorm: whether to apply quantile normalization per CRF (default: FALSE)

clustering_single_crf <- function(cm_paths, sample_names, row_km, out_dir, crf_names = NULL, seed = 42, plot = TRUE, show_column_names = TRUE, apply_filter = TRUE, apply_qnorm = FALSE) {
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

  if (apply_qnorm) {
    crf_mats <- lapply(crf_mats, function(m) {
      transform_mat(m, transformations = c("qnorm"))
    })
    names(crf_mats) <- crf_names
  }

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
      bad_rows <- apply(mat, 1, function(x) any(!is.finite(x)))
      if (any(bad_rows)) {
        message("Removing ", sum(bad_rows), " peaks with non-finite values")
        mat <- mat[!bad_rows, , drop = FALSE]
      }
      # save cleaned matrix
      tmp_df <- data.frame(pos = rownames(mat), as.data.frame(mat), check.names = FALSE)
      arrow::write_feather(tmp_df, file.path(crf_out_dir, paste0(crf, "_transformed.feather")))
    }

    # z-score
    mat <- transform_mat(mat, transformations = "zscore")

    # clustering
    result <- bidirectional_kmeans_clustering(
      mat = mat, row_k = row_km, col_k = NULL, seed = seed
    )
    row_letter <- result$row_letter

    # save cluster assignments
    df_row <- data.frame(
      region  = names(row_letter),
      cluster = unname(row_letter),
      stringsAsFactors = FALSE
    )
    row_path <- file.path(out_dir, paste0(crf, ".tsv"))
    write.table(df_row, row_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved cluster assignments to: ", row_path)
    row_paths[[crf]] <- row_path

    # heatmap
    if (plot) {
      message("Generating heatmap for CRF: ", crf)
      clustering_heatmap(
        mat                   = mat[names(row_letter), , drop = FALSE],
        row_cluster_file_path = row_path,
        out_dir               = crf_out_dir,
        pdf_name              = paste0(crf, ".pdf"),
        show_column_names     = show_column_names
      )
    }
  }

  return(row_paths)
}