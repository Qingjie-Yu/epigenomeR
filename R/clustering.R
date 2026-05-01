transform_crf_mat <- function(
  mat,          # peaks × samples，已经是 log2p1 值
  do_log2  = TRUE,
  do_qnorm = TRUE,
  do_zscore = TRUE
) {
  # Step 1: log2 
  if (do_log2) {
    mat <- log2(mat + 1)
  }

  # Step 2: quantile normalization across samples（列方向）
  if (do_qnorm) {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("请安装 preprocessCore: BiocManager::install('preprocessCore')")
    }
    rn <- rownames(mat)
    cn <- colnames(mat)
    mat <- preprocessCore::normalize.quantiles(mat)
    rownames(mat) <- rn
    colnames(mat) <- cn
  }

  # Step 3: row-wise z-score（每个 peak 跨样本标准化，用于画图）
  if (do_zscore) {
    mat <- t(scale(t(mat)))
  }

  mat
}

# Per-CRF Independent Biclustering
# Post: For each specified CRF, independently perform row k-means clustering
#       on the peaks x samples z-score matrix, save cluster assignments,
#       and optionally generate a heatmap.
# Parameters: cm_paths: character vector of feather file paths, one per sample
#             sample_names: character vector of sample names (same order as cm_paths)
#             row_km: number of row clusters for k-means
#             out_dir: root output directory; each CRF gets a subdirectory
#             crf_names: character vector of CRF names to analyze (default: NULL = all common CRFs)
#             seed: random seed for k-means reproducibility (default: 42)
#             plot: whether to generate heatmaps (default: TRUE)
#             show_column_names: whether to show sample names on heatmap (default: TRUE)
#             tsv_name: filename for cluster assignment .tsv (default: "row_clusters.tsv")
# Output: named list of row cluster .tsv file paths, one per CRF

clustering_single_crf <- function(cm_paths, sample_names, row_km, out_dir, crf_names = NULL, seed = 42, plot = TRUE, show_column_names = TRUE) {
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
    crf_names <- Reduce(intersect, lapply(sample_list, colnames))
  } else {
    # validate specified crf_names all exist
    common_crfs <- Reduce(intersect, lapply(sample_list, colnames))
    missing <- setdiff(crf_names, common_crfs)
    if (length(missing) > 0) {
      stop("The following CRFs are not found in all samples: ", paste(missing, collapse = ", "))
    }
  }

  # peaks × samples
  crf_mats <- lapply(crf_names, function(crf) {
    m <- sapply(sample_list, function(s) s[, crf])
    colnames(m) <- sample_names
    m
  })
  names(crf_mats) <- crf_names

  # Transform：log2 → qnorm → row z-score
  crf_mats_transformed <- lapply(crf_mats, function(m) {
    transform_crf_mat(m, do_log2 = TRUE, do_qnorm = TRUE, do_zscore = TRUE)
  })
  names(crf_mats_transformed) <- crf_names
  
  bad_rows <- Reduce("|", lapply(crf_mats_transformed, function(m) {
    apply(m, 1, function(x) any(!is.finite(x)))
  }))
  if (any(bad_rows)) {
    message("Removing ", sum(bad_rows), " peaks with non-finite values across CRFs")
    crf_mats_transformed <- lapply(crf_mats_transformed, function(m) m[!bad_rows, , drop = FALSE])
  }

  # per-CRF clustering + heatmap
  row_paths <- list()

  for (crf in crf_names) {
    message("Processing CRF: ", crf)
    mat <- crf_mats_transformed[[crf]]
    crf_out_dir <- file.path(out_dir, crf)
    if (!dir.exists(crf_out_dir)) dir.create(crf_out_dir, recursive = TRUE, showWarnings = FALSE)

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
      biclustering_heatmap(
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