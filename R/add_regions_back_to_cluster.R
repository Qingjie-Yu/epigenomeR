#' Add Non-highly-variable Regions Back to Clusters Based on Correlation
#'
#' Assigns cluster labels to regions excluded from the highly-variable set by
#' correlating them with existing cluster signatures.
#'
#' @param orig_cm_path Path to feather file from build_count_matrix
#' @param transformed_cm_path Path to feather file with transformed counts.
#' @param filtered_cm_path Path to feather file with highly-variable regions.
#' @param row_cluster_path Path to TSV with 'region' and 'cluster' columns.
#' @param out_dir Output directory where results and plots will be written.
#' @param cutoff_non_zero Min non-zero samples per region (default: 10).
#' @param quantile_threshold Quantile threshold for correlation filtering (0-1). Only regions with correlation above this quantile are assigned to clusters. (default: 0.75).
#' @param plot Save correlation histogram? (default: FALSE).
#'
#' @return Data frame with 'region' and 'cluster' columns. Labels follow priority:
#'         row_cluster_path > correlation-based > CRF_specific > Background.

add_regions_back_to_cluster <- function(orig_cm_path,
                                        transformed_cm_path,
                                        filtered_cm_path,
                                        row_cluster_path,
                                        out_dir,
                                        cutoff_non_zero = 10,
                                        quantile_threshold = 0.75,
                                        plot = FALSE) {
  suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(tibble)
  })

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Set up parallel processing
  BPPARAM <- get_BPPARAM(backend = "multicore")

  # Load data
  load_matrix <- function(path) {
    df <- read_feather(path)
    if ("pos" %in% colnames(df)) df <- tibble::column_to_rownames(df, "pos")
    as.matrix(df)
  }

  informative <- load_matrix(filtered_cm_path)
  transformed <- load_matrix(transformed_cm_path)

  clusters <- read.table(row_cluster_path, header = TRUE)
  if (!all(c("region", "cluster") %in% colnames(clusters))) {
    stop("Error: row_cluster_path must contain 'region' and 'cluster' columns")
  }
  cluster_vec <- setNames(clusters$cluster, clusters$region)

  orig <- read_feather(orig_cm_path)
  if ("pos" %in% colnames(orig)) orig <- tibble::column_to_rownames(orig, "pos")

  row_sums_orig <- rowSums(orig, na.rm = TRUE)
  nozero <- orig[row_sums_orig != 0, , drop = FALSE]
  if (nrow(nozero) == 0) {
    stop("All regions in orig_cm_path have zero counts after removing all-zero rows.")
  }

  # Filter non-informative regions
  noninformative <- nozero[setdiff(rownames(nozero), rownames(informative)), , drop = FALSE]
  nonzero_counts <- rowSums(noninformative > 0, na.rm = TRUE)
  filtered <- noninformative[nonzero_counts > cutoff_non_zero, , drop = FALSE]
  filtered_transformed <- transformed[rownames(transformed) %in% rownames(filtered), , drop = FALSE]

  # Compute cluster signatures (mean per cluster)
  common <- intersect(rownames(informative), names(cluster_vec))
  if (length(common) == 0) stop("No matching regions between informative and cluster assignments")

  cluster_means <- as.data.frame(informative[common, ]) |>
    mutate(cluster = cluster_vec[common]) |>
    group_by(cluster) |>
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
    tibble::column_to_rownames("cluster")

  # Correlate filtered regions with cluster signatures
  common_samples <- intersect(colnames(filtered_transformed), colnames(cluster_means))
  if (length(common_samples) == 0) stop("No common samples for correlation")

  region_names <- rownames(filtered_transformed)
  n_regions <- length(region_names)
  n_correlations <- n_regions * nrow(cluster_means)

  if (n_correlations < 1e7) {
    # Use serial for small datasets (< 10 million correlations)
    cor_mat <- cor(
      t(filtered_transformed[, common_samples, drop = FALSE]),
      t(cluster_means[, common_samples, drop = FALSE]),
      use = "complete.obs"
    )
  } else {
    # Parallel computation for large datasets
    filtered_t <- t(filtered_transformed[, common_samples, drop = FALSE])
    cluster_t <- t(cluster_means[, common_samples, drop = FALSE])

    cor_list <- BiocParallel::bplapply(seq_len(n_regions), function(i) {
      region_data <- filtered_t[, i]
      cors <- sapply(seq_len(ncol(cluster_t)), function(j) {
        cor(region_data, cluster_t[, j], use = "complete.obs")
      })
      return(cors)
    }, BPPARAM = BPPARAM)

    cor_mat <- matrix(
      unlist(cor_list, use.names = FALSE),
      nrow = n_regions,
      ncol = nrow(cluster_means),
      byrow = TRUE
    )
    rownames(cor_mat) <- region_names
    colnames(cor_mat) <- rownames(cluster_means)
  }

  # Find maximum correlation and best matching cluster for each region
  max_cor <- apply(cor_mat, 1, max, na.rm = TRUE)
  best_cluster <- colnames(cor_mat)[apply(cor_mat, 1, which.max)]

  # Filter by quantile threshold
  threshold <- quantile(max_cor, quantile_threshold, na.rm = TRUE)
  cat(sprintf("Correlation threshold (%.0f%%): %.4f\n", quantile_threshold * 100, threshold))

  if (plot) {
    histogram_path <- file.path(out_dir, "correlation_histogram.png")

    png(histogram_path, width = 800, height = 600)
    hist(max_cor,
      main = "Max Correlation Distribution", xlab = "Max Correlation",
      col = "lightblue", border = "black", breaks = 50
    )
    abline(v = threshold, col = "red", lty = 2, lwd = 2)
    legend("topright", sprintf("%.0f%% = %.3f", quantile_threshold * 100, threshold),
      col = "red", lty = 2
    )
    dev.off()
    cat("Saved histogram:", histogram_path, "\n")
  }

  high_cor_idx <- max_cor > threshold
  high_cor_regions <- data.frame(
    region = names(max_cor)[high_cor_idx],
    cor_mat[high_cor_idx, , drop = FALSE],
    max_correlation = max_cor[high_cor_idx],
    best_cluster = best_cluster[high_cor_idx],
    row.names = names(max_cor)[high_cor_idx]
  )

  cat(sprintf(
    "Regions: %d total -> %d high correlation (%.1f%%)\n",
    length(max_cor), sum(high_cor_idx), 100 * mean(high_cor_idx)
  ))

  # Build result with priority-based labels
  row_sums <- rowSums(transformed, na.rm = TRUE)
  valid_regions <- rownames(transformed)[row_sums != 0]

  result <- data.frame(region = valid_regions, cluster = "Background", stringsAsFactors = FALSE)
  result$cluster[result$region %in% rownames(filtered)] <- "CRF_specific"

  match_cor <- match(result$region, high_cor_regions$region)
  matched <- !is.na(match_cor)
  result$cluster[matched] <- high_cor_regions$best_cluster[match_cor[matched]]

  match_cluster <- match(result$region, clusters$region)
  matched <- !is.na(match_cluster)
  result$cluster[matched] <- as.character(clusters$cluster[match_cluster[matched]])

  # Save output
  tsv_all_path <- file.path(out_dir, "row_table_all.tsv")
  tsv_clean_path <- file.path(out_dir, "row_table_clean.tsv")

  write.table(result, tsv_all_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved full labels:", tsv_all_path, "\n")

  result_clean <- result %>% filter(!cluster %in% c("Background", "CRF_specific"))
  write.table(result_clean, tsv_clean_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved clean labels:", tsv_clean_path, "\n")
}
