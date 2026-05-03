# Bidirectional K-means Clustering with Hierarchical Ordering
#
# Performs consensus k-means clustering on matrix rows and columns, then hierarchically
# orders the resulting cluster groups to create ordered cluster assignments similar to
# ComplexHeatmap package output.
#
# Parameters:
#   mat: Numeric matrix for clustering (must have row and column names)
#   row_k: Number of row clusters
#   col_k: Number of column clusters
#   row_repeats: Number of k-means repetitions for row consensus clustering (default: 1)
#   col_repeats: Number of k-means repetitions for column consensus clustering (default: 1)
#   seed: Random seed for reproducibility (default: 42)
#   cluster_method: Clustering method for rows and columns (default: "kmeans")
#                   "kmeans": consensus k-means with euclidean distance
#                   "correlation": hierarchical clustering with 1 - pearson correlation distance
#   do_order_clusters: Whether to hierarchically order clusters (default: TRUE)
#                      If FALSE, clusters are ordered by mean expression
#   do_reorder_within_clusters: Whether to reorder features within each cluster (default: TRUE)
#                               If FALSE, features maintain their original order within clusters
#   cluster_linkage: Linkage method for cluster-level hierarchical clustering
#                    (default: "complete"; options: "average", "single", "ward.D2", etc.)
#   feature_linkage: Linkage method for within-cluster feature reordering (default: NULL)
#                    If NULL, uses the same as cluster_linkage
#
# Returns:
#   List with two named vectors:
#     - row_letter: Named character vector with cluster labels (A, B, C, ...) for each row
#                   Names are row names, ordered by optimal display order
#     - col_num: Named integer vector with cluster numbers (1, 2, 3, ...) for each column
#                Names are column names, ordered by optimal display order

bidirectional_clustering <- function(mat, row_k, col_k = NULL, row_repeats = 1, col_repeats = 1, seed = 42, cluster_method = c("kmeans", "correlation"),do_order_clusters = TRUE, do_reorder_within_clusters = TRUE, cluster_linkage = "complete", feature_linkage = NULL) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix, not ", class(mat)[1], call. = FALSE)
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have row names and column names. Please set them before clustering.")
  }

  cluster_method <- match.arg(cluster_method)
  get_dist <- function(m) {
    if (cluster_method == "kmeans") {
      dist(m, method = "euclidean")
    } else {
      cor_mat <- cor(t(m), method = "pearson", use = "pairwise.complete.obs")
      cor_mat[is.na(cor_mat)] <- 0
      as.dist(1 - cor_mat)
    }
  }

  if (is.null(feature_linkage)) {
    feature_linkage <- cluster_linkage
  }

  consensus_kmeans <- function(m, k, reps) {
    parts <- lapply(seq_len(reps), function(i) {
      clue::as.cl_hard_partition(stats::kmeans(m, k, iter.max = 50))
    })
    cons <- clue::cl_consensus(clue::cl_ensemble(list = parts))
    as.vector(clue::cl_class_ids(cons))
  }

  # Order clusters hierarchically
  order_clusters <- function(m, cl, do_order) {
    unique_cl <- sort(unique(cl))
    row_mean  <- sapply(unique_cl, function(i) colMeans(m[cl == i, , drop = FALSE], na.rm = TRUE))
    if (!is.matrix(row_mean)) row_mean <- matrix(row_mean, nrow = 1)
    if (!do_order) {
      order <- order(colMeans(row_mean))
    } else {
      hc    <- hclust(get_dist(t(row_mean)), method = cluster_linkage)
      dend  <- reorder(as.dendrogram(hc), colMeans(row_mean), mean)
      order <- order.dendrogram(dend)
    }
    match(cl, unique_cl[order])
  }

  # Reorder within clusters
  reorder_within_clusters <- function(m, cl, do_reorder) {
    weights     <- -rowMeans(m, na.rm = TRUE)
    final_order <- integer(0)
    for (i in sort(unique(cl))) {
      idx <- which(cl == i)
      if (!do_reorder || length(idx) <= 1) {
        final_order <- c(final_order, idx)
      } else {
        submat <- m[idx, , drop = FALSE]
        hc     <- hclust(get_dist(submat), method = feature_linkage)
        dend   <- reorder(as.dendrogram(hc), weights[idx], mean)
        final_order <- c(final_order, idx[order.dendrogram(dend)])
      }
    }
    cl[final_order]
  }

  # row clustering
  set.seed(seed)
  if (row_k >= nrow(mat)) {
    row_cl        <- seq_len(nrow(mat))
    names(row_cl) <- rownames(mat)
  } else {
    if (cluster_method == "kmeans") {
      row_cl <- consensus_kmeans(mat, row_k, row_repeats)
    } else {
      d      <- get_dist(mat)
      hc     <- hclust(d, method = cluster_linkage)
      row_cl <- cutree(hc, k = row_k)
    }
    row_cl <- order_clusters(mat, row_cl, do_order_clusters)
    names(row_cl) <- rownames(mat)
    row_cl <- reorder_within_clusters(mat, row_cl, do_reorder_within_clusters)
  }

  row_letter        <- LETTERS[row_cl]
  names(row_letter) <- names(row_cl)

  # col cluster
  if (is.null(col_k) || col_k >= ncol(mat)) {
    col_cl        <- seq_len(ncol(mat))
    names(col_cl) <- colnames(mat)
  } else {
    set.seed(seed + 1)
    if (cluster_method == "kmeans") {
      col_cl <- consensus_kmeans(t(mat), col_k, col_repeats)
    } else {
      d      <- get_dist(t(mat))
      hc     <- hclust(d, method = cluster_linkage)
      col_cl <- cutree(hc, k = col_k)
    }
    col_cl <- order_clusters(t(mat), col_cl, do_order_clusters)
    names(col_cl) <- colnames(mat)
    col_cl <- reorder_within_clusters(t(mat), col_cl, do_reorder_within_clusters)
  }
  col_num <- col_cl
  names(col_num) <- names(col_cl)
  
  invisible(list(row_letter = row_letter, col_num = col_num))
}
