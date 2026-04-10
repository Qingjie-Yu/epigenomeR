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
#   do_order_clusters: Whether to hierarchically order clusters (default: TRUE)
#                      If FALSE, clusters are ordered by mean expression
#   cluster_distance: Distance metric for cluster-level hierarchical clustering
#                     (default: "euclidean"; options: "manhattan", "maximum", etc.)
#   cluster_linkage: Linkage method for cluster-level hierarchical clustering
#                    (default: "complete"; options: "average", "single", "ward.D2", etc.)
#   do_reorder_within_clusters: Whether to reorder features within each cluster (default: TRUE)
#                               If FALSE, features maintain their original order within clusters
#   feature_distance: Distance metric for within-cluster feature reordering (default: NULL)
#                     If NULL, uses the same as cluster_distance
#   feature_linkage: Linkage method for within-cluster feature reordering (default: NULL)
#                    If NULL, uses the same as cluster_linkage
#
# Returns:
#   List with two named vectors:
#     - row_letter: Named character vector with cluster labels (A, B, C, ...) for each row
#                   Names are row names, ordered by optimal display order
#     - col_num: Named integer vector with cluster numbers (1, 2, 3, ...) for each column
#                Names are column names, ordered by optimal display order

bidirectional_kmeans_clustering <- function(mat, row_k, col_k, row_repeats = 1, col_repeats = 1, seed = 42, do_order_clusters = TRUE, cluster_distance = "euclidean", cluster_linkage = "complete", do_reorder_within_clusters = TRUE, feature_distance = NULL, feature_linkage = NULL) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix, not ", class(mat)[1], call. = FALSE)
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have row names and column names. Please set them before clustering.")
  }

  if (is.null(feature_distance)) {
    feature_distance <- cluster_distance
  }
  if (is.null(feature_linkage)) {
    feature_linkage <- cluster_linkage
  }

  # Consensus k-means
  consensus_kmeans <- function(m, k, reps) {
    if (k >= nrow(m)) return(seq_len(nrow(m)))
    parts <- lapply(seq_len(reps), function(i) {
      clue::as.cl_hard_partition(stats::kmeans(m, k, iter.max = 50))
    })
    cons <- clue::cl_consensus(clue::cl_ensemble(list = parts))
    as.vector(clue::cl_class_ids(cons))
  }

  # Order clusters hierarchically
  order_clusters <- function(m, cl, do_order) {
    unique_cl <- sort(unique(cl))
    row_mean <- sapply(unique_cl, function(i) {
      colMeans(m[cl == i, , drop = FALSE], na.rm = TRUE)
    })
    if (!is.matrix(row_mean)) row_mean <- matrix(row_mean, nrow = 1)
    if (!do_order) {
      order <- order(colMeans(row_mean))
    } else {
      hc <- hclust(dist(t(row_mean), method = cluster_distance), method = cluster_linkage)
      weights <- colMeans(row_mean)
      dend <- reorder(as.dendrogram(hc), weights, mean)
      order <- order.dendrogram(dend)
    }
    new_cl <- match(cl, unique_cl[order])
    new_cl
  }

  # Reorder within clusters
  reorder_within_clusters <- function(m, cl, do_reorder) {
    weights <- -rowMeans(m, na.rm = TRUE)
    final_order <- integer(0)
    if (!do_reorder) {
      for (i in sort(unique(cl))) {
        idx <- which(cl == i)
        final_order <- c(final_order, idx)
      }
    } else {
      for (i in sort(unique(cl))) {
        idx <- which(cl == i)
        if (length(idx) <= 1) {
          final_order <- c(final_order, idx)
        } else {
          submat <- m[idx, , drop = FALSE]
          hc <- hclust(dist(submat, method = feature_distance), method = feature_linkage)
          dend <- reorder(as.dendrogram(hc), weights[idx], mean)
          final_order <- c(final_order, idx[order.dendrogram(dend)])
        }
      }
    }
    cl <- cl[final_order]
    cl
  }

  # row cluster
  set.seed(seed)
  row_cl <- consensus_kmeans(mat, row_k, row_repeats)
  row_cl <- order_clusters(mat, row_cl, do_order_clusters)
  names(row_cl) <- rownames(mat)
  row_cl <- reorder_within_clusters(mat, row_cl, do_reorder_within_clusters)
  row_letter <- LETTERS[row_cl]
  names(row_letter) <- names(row_cl)

  # col cluster
  set.seed(seed + 1)
  col_cl <- consensus_kmeans(t(mat), col_k, col_repeats)
  col_cl <- order_clusters(t(mat), col_cl, do_order_clusters)
  names(col_cl) <- colnames(mat)
  col_cl <- reorder_within_clusters(t(mat), col_cl, do_reorder_within_clusters)

  list(row_letter = row_letter, col_num = col_cl)
}
