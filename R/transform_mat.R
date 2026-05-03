transform_mat <- function(mat, transformations = c("libnorm", "log2p1")
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

