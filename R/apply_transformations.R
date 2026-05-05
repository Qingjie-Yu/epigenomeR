transform_mat <- function(mat, transformations = c("libnorm", "log2p1")) {
  for (t in transformations) {
    if (t == "remove0") {
      zero_rows <- rowSums(mat) == 0
      n_zero <- sum(zero_rows)
      message(if (n_zero > 0) paste("Removed", n_zero, "all-zero row(s)") else "No all-zero rows found")
      if (n_zero > 0) mat <- mat[!zero_rows, , drop = FALSE]
      if (nrow(mat) == 0) stop("Error: all rows have zero counts; nothing left after filtering.")
    } else if (t == "libnorm") {
      col_sums <- colSums(mat, na.rm = TRUE)
      if (any(col_sums == 0)) stop("Error: one or more samples have zero total counts.")
      mat <- t(t(mat) / col_sums * 1e6)
    } else if (t == "log2p1") {
      mat <- log2(mat + 1)
    } else if (t == "qnorm") {
      rn <- rownames(mat); cn <- colnames(mat)
      mat <- preprocessCore::normalize.quantiles(mat, copy = TRUE)
      rownames(mat) <- rn; colnames(mat) <- cn
    } else if (t == "zscore") {
      mat <- t(scale(t(mat)))
    } else if (t == "minmaxnorm") {
      col_mins   <- apply(mat, 2, min)
      col_ranges <- apply(mat, 2, max) - col_mins
      is_constant <- col_ranges == 0
      mat <- sweep(mat, 2, col_mins, FUN = "-")
      if (any(is_constant))  { mat[, is_constant] <- 0; message(sum(is_constant), " constant column(s) set to 0") }
      if (any(!is_constant))   mat[, !is_constant] <- sweep(mat[, !is_constant, drop = FALSE], 2, col_ranges[!is_constant], FUN = "/")
    } else if (t == "sqrt") {
      mat <- sqrt(mat)
    } else {
      warning("Unrecognized transformation: '", t, "'. Skipping!")
    }
  }
  mat
}

# Apply Transformation
# Post: Apply a series of transformations to count matrix data.
# Supported transformations:
#   - "remove0": Remove all-zero rows
#   - "libnorm": Library size normalization (CPM-style)
#   - "log2p1": Log2(x+1) transformation
#   - "sqrt": Square root transformation
#   - "minmaxnorm": Min-max normalization to [0,1] range
#   - "qnorm": Quantile normalization across samples
#
# Parameters:
#   cm_path: Path to input count matrix (.feather file)
#   transformations: Character vector of transformation steps to apply in order. Default: c("remove0", "libnorm", "log2p1", "qnorm")
#   out_dir: Directory path for saving output files
#   save_each_step: Logical. If TRUE, save intermediate results after each transformation. Default: FALSE
#
# Output: Saves transformed count matrix as .feather file with "_transformed" suffix.  Returns the full output file path (character).

apply_transformations <- function(cm_path, out_dir = "./", transformations = c("libnorm", "log2p1"), save_each_step = FALSE) {
  suppressPackageStartupMessages({
    library(arrow); library(tibble)
    library(matrixStats); library(preprocessCore)
  })

  df <- as.data.frame(read_feather(cm_path), stringsAsFactors = FALSE, check.names = FALSE)
  input_prefix <- basename(tools::file_path_sans_ext(cm_path))

  pos_colname <- if ("pos" %in% colnames(df)) "pos" else {
    warning("'pos' column not found. Using first column as coordinate.")
    colnames(df)[1]
  }
  rownames(df) <- df[[pos_colname]]; df[[pos_colname]] <- NULL

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  SHARED_TRANSFORMS <- c("remove0", "libnorm", "log2p1", "qnorm", "zscore", "minmaxnorm", "sqrt")

  applied <- c()
  for (t in transformations) {
    if (t %in% SHARED_TRANSFORMS) {
      df <- as.data.frame(transform_mat(as.matrix(df), t))
    } else {
      warning("Unrecognized transformation: '", t, "'. Skipping!")
      next
    }

    applied <- c(applied, t)
    if (isTRUE(save_each_step)) {
      out <- rownames_to_column(as.data.frame(df), var = pos_colname)
      write_feather(out, file.path(out_dir, paste0(input_prefix, "_", paste(applied, collapse = "_"), ".feather")))
    }
  }

  out <- rownames_to_column(as.data.frame(df), var = pos_colname)
  output_path <- file.path(out_dir, paste0(input_prefix, "_transformed.feather"))
  write_feather(out, output_path)
  output_path
}