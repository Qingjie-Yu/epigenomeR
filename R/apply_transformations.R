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

apply_transformations <- function(cm_path, out_dir = "./", transformations = c("remove0", "libnorm", "log2p1", "qnorm"), save_each_step = FALSE) {
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(matrixStats)
    library(preprocessCore)
  })

  df <- read_feather(cm_path)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  input_prefix <- basename(tools::file_path_sans_ext(cm_path))

  pos_colname <- "pos"
  if (pos_colname %in% colnames(df)) {
    rownames(df) <- df[[pos_colname]]
    df[[pos_colname]] <- NULL
  } else {
    warning(pos_colname, " column not found in the feather file. Automatically set the first column as the chromosome coordinate information.")
    first_col <- colnames(df)[1]
    rownames(df) <- df[[first_col]]
    df[[first_col]] <- NULL
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  applied <- c()
  for (t in transformations) {
    if (t == "remove0") {
      zero_rows <- rowSums(df) == 0
      n_zero <- sum(zero_rows)
      if (n_zero > 0) {
        df <- df[!zero_rows, , drop = FALSE]
        message("Removed ", n_zero, " all-zero row(s)")
      } else {
        message("No all-zero rows found")
      }
      if (nrow(df) == 0) {
        stop("Error: all rows have zero counts; nothing left after filtering.")
      }
    } else if (t == "libnorm") {
      message("\n=== Perform Library Normalization ===")
      total_counts <- sum(df, na.rm = TRUE)
      if (total_counts == 0) {
        stop("Error: total counts across all targets is zero.")
      }
      df <- df / total_counts * 1e6
      message("Library normalization completed")
    } else if (t == "log2p1") {
      message("\n=== Apply log2(x+1) Transformation ===")
      df <- log2(df + 1)
      message("log2(x+1) transformation completed")
    } else if (t == "minmaxnorm") {
      message("\n=== Apply Min-Max Normalization ===")
      col_mins <- apply(df, 2, min)
      col_ranges <- apply(df, 2, max) - col_mins
      is_constant <- col_ranges == 0
      df <- sweep(df, 2, col_mins, FUN = "-")

      if (any(is_constant)) {
        df[, is_constant] <- 0
        message(sum(is_constant), " constant column(s) normalized to 0")
      }
      if (any(!is_constant)) {
        df[, !is_constant] <- sweep(df[, !is_constant, drop = FALSE], 2, col_ranges[!is_constant], FUN = "/")
      }
      message("Min-Max normalization completed")
    } else if (t == "sqrt") {
      message("\n=== Square Root ===")
      df <- sqrt(df)
    } else if (t == "qnorm") {
      message("\n=== Quantile Normalization ===")
      old_rownames <- rownames(df)
      old_colnames <- colnames(df)
      df <- normalize.quantiles(as.matrix(df), copy = TRUE)
      df <- as.data.frame(df)
      rownames(df) <- old_rownames
      colnames(df) <- old_colnames
      message("Quantile normalization completed")
    } else {
      warning("Unrecognized transformation: '", t, "'. Skipping!")
      next
    }

    applied <- c(applied, t)
    if (save_each_step == TRUE) {
      df_to_save <- rownames_to_column(as.data.frame(df), var = pos_colname)
      file_name <- paste0(input_prefix, "_", paste(applied, collapse = "_"), ".feather")
      save_path <- file.path(out_dir, file_name)
      write_feather(df_to_save, save_path)
    }
  }

  df_to_save <- rownames_to_column(as.data.frame(df), var = pos_colname)
  output_name <- paste0(input_prefix, "_transformed.feather")
  output_path <- file.path(out_dir, output_name)
  write_feather(df_to_save, output_path)
  return(output_path)
}
