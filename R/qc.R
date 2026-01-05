# Qualification Control
#' Perform QC by counting reads/peaks from input BAM/BED files
#'
#' @param file_path Character vector of BAM or BED file paths
#' @param filtered_percentile Percentile threshold (0-1) for filtering (default: 0.25)
#' @param out_dir Directory to save output CSVs if save == TRUE
#' @param save Logical. If TRUE, save count tables as CSV
#'
#' @return List with:
#'   - all_df: Data frame of all files and their read/peak counts
#'   - filtered_df: Data frame after filtering by read count threshold
#'   - filtered_crf: Vector of file names after filtering
#'   - total_reads: Total read/peak count across all files
qc <- function(file_path, filtered_percentile = 0.25, out_dir = "./", save = TRUE) {
  # load library
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(ComplexHeatmap)
    library(glue)
    library(latex2exp)
    library(Rsamtools)
    library(GenomicAlignments)
  })

  file_exts <- tools::file_ext(file_path)
  if (all(file_exts == "bam")) {
    ext <- "bam"
  } else if (all(file_exts == "bed")) {
    ext <- "bed"
  } else {
    stop("Error: file formats must be either 'bam' or 'bed'.")
  }

  df <- data.frame(pair = character(), read_count = numeric(), stringsAsFactors = FALSE)

  for (k in seq_along(file_path)) {
    file_path <- file_path[k]
    file_name <- tools::file_path_sans_ext(basename(file_path))
    if (!(file.exists(file_path) && file.size(file_path) > 0)) {
      warning(glue::glue("Skip invalid file: {file_path}"))
      next
    }

    if (ext == "bam") {
      temp <- readGAlignmentPairs(file_path)
      count <- length(temp)
    } else if (ext == "bed") {
      peak <- ChIPseeker::readPeakFile(file_path, as = "GRanges")   # Use ChIPseeker to read peak files
      count <- length(peak)
    }
    df <- rbind(df, data.frame(pair = file_name, read_count = count, stringsAsFactors = FALSE)) # return: total df
  }

  threshold <- quantile(as.numeric(df$read_count), probs = filtered_percentile, type = 3)
  filtered_df <- df[df$read_count >= threshold, ]
  all_df <- df
  filtered_df_vector <- filtered_df$pair # return: filtered vector
  total_reads <- sum(df$read_count) # return: total reads

  # save
  if (save == TRUE) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    write.table(all_df, file = file.path(out_dir, "all_read_count.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(filtered_df, file = file.path(out_dir, "filtered_read_count.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(list(all_df = all_df, filtered_df = filtered_df, filtered_crf = filtered_df_vector, total_reads = total_reads))
}