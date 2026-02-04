# Qualification Control
#' Perform QC by counting reads/peaks from input BAM/BED files
#'
#' @param file_path Character vector of BAM or BED file paths
#' @param out_dir Directory to save output TSV files if save == TRUE
#' @param filtered_percentile Percentile threshold (0-1) for filtering (default: 0.25)
#' @param save Logical. If TRUE, save count tables as TSV
#'
#' @return List with:
#'   - all_df: Data frame of all files and their read/peak counts
#'   - filtered_df: Data frame after filtering by read count threshold
#'   - filtered_crf: Vector of file names after filtering
#'   - total_reads: Total read/peak count across all files
#'   
#' @export
qc <- function(file_path, out_dir = "./", filtered_percentile = 0.25, save = TRUE) {
  # Load libraries
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(ComplexHeatmap)
    library(glue)
    library(latex2exp)
    library(Rsamtools)
    library(GenomicAlignments)
    library(BiocParallel)
  })

  # Set up parallel processing
  BPPARAM <- get_BPPARAM()

  # Validate file formats
  file_exts <- tools::file_ext(file_path)
  if (all(file_exts == "bam")) {
    ext <- "bam"
  } else if (all(file_exts == "bed")) {
    ext <- "bed"
  } else {
    stop("Error: All files must be either BAM or BED format. Mixed formats are not supported.")
  }

  # Parallel processing of files
  counts <- BiocParallel::bplapply(file_path, function(path) {
    result <- data.frame()
    pair <- tools::file_path_sans_ext(basename(path))
    
    # Check file validity
    if (!(file.exists(path) && file.size(path) > 0)) {
      attr(result, "warning") <- sprintf("Skip invalid file: %s", path)
      return(result)
    }

    # Count reads/peaks with error handling
    tryCatch({
      if (ext == "bam") {
        temp <- GenomicAlignments::readGAlignmentPairs(path)
        read_count <- length(temp)
      } else {  # ext == "bed"
        peak <- ChIPseeker::readPeakFile(path, as = "GRanges")
        read_count <- length(peak)
      }
      
      return(data.frame(pair = pair, read_count = read_count, stringsAsFactors = FALSE))
      
    }, error = function(e) {
      attr(result, "warning") <- sprintf("Error processing %s: %s", path, e$message)
      return(result)
    })
  }, BPPARAM = BPPARAM)

  # Extract and display warnings
  warnings_msgs <- sapply(counts, function(x) attr(x, "warning"))
  warnings_msgs <- warnings_msgs[!sapply(warnings_msgs, is.null)]
  if (length(warnings_msgs) > 0) {
    for (msg in warnings_msgs) {
      warning(msg, call. = FALSE)
    }
  }

  # Filter out empty data.frames (invalid files)
  counts <- counts[sapply(counts, nrow) > 0]
  if (length(counts) == 0) {
    stop("Error: No valid files were successfully processed")
  }
  df <- as.data.frame(data.table::rbindlist(counts))
  
  # Calculate threshold and filter
  threshold <- quantile(df$read_count, probs = filtered_percentile, type = 3)
  filtered_df <- df[df$read_count >= threshold, ]
  filtered_df_vector <- filtered_df$pair
  total_reads <- sum(df$read_count)

  # Save results
  if (save) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    write.table(df, 
                file = file.path(out_dir, "all_read_count.tsv"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(filtered_df, 
                file = file.path(out_dir, "filtered_read_count.tsv"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(list(
    all_df = df, 
    filtered_df = filtered_df, 
    filtered_crf = filtered_df_vector, 
    total_reads = total_reads
  ))
}