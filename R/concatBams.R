# frag decomposition
# Post: Compute and save fragment length distributions from paired-end BAM files.
# Parameters: bamFiles: A character vector of paths to BAM files (must be paired-end and coordinate-sorted).
#             save_dir: Directory to save intermediate RData file and PDF histogram output.
# Output: Otherwise:
#             * Computes fragment lengths from BAM files using `readGAlignmentPairs`.
#             * Saves RData to: <save_dir>/<scen>-premerge_all-qc_frag_lens.RData
#             * Plots and saves fragment length histogram to: <scen>-premerge_all-qc_frag_hist.pdf
#             * Returns `bamWidths` (numeric vector of fragment lengths).
concatBams <- function(bamFiles, save_dir) {
  # load library
  suppressPackageStartupMessages({
    library(GenomicAlignments)
  })

  if (is.null(save_dir)) {
    stop("Error: Need output dir path: 'save_dir' is NULL.")
  }

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  bamWidths_dir_filename <- file.path(save_dir, "premerge_all-qc_frag_lens.RData")

  if (file.exists(bamWidths_dir_filename)) {
    load(bamWidths_dir_filename)
  } else {
    bamWidths <- c()
    for (bamFile in bamFiles) {
      temp <- readGAlignmentPairs(bamFile)
      locus <- data.frame(
        first_start = start(temp@first),
        first_end = end(temp@first),
        last_start = start(temp@last),
        last_end = end(temp@last)
      )
      start <- rowMin(as.matrix(locus))
      end <- rowMax(as.matrix(locus))
      strand <- "*"
      seqnames <- as.vector(seqnames(temp))
      bam <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames, strand = strand, start = start, end = end))
      bamWidths <- append(bamWidths, width(bam))
    }
    save(bamWidths, file = bamWidths_dir_filename)
    # write.csv(bamWidths, saveDirectory)
  }
  return(bamWidths)
}
