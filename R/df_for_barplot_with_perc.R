# Post: Calculate fragment length distribution percentages across nucleosome categories (subnucleosomal, mononucleosomal, dinucleosomal) for multiple BAM files and save results as barplot-ready data frame.
# Parameter: bam_file_path: Vector of BAM file paths to analyze
#            save_dir: Directory path to save the output table
#            valley1: First valley position separating subnucleosomal from mononucleosomal fragments
#            valley2: Second valley position separating mononucleosomal from dinucleosomal fragments
# Output: List containing data frame with percentage distributions and file path to saved TSV file
df_for_barplot_with_perc <- function(bam_file_path, save_dir, valley1, valley2) {
  # load library
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(glue)
  })

  if (is.null(save_dir)) {
    stop("Error: Need output dir path: 'save_dir' is NULL.")
  }

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  bamWidths_list <- list()

  for (bamFile in bam_file_path) {
    temp <- readGAlignmentPairs(bamFile)
    locus <- data.frame(
      first_start = start(temp@first),
      first_end   = end(temp@first),
      last_start  = start(temp@last),
      last_end    = end(temp@last)
    )
    start <- rowMin(as.matrix(locus))
    end <- rowMax(as.matrix(locus))

    seqnames <- as.vector(seqnames(temp))
    bam <- makeGRangesFromDataFrame(
      data.frame(seqnames = seqnames, strand = "*", start = start, end = end)
    )

    crf_name <- tools::file_path_sans_ext(basename(bamFile))
    bamWidths_list[[crf_name]] <- width(bam)
  }

  result_df <- data.frame(
    crf_pair = character(),
    pct_low = numeric(),
    pct_mid = numeric(),
    pct_high = numeric()
  )

  for (i in seq_along(bamWidths_list)) {
    widths <- bamWidths_list[[i]]
    crf_name <- names(bamWidths_list)[i]

    total <- length(widths)
    pct_low <- sum(widths < valley1) / total
    pct_mid <- sum(widths >= valley1 & widths < valley2) / total
    pct_high <- sum(widths >= valley2) / total

    result_df <- rbind(result_df, data.frame(
      crf_pair = crf_name,
      pct_low = pct_low,
      pct_mid = pct_mid,
      pct_high = pct_high
    ))
  }

  colnames(result_df)[2] <- glue::glue("subnucleo")
  colnames(result_df)[3] <- glue::glue("monomer")
  colnames(result_df)[4] <- glue::glue("dimer")
  out_path <- file.path(save_dir, "df_for_barplot_with_perc_result.tsv")
  write.table(result_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  return(list(df_for_barplot_with_perc = result_df, file_path = out_path))
}
