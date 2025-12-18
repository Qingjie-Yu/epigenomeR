# TFBS Enrichment Analysis
# Description: This function analyzes the enrichment of transcription factor binding motifs in target genomic regions compared to control regions.
#
# Parameter: 
#   target_region: A GRanges object containing the target genomic regions for enrichment analysis.
#   control_region: A GRanges object containing the control genomic regions for background comparison.
#   regions: Size in bp to resize all regions around their center (default: 800)
#   ref_genome: Reference genome version - "hg38" or "mm10" (default: "hg38")
#   functional_region: A GRanges object or NULL. Optional functional regions (e.g., open chromatin, accessible regions) to restrict motif analysis. If provided, only motif sites overlapping these regions will be considered.
#   out_path: Output path to save enrichment results table (default: "./TFBS_enrichment.tsv")
#   style: A string specifying the desired seqlevelsStyle (e.g., "UCSC", "Ensembl"). If NULL, auto-detect from target_region; otherwise, enforce all objects to use the specified style. (default: NULL)
#
# Output: A data frame where each row corresponds to a motif with enrichment statistics (odds ratio, p-value, FDR) saved to output_path

TFBS_enrichment <- function(target_region, control_region, regions = NULL, out_path = "./TFBS_enrichment.tsv",  functional_region = NULL, ref_genome = "hg38", style = NULL) {
  # Load packages
  suppressPackageStartupMessages({
    library(IRanges)
    library(GenomicRanges)
    library(data.table)
    library(glue)
  })

  if (!is.null(regions)) {
    target_region <- resize(target_region, width = regions, fix='center')
    control_region <- resize(control_region, width = regions, fix='center')
  }

  # Load TFBS library
  TFBS_library_file <- switch(ref_genome,
    hg38 = "TFBS_lib_hg38.rds",
    mm10 = "TFBS_lib_mm10.rds",
    stop("Invalid ref_genome")
  )
  TFBS_library_file <- download_rds(TFBS_library_file)
  TFBS_library <- readRDS(TFBS_library_file)
  message(glue("Using reference genome {ref_genome} with {length(TFBS_library)} TFs"))

  if (is.null(style)) {
    style <- seqlevelsStyle(target_region)[1]
  } else {
    seqlevelsStyle(target_region) <- style
  }
  seqlevelsStyle(control_region) <- style
  seqlevelsStyle(TFBS_library) <- style

  # Filter by functional regions
  if(!is.null(functional_region)){
    seqlevelsStyle(functional_region) <- style
    message(glue("Filtering motif site using {length(functional_region)} functional regions"))
    TFBS_library <- endoapply(TFBS_library, subsetByOverlaps, functional_region)
    TFBS_library <- TFBS_library[lengths(TFBS_library) > 0]
    if (sum(lengths(TFBS_library)) == 0) {
      warning("No motif sites overlap with the provided functional regions.")
      return(NULL)
    }
  }

  # Count overlaps
  message("Counting overlaps and performing enrichment analysis...")
  target_overlap <- countOverlaps(TFBS_library, target_region)
  control_overlap <- countOverlaps(TFBS_library, control_region)
  n_target <- length(target_region)
  n_control <- length(control_region)
  tf_names <- names(TFBS_library)

  test_result <- data.table(
    TF = tf_names,
    target_hit = target_overlap,
    control_hit = control_overlap,
    target_off = n_target - target_overlap,
    control_off = n_control - control_overlap
  )

  test_result[, c("odds_ratio", "pvalue", "odds_ratio_se") := {
    # Add pseudocount to avoid zeros
    a <- target_hit + 1
    b <- control_hit + 1
    c <- target_off + 1
    d <- control_off + 1
    
    # Pre-allocate vectors
    odds_ratios <- numeric(.N)
    pvalues <- numeric(.N)
    odds_ses <- numeric(.N)
    
    for (i in seq_len(.N)) {
      test_table <- matrix(c(a[i], b[i], c[i], d[i]), nrow = 2,
                          dimnames = list(c("target", "control"), c("hit", "off")))
      
      odds_ses[i] <- sqrt(sum(1 / test_table))
      test_res <- fisher.test(test_table, alternative = "greater")
      odds_ratios[i] <- test_res$estimate
      pvalues[i] <- test_res$p.value
    }
    
    list(odds_ratios, pvalues, odds_ses)
  }]

  test_result[, FDR := p.adjust(pvalue, method = "BH")]
  setorder(test_result, pvalue)
  res <- as.data.frame(test_result)
  rownames(res) <- res$TF
  res <- res[, -1]
  write.table(res, file = out_path, sep = "\t", quote = FALSE, col.names = NA)

  n_significant <- sum(res$FDR < 0.05)
  message(glue("Results saved to {out_path}"))
  message(glue("Found {n_significant} significant motifs at FDR < 0.05"))
  invisible(out_path)
}
