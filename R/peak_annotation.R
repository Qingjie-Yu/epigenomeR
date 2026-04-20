read_peak_to_grl <- function(peak_path, pattern = "_peaks\\.narrowPeak$") {
  targets <- basename(peak_path) |> stringr::str_replace(pattern, "")
  grl <- GRangesList(lapply(peak_path, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t", comment.char = "#")
    gr <- GRanges(
      seqnames = df[[1]],
      ranges = IRanges(start = df[[2]], end = df[[3]])
    )
    gr
  }))
  names(grl) <- targets
  grl
}

# Peak Genomic Distribution Pipeline
#
# Annotates a set of peak BED files with distribution across genomic features.
#
# Parameters:
#   peak_path:
#     Path(s) to input PEAK files
#   out_dir:
#     Output directory for results. Default: "./"
#   pattern:
#     Regex pattern to match PEAK filenames (used for sample name extraction). Default: "_peaks\\.narrowPeak$"
#   distributions:
#     Annotation types to compute. Default: c("genic", "ccre")
#     Valid options: "genic", "ccre", "chromhmm", "repeat"
#   ref_genome:
#     Reference genome version. Default: "hg38". Supported: "hg38", "mm10"
#   ref_source:
#     Gene annotation source. Default: "knownGene"
#     Options: "knownGene", "GENCODE"
#   mode:
#     Assignment mode for overlapping annotations. Default: "nearest"
#     Options: "nearest", "weighted"
#   plot:
#     Whether to generate distribution plots. Default: TRUE
peak_genomic_distribution <- function(peak_path, out_dir = "./", pattern = "_peaks\\.narrowPeak$", distributions = c("genic", "ccre"), ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE) {
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
  })

  # Parameter validation
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  if (!ref_source %in% c("knownGene", "GENCODE")) {
    stop("ref_source must be 'knownGene' or 'GENCODE'")
  }
  if (!mode %in% c("nearest", "weighted")) {
    stop("mode must be 'nearest' or 'weighted'")
  }

  valid_distributions <- c("genic", "ccre", "chromhmm", "repeat")
  invalid_annos <- setdiff(distributions, valid_distributions)
  if (length(invalid_annos) > 0) {
    stop(
      "Invalid annotation types: ", paste(invalid_annos, collapse = ", "),
      "\nValid options are: ", paste(valid_distributions, collapse = ", ")
    )
  }
  if (length(distributions) == 0) {
    stop(
      "distributions parameter must contain at least one annotation type: ",
      paste(valid_distributions, collapse = ", ")
    )
  }

  # PEAK -> GRangesList
  grl <- read_peak_to_grl(peak_path=peak_path, pattern=pattern)

  # Genomic distribution for each target
  genomic_distribution(query = grl, out_dir = out_dir, distributions = distributions, ref_genome = ref_genome, ref_source = ref_source, mode = mode, plot = plot)
  message("Distribution annotation complete")
}

# Peak Pathway Annotation Pipeline
#
# Annotates a set of peak PEAK files with pathway and gene set enrichment.
#
# Parameters:
#   peak_path:
#     Path(s) to input PEAK files
#   out_dir
#     Output directory for results. Default: "./"
#   ref_genome:
#     Reference genome version. Default: "hg38". Supported: "hg38", "mm10"
#   pattern:
#     Regex pattern to match PEAK filenames (used for sample name extraction). Default: "_peaks\\.narrowPeak$"
#   plot:
#     Whether to generate annotation plots. Default: TRUE
peak_pathway_annotation <- function(peak_path, out_dir = "./", ref_genome = "hg38", msigdb_collection = "H", pattern = "_peaks\\.narrowPeak$", plot = TRUE) {
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
    library(ggplot2)
  })

  # Parameter validation
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  # PEAK -> GRangesList
  grl <- read_peak_to_grl(peak_path = peak_path, pattern = pattern)

  # Pathway annotation for each target
  pathway_annotation(query = grl, out_dir = out_dir, ref_genome = ref_genome, msigdb_collection = msigdb_collection, plot = plot)
  message("Pathway annotation complete")
}