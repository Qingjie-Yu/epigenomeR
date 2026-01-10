#!/usr/bin/env Rscript
# CLI for biclustering wrapper module

suppressPackageStartupMessages({
  library(argparse)
  library(multiEpiCore)
})

parser <- ArgumentParser(
  prog = "biclustering_wrapper_cli.R",
  description = "Performs biclustering analysis on chromatin accessibility count matrices,
  include optional filtering, k-means clustering, and gene annotation of genomic regions"
)

# Required arguments
parser$add_argument("--cm_paths", required = TRUE,
                    help = "Comma-separated .feather count-matrix path(s)")
parser$add_argument("--out_dir", required = TRUE,
                    help = "Output directory")

# Optional arguments
parser$add_argument("--no_filter", action = "store_true", default = FALSE,
                    help = "Disable HVR filtering (sets apply_filter=FALSE)")
parser$add_argument("--no_annotation", action = "store_true", default = FALSE,
                    help = "Disable annotation/TFBS steps (sets apply_annotation=FALSE)")
parser$add_argument("--no_plot", action = "store_true", default = FALSE,
                    help = "Disable plotting (sets plot=FALSE)")

parser$add_argument("--ref_genome", default = "hg38",
                    help = "Reference genome: hg38 or mm10 (default: hg38)")
parser$add_argument("--ref_source", default = "knownGene",
                    help = "Annotation source: knownGene or GENCODE (default: knownGene)")
parser$add_argument("--distributions", default = "genic,ccre",
                    help = "Comma-separated distributions (default: genic,ccre)")

parser$add_argument("--row_km", type = "integer", default = 15,
                    help = "Number of k-means clusters for rows (default: 15)")
parser$add_argument("--col_km", type = "integer", default = 3,
                    help = "Number of k-means clusters for cols (default: 3)")

args <- parser$parse_args()

tryCatch({

  # Parse cm paths
  cm_paths <- trimws(strsplit(args$cm_paths, ",")[[1]])
  if (length(cm_paths) == 1L) cm_paths <- cm_paths[[1]]

  # Parse distributions
  distributions <- trimws(strsplit(args$distributions, ",")[[1]])

  # Map negative flags to function booleans
  apply_filter <- !isTRUE(args$no_filter)
  apply_annotation <- !isTRUE(args$no_annotation)
  plot <- !isTRUE(args$no_plot)

  biclustering_wrapper(
    cm_path = cm_paths,
    out_dir = args$out_dir,
    apply_filter = apply_filter,
    row_km = args$row_km,
    col_km = args$col_km,
    apply_annotation = apply_annotation,
    ref_genome = args$ref_genome,
    ref_source = args$ref_source,
    distributions = distributions,
    plot = plot
  )

}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})