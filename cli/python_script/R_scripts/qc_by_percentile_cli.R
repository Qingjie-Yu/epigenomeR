#!/usr/bin/env Rscript
# CLI for QC-by-percentile module

suppressPackageStartupMessages({
  library(argparse)
  library(multiEpiCore)
})


parser <- ArgumentParser(
  prog = "qc_by_percentile_cli.R",
  description = "Perform QC on BAM/BED files using percentile filtering and optionally generate heatmap"
)

# Required arguments
parser$add_argument(
  "--file_paths",
  required = TRUE,
  help = "Comma-separated BAM or BED file paths"
)

parser$add_argument(
  "--out_dir",
  required = TRUE,
  help = "Output directory"
)

# Optional arguments
parser$add_argument(
  "--filtered_percentile",
  type = "double",
  default = 0.25,
  help = "Percentile threshold for QC filtering (default: 0.25)"
)

parser$add_argument(
  "--no_plot",
  action = "store_true",
  default = FALSE,
  help = "Disable heatmap generation"
)

parser$add_argument(
  "--split_pair_by",
  default = "-",
  help = "Delimiter for splitting CRF pairs (default: '-')"
)

parser$add_argument(
  "--group_csv",
  default = NULL,
  help = "Optional CSV file with CRF grouping info"
)

parser$add_argument(
  "--crf_col",
  default = "crf",
  help = "Column name for CRF identifiers in group CSV (default: crf)"
)

parser$add_argument(
  "--category_col",
  default = "category",
  help = "Column name for CRF categories in group CSV (default: category)"
)

args <- parser$parse_args()


tryCatch({

  # Parse input file paths
  file_path <- trimws(strsplit(args$file_paths, ",")[[1]])

  if (length(file_path) == 0) {
    stop("No input files provided")
  }

  qc_by_percentile(
    file_path = file_path,
    out_dir = args$out_dir,
    filtered_percentile = args$filtered_percentile,
    plot = !args$no_plot,
    split_pair_by = args$split_pair_by,
    group_csv = args$group_csv,
    crf_col = args$crf_col,
    category_col = args$category_col
  )

}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})
