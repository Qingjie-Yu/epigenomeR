read_peak_to_grl <- function(peak_path, pattern = "_peaks\\.narrowPeak$") {
  targets <- basename(peak_path) |> stringr::str_replace(pattern, "")

  grl_list <- lapply(seq_along(peak_path), function(i) {
    f <- peak_path[i]

    if (is.null(f) || !file.exists(f) || file.size(f) == 0L) {
      message("Skipping ", targets[i], ": no peak file or empty file.")
      return(NULL)
    }

    df <- read.table(f, header = FALSE, sep = "\t", comment.char = "#")
    if (nrow(df) == 0L) {
      message("Skipping ", targets[i], ": 0 peaks.")
      return(NULL)
    }

    GRanges(
      seqnames = df[[1]],
      ranges = IRanges(start = df[[2]], end = df[[3]])
    )
  })

  keep <- !vapply(grl_list, is.null, logical(1))
  grl <- GRangesList(grl_list[keep])
  names(grl) <- targets[keep]
  grl
}

read_tsv_to_grl <- function(tsv_path, pattern = "_sig\\.tsv$") {
  targets <- basename(tsv_path) |> stringr::str_replace(pattern, "")

  grl_list <- lapply(seq_along(tsv_path), function(i) {
    f <- tsv_path[i]

    if (is.null(f) || !file.exists(f) || file.size(f) == 0L) {
      message("Skipping ", targets[i], ": no tsv file or empty file.")
      return(NULL)
    }

    df <- read.table(f, header = TRUE, sep = "\t", comment.char = "#")
    if (nrow(df) == 0L) {
      message("Skipping ", targets[i], ": 0 regions.")
      return(NULL)
    }

    region <- df[[1]]
    parts <- stringr::str_split_fixed(region, "_", 3)
    GRanges(
      seqnames = parts[, 1],
      ranges = IRanges(start = as.integer(parts[, 2]), end = as.integer(parts[, 3]))
    )
  })

  keep <- !vapply(grl_list, is.null, logical(1))
  grl <- GRangesList(grl_list[keep])
  names(grl) <- targets[keep]
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
  if (length(grl) == 0L) {
    warning("No pairs with peaks remaining after reading peak files; skipping pathway annotation.")
    return(invisible(NULL))
  }

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
  if (length(grl) == 0L) {
    warning("No pairs with peaks remaining after reading peak files; skipping pathway annotation.")
    return(invisible(NULL))
  }

  # Pathway annotation for each target
  pathway_annotation(query = grl, out_dir = out_dir, ref_genome = ref_genome, msigdb_collection = msigdb_collection, plot = plot)
  message("Pathway annotation complete")
}

# Read narrowPeak files into a GRangesList, anchored on each peak's summit
# (column 10, "peak" = offset from chromStart) rather than the peak's
# geometric center. Each region is expanded to `width` bp around the summit.
read_summit_to_grl <- function(peak_path, pattern = "_peaks\\.narrowPeak$", width = 800) {
  targets <- basename(peak_path) |> stringr::str_replace(pattern, "")

  grl_list <- lapply(seq_along(peak_path), function(i) {
    f <- peak_path[i]

    if (is.null(f) || !file.exists(f) || file.size(f) == 0L) {
      message("Skipping ", targets[i], ": no peak file or empty file.")
      return(NULL)
    }

    df <- read.table(f, header = FALSE, sep = "\t", comment.char = "#")
    if (nrow(df) == 0L) {
      message("Skipping ", targets[i], ": 0 peaks.")
      return(NULL)
    }

    # narrowPeak columns: chrom, chromStart, chromEnd, name, score, strand,
    # signalValue, pValue, qValue, peak (summit offset from chromStart)
    chrom       <- df[[1]]
    chromStart  <- df[[2]]
    summit_pos  <- chromStart + df[[10]]   # absolute summit coordinate (0-based)

    half_width <- floor(width / 2)
    GRanges(
      seqnames = chrom,
      ranges   = IRanges(start = summit_pos - half_width + 1L,
                          end   = summit_pos + half_width)
    )
  })

  keep <- !vapply(grl_list, is.null, logical(1))
  grl <- GRangesList(grl_list[keep])
  names(grl) <- targets[keep]
  grl
}

# TFBS enrichment pipeline
peak_TFBS_enrichment <- function(peak_path, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", pattern = "_peaks\\.narrowPeak$", control_rep = 1, regions = 800, plot = TRUE, plot_n_top = 20, seed = 42) {
  # Load packages
  suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ComplexHeatmap)
    library(circlize)
    library(rtracklayer)
  })

  # Validate inputs
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }
  if (!ref_source %in% c("knownGene", "GENCODE")) {
    stop("Unsupported ref_source. Please use 'knownGene' or 'GENCODE'.")
  }
  if (control_rep < 1) {
    stop("control_rep must be at least 1")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # PEAK -> GRangesList, anchored on the summit (not the geometric center)
  grl <- read_summit_to_grl(peak_path = peak_path, pattern = pattern, width = regions)
  if (length(grl) == 0L) {
    warning("No pairs with peaks remaining after reading peak files; skipping TFBS enrichment.")
    return(invisible(NULL))
  }

  # Generate matched control regions
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Generating matched control regions")
  cat("\n", strrep("=", 40), "\n", sep = "")
  control_gr <- get_matched_control(query = grl, ref_genome = ref_genome, ref_source = ref_source, n_rep = control_rep, regions = regions)
  # Eliminating bias caused by overlap
  control_gr_reduced <- reduce(control_gr)
  control_gr <- resize(control_gr_reduced, width = regions, fix = "center")
  export(control_gr, file.path(out_dir, "all_controls.bed"), format = "bed")
  cat("Generated", length(control_gr), "unique control regions\n")
  cat("\n", "Control regions saved to:", file.path(out_dir, "all_controls.bed"), "\n")

  # TFBS enrichment for each pair
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  TFBS Enrichment Analysis")
  cat("\n", strrep("=", 40), "\n", sep = "")
  tsv_paths <- TFBS_enrichment(query = grl, control = control_gr, out_dir = out_dir, ref_genome = ref_genome)

  if (plot) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  TFBS Enrichment Heatmap Visualization")
    cat("\n", strrep("=", 40), "\n", sep = "")
    TFBS_enrichment_heatmap(tsv_path = tsv_paths, label = names(grl), out_dir = out_dir, top_n = plot_n_top)
  }
}