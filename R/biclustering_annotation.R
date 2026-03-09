

extract_region_to_grl <- function(row_cluster_file_path, region_col = "region", cluster_col = "cluster") {
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c(region_col, cluster_col) %in% colnames(row_cluster))) {
    stop(glue("Input file must contain '{region_col}' and '{cluster_col}' columns"))
  }
  pos_df <- do.call(rbind, (strsplit(row_cluster[[region_col]], "_")))
  colnames(pos_df) <- c("seqnames", "start", "end")
  row_cluster <- cbind(pos_df, row_cluster)
  row_gr <- makeGRangesFromDataFrame(
    row_cluster,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  row_grl <- split(row_gr, mcols(row_gr)[[cluster_col]])
  message("Loaded ", length(row_gr), " regions across ", length(row_grl), " clusters")
  row_grl
}

# Biclustering Genomic Distribution Pipeline
#
# Annotates biclustered genomic regions with distribution across genomic features.
#
# Parameters:
#   row_cluster_file_path:
#     Path to a tab-delimited file with clustered regions.
#     Required columns: 'region' (format: chr_start_end), 'cluster' (cluster assignment).
#   out_dir:
#     Output directory for results. Default: "./"
#   distributions:
#     Annotation types to compute. Default: c("genic", "ccre")
#     Valid options: "genic", "ccre", "chromhmm", "repeat"
#   ref_genome:
#     Reference genome version. Default: "hg38"
#     Supported: "hg38", "mm10"
#   ref_source:
#     Gene annotation source. Default: "knownGene"
#     Options: "knownGene", "GENCODE"
#   mode:
#     Assignment mode for overlapping annotations. Default: "nearest"
#     Options: "nearest", "weighted"
#   plot:
#     Whether to generate distribution plots. Default: TRUE
biclustering_genomic_distribution <- function(row_cluster_file_path, out_dir = "./", distributions = c("genic", "ccre"), ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE) {
  # Validate parameters
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

  # Read row cluster file and convert to GRangesList
  row_grl <- extract_region_to_grl(row_cluster_file_path=row_cluster_file_path)

  # Genomic distribution for each cluster
  genomic_distribution(query = row_grl, out_dir = out_dir, distributions = distributions, ref_genome = ref_genome, ref_source = ref_source, mode = mode, plot = plot)
  message("Distribution annotation complete")
}

# Biclustering TFBS Enrichment Pipeline
#
# Performs TFBS enrichment analysis for biclustered genomic regions,
# including matched control generation, enrichment testing, and heatmap visualization.
#
# Parameters:
#   row_cluster_file_path:
#     Path to a tab-delimited file with clustered regions.
#     Required columns: 'region' (format: chr_start_end), 'cluster' (cluster assignment).
#   out_dir:
#     Output directory for all results. Created if doesn't exist. Default: "./"
#   ref_genome:
#     Reference genome version. Default: "hg38". Supported: "hg38", "mm10"
#   ref_source:
#     Gene annotation source for control region generation. Default: "knownGene"
#     Options: "knownGene", "GENCODE"
#   control_rep:
#     Number of matched control regions per input region. Default: 1
#   regions:
#     Width of regions in base pairs. Default: 800
#   plot:
#     Whether to generate enrichment heatmaps. Default: TRUE
#   plot_n_top:
#     Number of top enriched TFBSs to display in heatmap. Default: 20
#   seed:
#     Random seed for reproducible control region generation. Default: 42
#
# Output:
#   all_controls.bed: control regions BED file
#   TFBS_enrichment_cluster_<label>.tsv: enrichment results per cluster
#   TFBS_heatmap_all.pdf, odds_ratio_log2.csv, FDR.csv: heatmap and matrices (if plot=TRUE)
biclustering_TFBS_enrichment <- function(row_cluster_file_path, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", control_rep = 1, regions = 800, plot = TRUE, plot_n_top = 20, seed = 42) {
  # Load packages
  suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ComplexHeatmap)
    library(circlize)
    library(rtracklayer)
  })

  # Validate inputs
  if (!file.exists(row_cluster_file_path)) {
    stop("Input file does not exist: ", row_cluster_file_path)
  }
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

  # Read row cluster file and convert to GRanges
  row_grl <- extract_region_to_grl(row_cluster_file_path=row_cluster_file_path)

  # Generate matched control regions
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Generating matched control regions")
  cat("\n", strrep("=", 40), "\n", sep = "")
  control_gr <- get_matched_control(query = row_grl, ref_genome = ref_genome, ref_source = ref_source, n_rep = control_rep, regions = regions)
  # Eliminating bias caused by overlap
  control_gr_reduced <- reduce(control_gr)
  control_gr <- resize(control_gr_reduced, width = regions, fix = "center")
  export(control_gr, file.path(out_dir, "all_controls.bed"), format = "bed")
  cat("Generated", length(control_gr), "unique control regions\n")
  cat("\n", "Control regions saved to:", file.path(out_dir, "all_controls.bed"), "\n")

  # TFBS enrichment for each cluster
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  TFBS Enrichment Analysis")
  cat("\n", strrep("=", 40), "\n", sep = "")
  tsv_paths <- TFBS_enrichment(query = row_grl, contro = control_gr, out_dir = out_dir, ref_genome = ref_genome)

  if (plot) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  TFBS Enrichment Heatmap Visualization")
    cat("\n", strrep("=", 40), "\n", sep = "")
    TFBS_enrichment_heatmap(tsv_path = tsv_paths, label = names(row_grl), out_dir = out_dir, top_n = plot_n_top)
  }
}

# Biclustering Pathway Annotation Pipeline
#
# Performs pathway enrichment analysis on biclustered genomic regions using rGREAT.
#
# Parameters:
#   row_cluster_file_path:
#     Path to a tab-delimited file with clustered regions.
#     Required columns: 'region' (format: chr_start_end), 'cluster' (cluster assignment).
#   out_dir:
#     Output directory for results. Default: "./"
#   ref_genome:
#     Reference genome version. Default: "hg38"
#     Supported: "hg38", "mm10"
#   gene_sets:
#     Gene set collection for enrichment testing. Default: "MSigDB:H"
#   plot:
#     Whether to generate a bubble plot. Default: TRUE
biclustering_pathway_annotation <- function(row_cluster_file_path, out_dir = "./", ref_genome = "hg38", gene_sets = "MSigDB:H", plot = TRUE) {
  # Validate parameters
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }

  # Read row cluster file and convert to GRangesList
  row_grl <- extract_region_to_grl(row_cluster_file_path = row_cluster_file_path)

  # Pathway annotation for each cluster
  pathway_annotation(query = row_grl, out_dir = out_dir, ref_genome = ref_genome, gene_sets = gene_sets, plot = plot)
  message("Pathway annotation complete")
}