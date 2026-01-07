# Biclustering TFBS Enrichment Pipeline
#
# This function performs a complete TFBS (Transcription Factor Binding Site) enrichment
# analysis workflow for biclustered genomic regions, including matched control generation,
# enrichment testing, and heatmap visualization.
#
# Parameters:
#   row_cluster_file_path:
#     Path to a tab-delimited file with clustered regions.
#     Required columns: 'feature' (format: chr_start_end), 'label' (cluster assignment).
#     Example: chr1_1000_2000 \t ClusterA
#   out_dir:
#     Output directory for all results. Created if doesn't exist. Default: "./"
#   ref_genome:
#     Reference genome version. Default: "hg38"
#     Supported: "hg38", "mm10"
#   ref_source:
#     Gene annotation source for control region generation. Default: "knownGene"
#     Options:
#       - "knownGene": UCSC knownGene from TxDb packages
#       - "GENCODE": GENCODE annotations (v49 for hg38, vM23 for mm10)
#   control_rep:
#     Number of control regions per input region. Default: 1
#     E.g., control_rep = 2 generates 2 matched controls per target region
#   regions:
#     Width of regions in base pairs. Default: 800
#     All regions are resized to this width, centered on original midpoint
#   plot:
#     Whether to generate enrichment heatmaps. Default: TRUE
#   plot_n_top:
#     Number of top enriched TFBSs to display in heatmap. Default: 20
#     TFBSs ranked by minimum FDR across all clusters
#   seed:
#     Random seed for reproducible control region generation. Default: 42
#
# Output:
#   Control regions BED file: all_controls.bed
#   TFBS enrichment TSV files (one per cluster): TFBS_enrichment_cluster_<label>.tsv
#   Heatmap PDF (all filtered TFBS): TFBS_heatmap_all.pdf (if plot=TRUE)
#   log2 odds ratio matrix: odds_ratio_log2.csv (if plot=TRUE)
#   FDR matrix: FDR.csv (if plot=TRUE)

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
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
  pos_df <- do.call(rbind, (strsplit(row_cluster$feature, "_")))
  colnames(pos_df) <- c("seqnames", "start", "end")
  row_cluster <- cbind(pos_df, row_cluster)
  row_gr <- makeGRangesFromDataFrame(row_cluster, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  style <- seqlevelsStyle(row_gr)[1]
  row_grl <- split(row_gr, row_gr$label)

  # Generate matched control regions
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Generating matched control regions")
  cat("\n", strrep("=", 40), "\n", sep = "")
  control_gr <- get_matched_control(query = row_grl, ref_genome = ref_genome, ref_source = ref_source, style = style, n_rep = control_rep, regions = regions)
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
  tsv_paths <- TFBS_enrichment(query = row_grl, contro = control_gr, out_dir = out_dir, ref_genome = ref_genome, style = style)

  if (plot) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  TFBS Enrichment Heatmap Visualization")
    cat("\n", strrep("=", 40), "\n", sep = "")
    TFBS_enrichment_heatmap(tsv_path = tsv_paths, label = names(row_grl), out_dir = out_dir, top_n = plot_n_top)
  }
}
