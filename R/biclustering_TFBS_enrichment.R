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
#   out_dir:
#     Output directory for all results (default: "./").
#   ref_genome:
#     Reference genome version (default: "hg38").
#     Options: "hg38", "mm10".
#   ref_source:
#     Gene annotation source used to define gene models and TSS coordinates.
#     Supported options are:
#       - "knownGene": Uses the UCSC knownGene annotation obtained via
#         the Bioconductor package TxDb.Hsapiens.UCSC.hg38.knownGene.
#       - "GENCODE": Uses GENCODE gene annotations (GENCODE v49).
#   control_rep:
#     Multiplier for control region generation (default: 1).
#     E.g., control_rep = 2 generates 2x as many control regions as input regions.
#   regions:
#     Region size in base pairs for analysis (default: 800).
#     All regions resized to this width centered on original midpoint.
#   plot:
#     Whether to generate heatmaps (default: TRUE). If FALSE, only performs enrichment analysis.
#
# Output: 
#   Control regions BED file: all_controls.bed
#   TFBS enrichment TSV files (one per cluster): TFBS_enrichment_cluster_<label>.tsv
#   Heatmap PDF (all filtered TFBS): TFBS_heatmap_all.pdf (if plot=TRUE)
#   log2 odds ratio matrix: odds_ratio_log2.csv (if plot=TRUE)
#   FDR matrix: FDR.csv (if plot=TRUE)

biclustering_TFBS_enrichment <- function(row_cluster_file_path, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", control_rep = 1, regions = 800, plot = TRUE, plot_n_top = 20) {
    # Load packages
    suppressPackageStartupMessages({
        library(data.table)
        library(GenomicRanges)
        library(ComplexHeatmap)
        library(circlize)
        library(rtracklayer)
    })

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
    control_grl <- lapply(names(row_grl), function(label) {
        cat("\n", "Processing cluster:", label, "with", length(row_grl[[label]]), "regions\n")
        target_gr <- row_grl[[label]]
        control_gr <- get_matched_control(target_gr = target_gr, ref_genome = ref_genome, style = style, n_rep = control_rep, regions = regions)
        return(control_gr)
    })
    names(control_grl) <- names(row_grl)
    control_gr <- unlist(GRangesList(control_grl), use.names = FALSE)
    # Eliminating bias caused by overlap
    control_gr_reduced <- reduce(control_gr)
    control_gr <- resize(control_gr_reduced, width = regions, fix = "center")

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    export(control_gr, file.path(out_dir, "all_controls.bed"))
    cat("\n", "Control regions saved to:", file.path(out_dir, "all_controls.bed"), "\n")

    # TFBS enrichment for each cluster
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  TFBS Enrichment Analysis")
    cat("\n", strrep("=", 40), "\n", sep = "")
    tsv_paths <- lapply(names(row_grl), function(label) {
        cat("\n", "Processing TFBS enrichment for:", label, "\n")
        out_path <- file.path(out_dir, paste0("TFBS_enrichment_cluster_", label, ".tsv"))
        TFBS_enrichment(target_region = row_grl[[label]], control_region = control_gr,out_path = out_path, ref_genome = ref_genome, style = style)
        out_path
    })
    tsv_paths <- unlist(tsv_paths)

    if (plot) {
        cat("\n", strrep("=", 40), "\n", sep = "")
        cat("  TFBS Enrichment Heatmap Visualization")
        cat("\n", strrep("=", 40), "\n", sep = "")
        TFBS_enrichment_heatmap(tsv_path = tsv_paths, label = names(row_grl), out_dir = out_dir, top_n = plot_n_top)
    }
}

