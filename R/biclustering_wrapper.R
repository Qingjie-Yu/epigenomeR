# Biclustering Analysis Wrapper
# 
# Performs biclustering analysis on chromatin accessibility count matrices,
# with optional filtering, k-means clustering, and gene annotation of genomic regions.
# 
# Parameters:
#   cm_path: Path to the count matrix `.feather` file or a vector of paths to 
#            multiple count matrix files to be merged.
#   out_dir: Directory to save all output files (cluster tables, filtered matrices, 
#            annotation plots).
#   apply_filter: Logical. Controls whether to further filter genomic regions.
#       - When the genome was segmented into equal-sized bins
#         (i.e., the count matrix was built using a numeric `regions` argument),
#         this should be TRUE so that low-information or uninformative bins can be removed.
#       - When the user supplied specific genomic intervals of interest
#         (i.e., the count matrix was built using a region file path),
#         this should be FALSE because no additional filtering is needed.
#   row_km: Number of k-means clusters for rows (genomic regions). Default: 15
#   col_km: Number of k-means clusters for columns (CRF pairs). Default: 3
#   apply_annotation: Logical. Controls whether to annotate genomic regions to nearby genes.
#       - When the genome was segmented into bins (numeric `regions`), this should be TRUE,
#         since bins lack inherent biological meaning and benefit from gene-level annotation.
#       - When the user provided specific regions of interest (region file path),
#         this should be FALSE, as those regions are already meaningful and do not require annotation.
#   ref_genome: Reference genome version for annotation and control region generation.
#               Default: "hg38". Options: "hg38", "mm10"
#   ref_source: Gene annotation source for control region generation and gene annotation. Default: "knownGene"
#        - "knownGene": UCSC knownGene from TxDb packages
#        - "GENCODE": GENCODE annotations (v49 for hg38, vM23 for mm10)
#   distributions: Character vector specifying genomic feature distributions to analyze for TFBS enrichment. Default: c("genic", "ccre")
#        - Options include: "genic", "ccre", "cpg", "repeat"
#   plot: Logical. Whether to generate diagnostic plots during filtering and biclustering steps. Default: TRUE



biclustering_wrapper <- function(cm_path, out_dir, apply_filter = TRUE, row_km = 15, col_km = 3, apply_annotation = TRUE, ref_genome = "hg38", ref_source = "knownGene", distributions = c("genic","ccre"), plot = TRUE) {
    # Step1: Merge all the count matrix files
    if (is.vector(cm_path) && length(cm_path) > 1) {
        cat("\n", strrep("=", 40), "\n", sep = "")
        cat("  Merge all count matrix files")
        cat("\n", strrep("=", 40), "\n", sep = "")
        merged_cm_path <- merge_count_matrices(cm_path = cm_path, out_dir = out_dir)
    } else {
        merged_cm_path <- cm_path
    }

    # Step2: Apply transformation
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Apply transformation")
    cat("\n", strrep("=", 40), "\n", sep = "")
    transformed_cm_path <- apply_transformations(cm_path = merged_cm_path, out_dir = out_dir)

    # Step3: Filter highly variable regions
    if (apply_filter) {
        cat("\n", strrep("=", 40), "\n", sep = "")
        cat("  Filter highly variable regions")
        cat("\n", strrep("=", 40), "\n", sep = "")
        f_cm_path <- detect_hvr(transformed_cm_path = transformed_cm_path, out_dir = out_dir, plot = plot)
    } else {
        f_cm_path <- transformed_cm_path
    }

    # Step3: Biclustering 
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Biclustering")
    cat("\n", strrep("=", 40), "\n", sep = "")
    cluster_list <- biclustering(cm_path = f_cm_path, row_km = row_km, col_km = col_km, out_dir = out_dir, plot =  plot)

    # Step4: Biclustering annotation
    if (apply_annotation) {
        cat("\n", strrep("=", 40), "\n", sep = "")
        cat("  Annotation")
        cat("\n", strrep("=", 40), "\n", sep = "")
        genomic_dir <- file.path(out_dir, "genomic_distribution")
        tfbs_dir <- file.path(out_dir, "TFBS_enrichment")
        dir.create(genomic_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(tfbs_dir, recursive = TRUE, showWarnings = FALSE)
        biclustering_genomic_distribution(row_cluster_file_path = cluster_list$row_table, out_dir = genomic_dir, distributions = distributions, ref_genome = ref_genome, ref_source = ref_source)
        biclustering_TFBS_enrichment(row_cluster_file_path = cluster_list$row_table, out_dir = tfbs_dir, ref_genome = ref_genome, ref_source = ref_source)
    }
}