differential_analysis <- function(cm_path, conditions, sample_names = NULL, out_dir = "./", col_cluster_file_path = NULL, ref_genome = "hg38", apply_annotation = TRUE, plot = TRUE) {
  # Step1: get differential regions
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Differential Regions")
  cat("\n", strrep("=", 40), "\n", sep = "")

  dr_result <- differential_regions(cm_path = cm_path, conditions = conditions, sample_names = sample_names, out_dir = out_dir, col_cluster_file_path = col_cluster_file_path)

  # Step2: generate differential heatmap
  if (plot) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Differential Heatmap")
    cat("\n", strrep("=", 40), "\n", sep = "")
    differential_heatmap(dr_result = dr_result, out_dir = out_dir)
  }

  if (apply_annotation) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Pathway Annotation")
    cat("\n", strrep("=", 40), "\n", sep = "")
    pathway_dir <- file.path(out_dir, "pathway_enrichment")
    differential_pathway_annotation(dr_result = dr_result, out_dir = pathway_dir, ref_genome = ref_genome, plot = plot)
  }
  invisible(dr_result)
}