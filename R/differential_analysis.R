differential_analysis <- function(cm_path, conditions, sample_names = NULL, out_dir = "./", col_cluster_file_path = NULL, ref_genome = "hg38", min_support = 2, lfc_threshold = 0.5, p_threshold = 0.05, p_type = c("fdr", "nominal", "bonferroni"), norm_method = c("tmm", "libsize"), apply_annotation = TRUE, plot = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Step1: get differential regions
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Differential Regions")
  cat("\n", strrep("=", 40), "\n", sep = "")

  dr_result <- differential_regions(cm_path = cm_path, conditions = conditions, sample_names = sample_names, out_dir = out_dir, col_cluster_file_path = col_cluster_file_path, min_support = min_support, lfc_threshold = lfc_threshold, p_threshold = p_threshold, p_type = p_type, norm_method = norm_method)

  # Step2: generate differential heatmap
  if (plot) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Differential Heatmap")
    cat("\n", strrep("=", 40), "\n", sep = "")
    differential_heatmap(dr_result = dr_result, out_dir = out_dir)
  }

  # Step3: pathway annotation
  if (apply_annotation) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Pathway Annotation")
    cat("\n", strrep("=", 40), "\n", sep = "")
    differential_pathway_annotation(dr_result = dr_result, out_dir = out_dir, ref_genome = ref_genome, plot = plot)
  }
  invisible(dr_result)
}