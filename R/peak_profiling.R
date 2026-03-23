peak_profiling <- function(bedgraph_path, out_dir = "./", ref_genome = "hg38", auc_top_pct = 0.1, qvalue_cutoff = 0.05, fc_cutoff = 2, apply_annotation = TRUE, plot = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Step1: peak calling
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Peak Calling")
  cat("\n", strrep("=", 40), "\n", sep = "")
  peak_result <- peak_calling(bedgraph_path = bedgraph_path, out_dir = out_dir, ref_genome = ref_genome, auc_top_pct = auc_top_pct, qvalue_cutoff = qvalue_cutoff, fc_cutoff = fc_cutoff)

  # Step2: pathway annotation
  if (apply_annotation) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Pathway Annotation")
    cat("\n", strrep("=", 40), "\n", sep = "")
    peak_path <- unlist(peak_result)
    pathway_dir <- file.path(out_dir, "pathway_annotation")
    peak_pathway_annotation(peak_path = peak_path, out_dir = pathway_dir, ref_genome = ref_genome, plot = plot)
  }
}