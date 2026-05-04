clustering_wrapper <- function(cm_paths, row_km, out_dir, sample_names, crf_names = NULL, apply_filter = TRUE, seed = 42, order_clusters = TRUE, cluster_linkage = "complete", order_within_clusters = TRUE, feature_distance = NULL, feature_linkage = NULL, plot = TRUE, apply_zscore = TRUE, show_column_names = TRUE, lower_range = NULL, upper_range = NULL, apply_annotation = TRUE, ref_genome = "hg38", ref_source = "UCSC") {
  # Clustering
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Clustering")
  cat("\n", strrep("=", 40), "\n", sep = "")

  row_paths <- clustering(
    cm_paths              = cm_paths,
    sample_names          = sample_names,
    row_km                = row_km,
    out_dir               = out_dir,
    crf_names             = crf_names,
    seed                  = seed,
    apply_filter          = apply_filter,
    order_clusters        = order_clusters,
    cluster_linkage       = cluster_linkage,
    order_within_clusters = order_within_clusters,
    feature_distance      = feature_distance,
    feature_linkage       = feature_linkage,
    plot                  = plot,
    apply_zscore          = apply_zscore,
    show_column_names     = show_column_names,
    lower_range           = lower_range,
    upper_range           = upper_range
  )

  # Annotation
  if (apply_annotation) {
    for (crf in names(row_paths)) {
      cat("\n", strrep("=", 40), "\n", sep = "")
      cat("  Annotation — CRF: ", crf)
      cat("\n", strrep("=", 40), "\n", sep = "")

      crf_out_dir  <- file.path(out_dir, crf)
      genomic_dir  <- file.path(crf_out_dir, "genomic_distribution")
      tfbs_dir     <- file.path(crf_out_dir, "TFBS_enrichment")
      dir.create(genomic_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(tfbs_dir,    recursive = TRUE, showWarnings = FALSE)

      clustering_genomic_distribution(
        row_cluster_file_path = row_paths[[crf]],
        out_dir               = genomic_dir,
        distributions         = distributions,
        ref_genome            = ref_genome,
        ref_source            = ref_source
      )

      clustering_TFBS_enrichment(
        row_cluster_file_path = row_paths[[crf]],
        out_dir               = tfbs_dir,
        ref_genome            = ref_genome,
        ref_source            = ref_source
      )
    }
  }
  invisible(row_paths)
}