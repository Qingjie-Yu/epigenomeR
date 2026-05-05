clustering_wrapper <- function(cm_paths, row_km, out_dir, sample_names, crf_names = NULL, transformations  = c("libnorm", "log2p1"), cluster_mode = "multi", apply_filter = TRUE, filter_regions = NULL, filter_mode = "union", seed = 42, feature_distance = NULL, apply_zscore = TRUE, apply_annotation = TRUE, ref_genome = "hg38", ref_source = "knownGene", distributions = c("genic", "ccre")) {
  # ── Transform ──
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Applying Transformations")
  cat("\n", strrep("=", 40), "\n", sep = "")

  transform_dir <- file.path(out_dir, "transformed")
  transformed_paths <- vapply(cm_paths, function(p) {
    apply_transformations(
      cm_path         = p,
      out_dir         = transform_dir,
      transformations = transformations
    )
  }, character(1))
  names(transformed_paths) <- sample_names

  # ── Clustering ───
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Clustering (mode: ", cluster_mode, ")")
  cat("\n", strrep("=", 40), "\n", sep = "")

  if (cluster_mode == "single") {
    row_paths <- clustering_single_crf(
      cm_paths              = transformed_paths,
      sample_names          = sample_names,
      row_km                = row_km,
      out_dir               = out_dir,
      crf_names             = crf_names,
      seed                  = seed,
      apply_filter          = apply_filter,
      filter_regions        = filter_regions,
      feature_distance      = feature_distance,
      apply_zscore          = apply_zscore
    )
  } else if (cluster_mode == "multi") {
    row_paths <- clustering_multi_crf(
      cm_paths              = transformed_paths,
      sample_names          = sample_names,
      row_km                = row_km,
      out_dir               = out_dir,
      crf_names             = crf_names,
      apply_filter          = apply_filter,
      filter_regions        = filter_regions,
      filter_mode           = filter_mode,
      seed                  = seed,
      order_clusters        = TRUE,
      cluster_linkage       = "complete",
      order_within_clusters = TRUE,
      feature_distance      = feature_distance,
      apply_zscore          = apply_zscore
    )
  } else {
    stop("cluster_mode must be 'single' or 'multi'")
  }

  # Annotation
  if (apply_annotation) {
    if (cluster_mode == "single") {
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
    } else {
      cat("\n", strrep("=", 40), "\n", sep = "")
      cat("  Annotation — Multi-CRF")
      cat("\n", strrep("=", 40), "\n", sep = "")

      genomic_dir  <- file.path(out_dir, "genomic_distribution")
      tfbs_dir     <- file.path(out_dir, "TFBS_enrichment")
      dir.create(genomic_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(tfbs_dir,    recursive = TRUE, showWarnings = FALSE)

      clustering_genomic_distribution(
        row_cluster_file_path = row_paths,
        out_dir               = genomic_dir,
        distributions         = distributions,
        ref_genome            = ref_genome,
        ref_source            = ref_source
      )

      clustering_TFBS_enrichment(
        row_cluster_file_path = row_paths,
        out_dir               = tfbs_dir,
        ref_genome            = ref_genome,
        ref_source            = ref_source
      )
    }

  }
  invisible(row_paths)
}