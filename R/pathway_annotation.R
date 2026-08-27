pathway_annotation_plot_bubble <- function(table_list, out_dir, min_padj_for_color = 1e-300, color_cap_pct = 0.99) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scales)
  })

  message(sprintf("Plotting bubble chart: %d targets, pathways from table_list", length(table_list)))

  df_long <- dplyr::bind_rows(lapply(names(table_list), function(nm) {
    tb <- table_list[[nm]]
    tb$target <- nm
    tb
  }))

  message(sprintf("  df_long: %d rows, %d unique pathways", nrow(df_long), length(unique(df_long$pathway))))

  df_long <- df_long |>
    dplyr::mutate(
      padj_safe     = dplyr::if_else(is.na(padj), NA_real_, pmax(padj, min_padj_for_color)),
      neglog10_padj = -log10(padj_safe),
      size_val      = log2(1 + pmax(fold, 0))
    )

  # determine color scale upper bound
  if (is.null(color_cap_pct)) {
    cap_neglog10 <- max(df_long$neglog10_padj, na.rm = TRUE)
    message(sprintf("  color_cap_pct not set -> using max(neglog10_padj) = %.2f as color cap", cap_neglog10))
  } else {
    stopifnot(is.numeric(color_cap_pct), color_cap_pct > 0, color_cap_pct <= 1)
    cap_neglog10 <- stats::quantile(df_long$neglog10_padj, probs = color_cap_pct, na.rm = TRUE, names = FALSE)
    message(sprintf("  color_cap_pct = %.3g -> cap_neglog10 set to %.2f", color_cap_pct, cap_neglog10))
  }
  if (!is.finite(cap_neglog10) || cap_neglog10 <= 0) {
    cap_neglog10 <- 1
    message("  cap_neglog10 was non-finite/<=0, falling back to 1")
  }

  df_long <- df_long |>
    dplyr::mutate(neglog10_cap = pmin(neglog10_padj, cap_neglog10))

  pathway_levels <- df_long |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(best_padj = suppressWarnings(min(padj, na.rm = TRUE)), .groups = "drop") |>
    dplyr::arrange(best_padj) |>
    dplyr::pull(pathway)

  df_long <- df_long |>
    dplyr::mutate(
      target  = factor(target,  levels = names(table_list)),
      pathway = factor(pathway, levels = rev(pathway_levels))
    )

  n_targets  <- length(table_list)
  n_pathways <- length(unique(df_long$pathway))

  max_pathway_nchar <- max(nchar(as.character(pathway_levels)))
  yaxis_margin      <- max_pathway_nchar * 0.07

  max_target_nchar <- max(nchar(names(table_list)))
  xaxis_margin     <- max_target_nchar * 0.07

  col_width  <- 0.5
  row_height <- 0.25
  legend_width <- 2.5

  width  <- yaxis_margin + n_targets * col_width + legend_width
  height <- xaxis_margin + n_pathways * row_height + 0.5

  width  <- max(width,  5)
  height <- max(height, 4)

  message(sprintf("  PDF size: %.1f x %.1f inch (%d targets × %d pathways, y-margin: %.1f, x-margin: %.1f)",
                  width, height, n_targets, n_pathways, yaxis_margin, xaxis_margin))

  p <- ggplot(df_long, aes(x = target, y = pathway)) +
    geom_point(aes(size = size_val, color = neglog10_cap)) +
    scale_color_gradient(
      name   = expression(-log[10](padj)),
      low    = "grey90",
      high   = "blue4",
      limits = c(0, cap_neglog10),
      oob    = scales::squish
    ) +
    scale_size_continuous(name = "log2(1 + fold_enrichment)") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position    = "right",
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(linewidth = 0.3),
      panel.grid.minor   = element_blank()
    )

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  plot_path <- file.path(out_dir, "pathway_annotation.pdf")
  ggsave(plot_path, plot = p, width = width, height = height, limitsize = FALSE)
  message(sprintf("  Saved: %s", plot_path))
}


pathway_annotation <- function(query, out_dir = "./", ref_genome = "hg38", msigdb_collection = "H", plot = TRUE, color_cap_pct = 0.99) {
  # Load required packages
  suppressPackageStartupMessages({
    library(rGREAT)
    library(rtracklayer)
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(msigdbr)
  })

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Set up parallel processing
  BPPARAM <- get_BPPARAM()

  # Parameter validation
  genome_map <- list(
    hg38 = list(species = "Homo sapiens", db_species = "HS",
                tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene"),
    mm10 = list(species = "Mus musculus", db_species = "MM",
                tss_source = "TxDb.Mmusculus.UCSC.mm10.knownGene")
  )
  if (!ref_genome %in% names(genome_map)) {
    stop("Unsupported genome: ", ref_genome)
  }
  species    <- genome_map[[ref_genome]]$species
  db_species <- genome_map[[ref_genome]]$db_species
  tss_source <- genome_map[[ref_genome]]$tss_source

  # collection mapping
  msigdb_collection_map <- data.frame(
    HS = c("H",  "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"),
    MM = c("MH", "M1", "M2", "M3", NA,   "M5", NA,   NA,   "M8"),
    stringsAsFactors = FALSE
  )

  resolve_msigdb_collection <- function(db_species, requested_collection) {
    valid <- msigdbr_collections(db_species = db_species)$gs_collection
    if (requested_collection %in% valid) {
      return(requested_collection)
    }

    other_species <- setdiff(c("HS", "MM"), db_species)
    idx <- match(requested_collection, msigdb_collection_map[[other_species]])
    if (is.na(idx)) {
      stop(sprintf(
        "Collection '%s' is not available for db_species '%s', and no mapping exists in msigdb_collection_map. Available collections: %s",
        requested_collection, db_species, paste(valid, collapse = ", ")
      ))
    }

    mapped <- msigdb_collection_map[[db_species]][idx]
    if (is.na(mapped)) {
      stop(sprintf(
        "Collection '%s' has no %s-native equivalent (checked msigdb_collection_map). Available collections: %s",
        requested_collection, db_species, paste(unique(valid), collapse = ", ")
      ))
    }

    warning(sprintf(
      "Collection '%s' not found natively for db_species '%s'; substituting mapped equivalent '%s' instead.",
      requested_collection, db_species, mapped
    ))
    mapped
  }

  collection <- resolve_msigdb_collection(db_species, msigdb_collection)
  message(sprintf("  Using collection '%s' (requested: '%s', db_species: '%s')",
                  collection, msigdb_collection, db_species))

  # Run rGREAT
  gene_sets_data <- msigdbr(species = species, db_species = db_species, collection = collection) |>
  (\(df) split(as.character(df$ncbi_gene), df$gs_name))()
  message(sprintf("  Loaded %d gene sets from MSigDB", length(gene_sets_data)))

  gr_names <- names(query)
  table_list <- BiocParallel::bplapply(gr_names, function(nm) {
    gr <- query[[nm]]
    res <- great(gr, gene_sets = gene_sets_data, tss_source = tss_source)
    tb <- res@table |>
      dplyr::transmute(
        pathway = if (collection %in% c("H", "MH")) str_replace(id, "^HALLMARK_", "") else id,
        hits_region = observed_region_hits,
        fold = fold_enrichment,
        p = p_value,
        padj = p_adjust,
        hits_gene = observed_gene_hits
      )
    tb_path <- file.path(out_dir, paste0("pathway_annotation_", nm, ".tsv"))
    write.table(tb, tb_path, sep = "\t", quote = FALSE, row.names = FALSE)
    tb
  }, BPPARAM = BPPARAM)
  names(table_list) <- gr_names

  # Dot/Bubble plot
  if (plot) {
    pathway_annotation_plot_bubble(table_list=table_list, out_dir=out_dir, color_cap_pct=color_cap_pct)
  }
}