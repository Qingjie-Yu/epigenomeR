pathway_annotation_plot_bubble <- function(table_list, out_dir, min_padj_for_color = 1e-300, cap_neglog10 = 50) {
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
      neglog10_cap  = pmin(neglog10_padj, cap_neglog10),
      size_val      = log2(1 + pmax(fold, 0))
    )

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

  # longest pathway label → estimate y-axis margin
  max_pathway_nchar <- max(nchar(as.character(pathway_levels)))
  yaxis_margin      <- max_pathway_nchar * 0.07  # ~0.07 inch per character

  # longest target label → estimate x-axis margin (rotated 45°)
  max_target_nchar <- max(nchar(names(table_list)))
  xaxis_margin     <- max_target_nchar * 0.07

  # plot body: per-column and per-row cell size
  col_width  <- 0.5   # inch per target column
  row_height <- 0.25  # inch per pathway row
  legend_width <- 2.5 # inch for right-side legend

  width  <- yaxis_margin + n_targets * col_width + legend_width
  height <- xaxis_margin + n_pathways * row_height + 0.5  # 0.5 for top margin

  # floor values to avoid tiny PDFs
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

pathway_annotation <- function(query, out_dir = "./", ref_genome = "hg38", msigdb_collection = "H", plot = TRUE) {
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
  species <- switch(ref_genome,
    hg38 = "Homo sapiens",
    mm10 = "Mus musculus",
    stop("Unsupported genome: ", ref_genome)
  )
  tss_source <- switch(ref_genome,
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )

  # Run rGREAT
  gene_sets_data <- msigdbr(species = species, collection = msigdb_collection) |>
  (\(df) split(as.character(df$ncbi_gene), df$gs_name))()
  message(sprintf("  Loaded %d gene sets from MSigDB", length(gene_sets_data)))

  gr_names <- names(query)
  table_list <- BiocParallel::bplapply(gr_names, function(nm) {
    gr <- query[[nm]]
    res <- great(gr, gene_sets = gene_sets_data, tss_source = tss_source)
    tb <- res@table |>
      dplyr::transmute(
        pathway = str_replace(id, "^HALLMARK_", ""),
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
    pathway_annotation_plot_bubble(table_list=table_list, out_dir=out_dir)
  }
}