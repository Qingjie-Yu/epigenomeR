dr_result_to_grl <- function(dr_result, out_dir = "./") {
  gr_list <- list()

  for (grp_name in names(dr_result)) {
    for (cmp_tag in names(dr_result[[grp_name]])) {
      counts_path <- dr_result[[grp_name]][[cmp_tag]]
      if (is.null(counts_path) || !file.exists(counts_path)) next

      dt  <- data.table::fread(counts_path, select = "pos")
      pos <- as.character(dt$pos)

      chr   <- sub("_(\\d+)_(\\d+)$", "", pos)
      start <- as.integer(sub(".*_(\\d+)_\\d+$", "\\1", pos))
      end   <- as.integer(sub(".*_(\\d+)$", "\\1", pos))

      gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(start = start, end = end)
      )

      nm <- paste0(cmp_tag, "_", grp_name)
      gr_list[[nm]] <- gr
    }
  }

  GenomicRanges::GRangesList(gr_list)
}


differential_pathway_annotation <- function(dr_result, out_dir = "./", ref_genome = "hg38", gene_sets = "MSigDB:H", plot = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  grl <- dr_result_to_grl(dr_result)
  pathway_annotation(query = grl, out_dir = out_dir, ref_genome = ref_genome, gene_sets = gene_sets, plot = plot)
}