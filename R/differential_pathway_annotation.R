dr_result_to_grl <- function(dr_result_cmp) {
  gr_list <- list()

  for (grp_name in names(dr_result_cmp)) {
    counts_path <- dr_result_cmp[[grp_name]]
    if (is.null(counts_path) || !file.exists(counts_path)) next

    dt  <- data.table::fread(counts_path, select = "pos")
    pos <- as.character(dt$pos)

    chr   <- sub("_(\\d+)_(\\d+)$", "", pos)
    start <- as.integer(sub(".*_(\\d+)_\\d+$", "\\1", pos))
    end   <- as.integer(sub(".*_(\\d+)$", "\\1", pos))

    gr_list[[grp_name]] <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges   = IRanges::IRanges(start = start, end = end)
    )
  }

  GenomicRanges::GRangesList(gr_list)
}

differential_pathway_annotation <- function(dr_result, out_dir = "./", ref_genome = "hg38", msigdb_collection = "H", plot = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (cmp_tag in names(dr_result)) {
    grl <- dr_result_to_grl(dr_result[[cmp_tag]])
    cmp_dir <- file.path(out_dir, cmp_tag, "pathway_annotation")
    if (!length(grl)) next
    pathway_annotation(query = grl, out_dir = cmp_dir, ref_genome = ref_genome, msigdb_collection = msigdb_collection, plot = plot)
  }
}