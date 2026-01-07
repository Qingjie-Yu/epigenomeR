fix_region_colnames <- function(region_df) {
  cn <- colnames(region_df)
  cn_lower <- tolower(cn)

  pick_first <- function(cands) {
    idx <- match(cands, cn_lower)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) {
      return(NA)
    }
    return(idx[1])
  }

  idx_seq <- pick_first(c("seqnames", "seqname", "chr", "chrom", "chromosome"))
  idx_start <- pick_first(c("start", "start_pos", "chromstart", "begin", "pos", "position"))
  idx_end <- pick_first(c("end", "end_pos", "chromend", "stop"))
  idx_gene <- pick_first(c("gene_id", "geneid", "gene", "gene_name", "symbol", "gene_symbol", "ensembl_id", "ensembl_gene_id", "target_id", "id"))

  if (any(is.na(c(idx_seq, idx_start, idx_end)))) {
    stop("Error: cannot infer seqnames/start/end columns from region file.")
  }
  cn[idx_seq] <- "seqnames"
  cn[idx_start] <- "start"
  cn[idx_end] <- "end"
  cn[idx_gene] <- "gene_id"
  colnames(region_df) <- cn
  return(region_df)
}
