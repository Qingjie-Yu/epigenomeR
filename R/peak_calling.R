# Extract blocks for one chromosome from bedGraph intervals
extract_signal_blocks_vec <- function(starts0, ends0, values, chr, chr_size) {
  n <- length(values)
  if (n == 0L) return(NULL)

  # Keep non-zero & non-NA
  keep <- !is.na(values) & values > 0
  if (!any(keep)) return(NULL)

  starts0 <- as.integer(starts0[keep])
  ends0   <- as.integer(ends0[keep])
  values  <- as.numeric(values[keep])

  # Clamp to [0, chr_size]
  chr_size <- as.integer(chr_size)
  starts0 <- pmax.int(0L, starts0)
  ends0   <- pmin.int(chr_size, ends0)

  len0 <- ends0 - starts0
  keep2 <- len0 > 0L
  if (!any(keep2)) return(NULL)

  starts0 <- starts0[keep2]
  ends0   <- ends0[keep2]
  values  <- values[keep2]
  len0    <- len0[keep2]

  # Merge consecutive intervals into blocks: new block if current start != previous end
  new_block <- starts0 != data.table::shift(ends0, fill = starts0[1L])
  block_id <- cumsum(as.integer(new_block))

  dt_tmp <- data.table::data.table(
    block_id = block_id,
    chromStart = starts0,
    chromEnd   = ends0,
    value      = values,
    len        = len0
  )

  blocks <- dt_tmp[, {
    max_val <- max(value)
    is_max  <- value == max_val
    # Check whether all max-value intervals are contiguous in the original
    # interval order (i.e. form a single uninterrupted plateau). If they are
    # scattered with lower-value intervals in between, taking min/max across
    # all of them would incorrectly stretch the summit across the gap.
    max_idx <- which(is_max)
    is_contiguous <- length(max_idx) == 1L || all(diff(max_idx) == 1L)

    if (is_contiguous) {
      s_start <- min(chromStart[is_max])
      s_end   <- max(chromEnd[is_max])
    } else {
      # Scattered ties: fall back to the tie closest to the block's
      # geometric center, rather than spanning across the gap.
      block_mid  <- (min(chromStart) + max(chromEnd)) / 2
      cand_mid   <- (chromStart[max_idx] + chromEnd[max_idx]) / 2
      pick       <- max_idx[which.min(abs(cand_mid - block_mid))]
      s_start    <- chromStart[pick]
      s_end      <- chromEnd[pick]
    }

    .(
      chr          = chr,
      start        = min(chromStart),
      end          = max(chromEnd),
      auc          = sum(value * len),
      length       = sum(len),
      summit_start = s_start,
      summit_end   = s_end
    )
  }, by = block_id]

  # Background window definition (0-based)
  blocks[, ext := as.integer(ceiling(4.5 * length))]
  blocks[, `:=`(
    up_start = pmax.int(0L, start - ext),
    dn_end   = pmin.int(chr_size, end + ext)
  )]

  # Split assignment for robustness
  blocks[, `:=`(
    up_len = pmax.int(0L, start - up_start),
    dn_len = pmax.int(0L, dn_end - end)
  )]
  blocks[, bg_length := length + up_len + dn_len]

  # ---- IRanges overlap AUC ----
  # Convert 0-based half-open [s0, e0) to IRanges 1-based inclusive [s0+1, e0]
  subj_gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = starts0 + 1L, end = ends0),
    strand   = "*"
  )
  subj_val <- values

  auc_for_queries <- function(q0_start, q0_end) {
    out <- numeric(length(q0_start))
    keepq <- q0_start < q0_end
    if (!any(keepq)) return(out)

    qry_gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges   = IRanges::IRanges(start = q0_start[keepq] + 1L, end = q0_end[keepq]),
      strand   = "*"
    )

    hits <- GenomicRanges::findOverlaps(qry_gr, subj_gr, ignore.strand = TRUE)
    if (length(hits) == 0L) return(out)

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # compute overlap width in 1-based inclusive coords
    inter <- IRanges::pintersect(GenomicRanges::ranges(qry_gr)[qh],
                                GenomicRanges::ranges(subj_gr)[sh])
    w <- IRanges::width(inter)

    contrib <- w * subj_val[sh]

    acc <- rowsum(contrib, group = qh, reorder = FALSE)
    out_keep <- numeric(length(qry_gr))
    out_keep[as.integer(rownames(acc))] <- acc[, 1L]

    out[keepq] <- out_keep
    out
  }

  up_auc <- auc_for_queries(blocks$up_start, blocks$start)  # [up_start, start)
  dn_auc <- auc_for_queries(blocks$end, blocks$dn_end)      # [end, dn_end)
  blocks[, bg_auc := auc + up_auc + dn_auc]

  blocks[, c("ext", "up_start", "dn_end", "up_len", "dn_len", "block_id") := NULL]
  blocks[]
}

# Write a data.table of peak blocks (possibly empty) to narrowPeak format
write_narrowpeak <- function(blocks, pair, out_dir) {
  output_file <- file.path(out_dir, paste0(pair, "_peaks.narrowPeak"))
  n_peaks <- nrow(blocks)

  if (n_peaks == 0L) {
    bed <- data.table::data.table(
      chrom = character(0), chromStart = integer(0), chromEnd = integer(0),
      name = character(0), score = integer(0), strand = character(0),
      signalValue = numeric(0), pValue = numeric(0), qValue = numeric(0),
      peak = integer(0), length = integer(0), auc = numeric(0), cov = numeric(0)
    )
  } else {
    blocks$pValue <- -log10(pmax(blocks$p_value, .Machine$double.xmin))
    blocks$qValue <- -log10(pmax(blocks$q_value, .Machine$double.xmin))
    blocks$score  <- pmin(as.integer(round(blocks$qValue * 10)), 1000L)
    blocks$name   <- paste0(blocks$chr, ":", blocks$start, "-", blocks$end)
    blocks$strand <- "."
    blocks$chromStart <- blocks$start
    blocks$peak   <- as.integer(round((blocks$summit_start + blocks$summit_end) / 2 - blocks$start))

    bed <- blocks[, c("chr", "chromStart", "end", "name", "score", "strand",
                       "fc", "pValue", "qValue", "peak", "length", "auc", "cov")]
    data.table::setnames(bed, c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak", "length", "auc", "cov"))
  }

  utils::write.table(bed, file = output_file, sep = "\t", quote = FALSE,  row.names = FALSE, col.names = FALSE)

  message(sprintf("[%s] %d peak%s called\n", pair, n_peaks, if (n_peaks == 1L) "" else "s"))

  output_file
}

# Extract signal blocks from BEDGRAPH files
peak_calling <- function(bedgraph_path, out_dir = "./", ref_genome = "hg38", min_cov = 2, qvalue_cutoff = 0.05, fc_cutoff = 2, auc_top_pct = 0.1) {
  suppressPackageStartupMessages({
    library(GenomeInfoDb)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
  })
  # Create output directory
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Parameter Check
  if (min_cov < 0) {
    stop("min_cov must be non-negative.")
  }
  if (auc_top_pct <= 0 || auc_top_pct > 1) {
    stop("auc_top_pct must be between 0 (exclusive) and 1 (inclusive).")
  }

  # Detect chromosome naming style from first file
  dt_preview <- data.table::fread(bedgraph_path[[1]], header = FALSE, sep = "\t", nrows = 200, showProgress = FALSE)
  bed_style <- if (any(grepl("^chr", dt_preview$V1))) "UCSC" else "NCBI"

  # Get reference genome object and sizes
  genome <- switch(
    ref_genome,
    "hg38" = BSgenome.Hsapiens.UCSC.hg38,
    "mm10" = BSgenome.Mmusculus.UCSC.mm10,
    stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
  )
  GenomeInfoDb::seqlevelsStyle(genome) <- bed_style
  chr_names <- GenomeInfoDb::standardChromosomes(genome)
  chr_names <- chr_names[!tolower(chr_names) %in% c("mt", "chrm", "m", "mito")]
  chr_sizes <- GenomeInfoDb::seqlengths(genome)[chr_names]

  # Parallel over files
  BPPARAM <- get_BPPARAM(backend = "multicore")

  # helper function
  is_integerish <- function(x, n_check = 200L, tol = 1e-6) {
    n <- min(as.integer(n_check), length(x))
    xx <- x[seq_len(n)]
    all(abs(xx - round(xx)) < tol)
  }

  result_list <- BiocParallel::bplapply(bedgraph_path, function(bg) {
    pair <- tools::file_path_sans_ext(basename(bg))
    # Read bedGraph
    dt_all <- data.table::fread(bg, header = FALSE, sep = "\t", select = 1:4, showProgress = FALSE)
    if (ncol(dt_all) < 4) {
      message("Skipping ", pair, ": invalid bedGraph format.")
      return(write_narrowpeak(data.table::data.table(), pair, out_dir))
    }

    dt_all <- dt_all[, 1:4]
    data.table::setnames(dt_all, c("chrom", "chromStart", "chromEnd", "value"))
    dt_all[, chromStart := as.integer(chromStart)]
    dt_all[, chromEnd   := as.integer(chromEnd)]
    dt_all[, value      := as.numeric(value)]

    # Ensure grouped by chrom for run slicing
    data.table::setorder(dt_all, chrom, chromStart, chromEnd)

    # Build chrom run ranges once; no dt_chr subsetting
    chrom_vec <- dt_all$chrom
    run_id <- data.table::rleid(chrom_vec)
    runs <- dt_all[, .(chrom = chrom[1L], i = .I[1L], j = .I[.N]), by = run_id]
    runs <- runs[chrom %in% chr_names]
    if (nrow(runs) == 0L) {
      return(write_narrowpeak(data.table::data.table(), pair, out_dir))
    }

    starts_all <- dt_all$chromStart
    ends_all   <- dt_all$chromEnd
    values_all <- dt_all$value

    blocks_list <- lapply(seq_len(nrow(runs)), function(k) {
      chr <- runs$chrom[k]
      chr_size <- unname(chr_sizes[chr])
      if (is.na(chr_size) || chr_size <= 0) return(NULL)

      i <- runs$i[k]
      j <- runs$j[k]

      extract_signal_blocks_vec(
        starts0 = starts_all[i:j],
        ends0   = ends_all[i:j],
        values  = values_all[i:j],
        chr     = chr,
        chr_size = chr_size
      )
    })

    blocks_list <- blocks_list[!vapply(blocks_list, is.null, logical(1))]
    if (!length(blocks_list)) {
      return(write_narrowpeak(data.table::data.table(), pair, out_dir))
    }

    blocks <- data.table::rbindlist(blocks_list)

    # FC
    blocks$fc <- blocks$auc * blocks$bg_length / (blocks$bg_auc * blocks$length)

    # Minimum mean coverage pre-filter: drop "wide-and-shallow" blocks
    blocks$cov <- blocks$auc / blocks$length
    blocks <- blocks[blocks$cov >= min_cov, ]
    if (nrow(blocks) == 0L) {
      warning("No blocks remaining after min_cov (", min_cov, ") filter for: ", pair)
      return(write_narrowpeak(blocks, pair, out_dir))
    }

    # P-value
    if (is_integerish(blocks$auc) && is_integerish(blocks$bg_auc)) {
      blocks$p_value <- stats::pbinom(
        round(blocks$auc) - 1,
        size = round(blocks$bg_auc),
        prob = blocks$length / blocks$bg_length,
        lower.tail = FALSE
      )
    } else {
      k <- pmax(floor(blocks$auc), 0)
      lambda <- (blocks$bg_auc / blocks$bg_length) * blocks$length
      blocks$p_value <- stats::ppois(k - 1, lambda = lambda, lower.tail = FALSE)
    }

    # Q-value
    blocks$q_value <- stats::p.adjust(blocks$p_value, method = "BH")

    # Filter
    blocks <- blocks[blocks$q_value < qvalue_cutoff & blocks$fc >= fc_cutoff, ]
    if (nrow(blocks) == 0L) {
      warning("No blocks passed qvalue_cutoff (", qvalue_cutoff,
              ") and fc_cutoff (", fc_cutoff, ") for: ", pair)
      return(write_narrowpeak(blocks, pair, out_dir))
    }

    # AUC percentile AFTER statistical testing
    if (auc_top_pct < 1.0) {
      auc_thresh <- stats::quantile(blocks$auc, probs = 1 - auc_top_pct)
      blocks <- blocks[blocks$auc >= auc_thresh, ]
      if (nrow(blocks) == 0L) {
        warning("No blocks remaining after post-filter AUC top ",
                auc_top_pct * 100, "% for: ", pair)
        return(write_narrowpeak(blocks, pair, out_dir))
      }
    }

    write_narrowpeak(blocks, pair, out_dir)
  }, BPPARAM = BPPARAM)

  names(result_list) <- tools::file_path_sans_ext(basename(bedgraph_path))
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  invisible(result_list)
}
