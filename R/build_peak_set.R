read_peak_file <- function(f) {
  dt <- data.table::fread(f, head = FALSE)
  if (ncol(dt) < 3) {
    stop("Peak file has fewer than 3 columns: ", f)
  }
  dt <- dt[, 1:3]
  data.table::setnames(dt, c("chrom", "chromStart", "chromEnd"))
  dt[, chrom := as.character(chrom)]
  dt[, chromStart := as.integer(chromStart)]
  dt[, chromEnd := as.integer(chromEnd)]

  GRanges(
    seqnames = dt$chrom,
    ranges = IRanges::IRanges(start = dt$chromStart + 1L, end = dt$chromEnd),
    strand = "*"
  )
}

build_master_regions <- function(gr_list, conditions, min_support = 2) {
  sample_to_cond <- setNames(conditions, as.character(seq_along(conditions)))
  cond_levels    <- unique(conditions)

  all_gr <- unlist(GRangesList(gr_list), use.names = FALSE) 
  sample_idx <- rep(seq_along(gr_list),
                    vapply(gr_list, length, integer(1)))

  pieces <- disjoin(all_gr, ignore.strand = TRUE)
  ov     <- findOverlaps(pieces, all_gr, ignore.strand = TRUE)

  piece2sample <- rep(list(integer(0)), length(pieces))
  if (length(ov) > 0) {
    piece2sample_raw <- split(
      sample_idx[subjectHits(ov)],
      queryHits(ov)
    )
    piece2sample[as.integer(names(piece2sample_raw))] <- piece2sample_raw
  }

  keep <- vapply(piece2sample, function(samp_ids) {
    if (length(samp_ids) == 0L) return(FALSE)
    unique_samp <- as.character(unique(samp_ids))
    any(vapply(cond_levels, function(cond) {
      sum(sample_to_cond[unique_samp] == cond, na.rm = TRUE) >= min_support
    }, logical(1)))
  }, logical(1))

  pieces_keep <- pieces[keep]
  if (length(pieces_keep) == 0L) return(pieces_keep)
  reduce(pieces_keep, ignore.strand = TRUE)
}

tile_region_with_tail <- function(gr, window_size, overlap = 0L) {
  # overlap: number of bases shared between adjacent windows (0 = non-overlapping)
  if (overlap < 0L || overlap >= window_size) {
    stop("'overlap' must be >= 0 and < window_size.")
  }
  step_size <- window_size - overlap

  out_list <- vector("list", length(gr))

  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr)[i])
    s <- start(gr)[i]
    e <- end(gr)[i]

    if (width(gr)[i] < window_size) next

    starts <- seq(from = s, to = e - window_size + 1L, by = step_size)

    # append tail window if the last regular window doesn't reach the region end
    tail_start <- e - window_size + 1L
    if (tail_start > tail(starts, 1L)) {
      starts <- c(starts, tail_start)
    }
    ends <- starts + window_size - 1L

    out_list[[i]] <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = starts, end = ends),
      strand = "*"
    )
  }

  out_list <- Filter(Negate(is.null), out_list)
  if (length(out_list) == 0L) return(gr[0])
  unlist(GRangesList(out_list), use.names = FALSE)
}


build_peak_set <- function(peak_path, pair, conditions, out_dir = "./", window_size = NULL, min_support = 2) {
  suppressPackageStartupMessages({
    library(data.table)
  })
  # parameter check
  if (length(peak_path) != length(conditions)) {
    stop("peak_path and conditions must have the same length.")
  }
  cond_table <- table(conditions)
  if (any(cond_table < min_support)) {
    warning("Some conditions have fewer samples than min_support (",
            min_support, "): ",
            paste(names(cond_table)[cond_table < min_support], collapse = ", "))
  }

  if (!is.null(window_size)) {
    if (!is.numeric(window_size) || window_size <= 0) {
      stop("'window_size' must be a positive integer.")
    } else {
      message("Tiling mode: window_size = ", window_size)
    }
  } else {
    message("Peak union mode: using master regions directly (no tiling).")
  }

  # build output path
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_bed <- file.path(out_dir, paste0(pair, ".bed"))

  # Read peak files
  gr_list <- lapply(seq_along(peak_path), function(i) {
    read_peak_file(peak_path[i])
  })
  valid <- vapply(gr_list, length, integer(1)) > 0
  gr_list <- gr_list[valid]
  conditions <- conditions[valid]
  if (length(gr_list) == 0) {
    stop("No valid peaks remain after filtering.")
  }

  # Build master regions
  master_regions <- build_master_regions(gr_list = gr_list, conditions = conditions, min_support = min_support)
  if (length(master_regions) == 0) {
    stop("No master regions remain after support filtering.")
  }

  if (is.null(window_size)) {
    peak_set <- master_regions
  } else {
    peak_set <- tile_region_with_tail(gr = master_regions, window_size = window_size)
  }
  peak_set <- sort(peak_set, ignore.strand = TRUE)

  # Write
  bed_dt <- data.table(
    chrom = as.character(seqnames(peak_set)),
    chromStart = start(peak_set) - 1L,
    chromEnd = end(peak_set)
  )
  fwrite(bed_dt, file = out_bed, sep = "\t", col.names = FALSE)
  message("Output peak set written to: ", out_bed)
  message("Final number of regions/windows: ", nrow(bed_dt))
  invisible(out_bed)
}