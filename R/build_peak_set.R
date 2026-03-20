read_peak_file <- function(f) {
  dt <- fread(f, head = FALSE)
  if (ncol(dt) < 3) {
    stop("Peak file has fewer than 3 columns: ", f)
  }
  dt <- dt[, 1:3]
  setnames(dt, c("chrom", "chromStart", "chromEnd"))
  dt[, chrom := as.character(chrom)]
  dt[, chromStart := as.integer(chromStart)]
  dt[, chromEnd := as.integer(chromEnd)]

  GRanges(
    seqnames = dt$chrom,
    ranges = IRanges(start = dt$chromStart + 1L, end = dt$chromEnd),
    strand = "*"
  )
}

build_master_regions <- function(all_gr, min_support = 2) {
  pieces <- disjoin(all_gr, ignore.strand = TRUE)
  ov <- findOverlaps(pieces, all_gr, ignore.strand = TRUE)

  support_n <- integer(length(pieces))
  if (length(ov) > 0) {
    piece2sample <- split(mcols(all_gr)$sample[subjectHits(ov)], queryHits(ov))
    support_n[as.integer(names(piece2sample))] <- vapply(
      piece2sample,
      function(x) length(unique(x)),
      integer(1)
    )
  }

  mcols(pieces)$support_n <- support_n
  pieces_keep <- pieces[mcols(pieces)$support_n >= as.integer(min_support)]

  if (length(pieces_keep) == 0) return(pieces_keep)
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


build_peak_set <- function(peak_path, pair, out_dir = "./", window_size = NULL, min_support = 2) {
  # parameter check
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
  out_bed <- file.path(out_dir, print0(pair, ".bed"))

  # Read peak files
  gr_list <- lapply(seq_along(peak_path), function(i) {
    f <- peak_path[i]
    nm <- tools::file_path_sans_ext(basename(f))
    gr <- read_peak_file(f)
    gr
  })

  gr_list <- gr_list[vapply(gr_list, length, integer(1)) > 0]
  if (length(gr_list) == 0) {
    stop("No valid peaks remain after filtering.")
  }

  # Build master regions
  all_gr <- unlist(GRangesList(gr_list), use.names = FALSE)
  master_regions <- build_master_regions(all_gr = all_gr, min_support = min_support)
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