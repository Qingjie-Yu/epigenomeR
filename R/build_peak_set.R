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
  stopifnot(length(conditions) == length(gr_list))
  
  cond_levels <- unique(conditions)
  pieces <- disjoin(
    unlist(GRangesList(gr_list), use.names = FALSE),
    ignore.strand = TRUE
  )

  mat <- sapply(gr_list, function(gr) {
    countOverlaps(pieces, gr, ignore.strand = TRUE)
  })
  mat <- matrix(mat, nrow = length(pieces))
  
  keep <- apply(mat, 1, function(row) {
    any(vapply(cond_levels, function(cond) {
      sum((row > 0)[conditions == cond]) >= min_support
    }, logical(1)))
  })
  reduce(pieces[keep], ignore.strand = TRUE)
}

tile_region_centered <- function(gr, chr_sizes, window_size, overlap = 0L) {
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
    w   <- e - s + 1L

    n <- ceiling(w / window_size)
    mid <- ceiling((s + e) / 2)

    if (n %% 2L == 1L) {
      anchor <- mid - (window_size %/% 2L)
      offsets <- seq(-(n %/% 2L), n %/% 2L, by = 1L)
    } else {
      anchor  <- mid
      offsets <- seq(-(n %/% 2L), n %/% 2L - 1L, by = 1L)
    }
    starts <- anchor + offsets * step_size
    ends   <- starts + window_size - 1L

    # discard windows that exceed chromosome boundaries
    chr_len <- chr_sizes[[chr]]
    in_bound <- starts >= 1L & ends <= chr_len
    starts <- starts[in_bound]
    ends <- ends[in_bound]
    if (length(starts) == 0L) next

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


build_peak_set <- function(peak_path, pair, conditions, ref_genome = "hg38", out_dir = "./", window_size = NULL, min_support = 2) {
  suppressPackageStartupMessages({
    library(data.table)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
  })
  # parameter check
  if (ref_genome == "hg38") {
    chr_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  } else if (ref_genome == "mm10") {
    chr_sizes <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
  } else {
    stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
  }

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

  # 
  peak_style <- seqlevelsStyle(gr_list[[1]])[1]
  names(chr_sizes) <- GenomeInfoDb::mapSeqlevels(names(chr_sizes), peak_style)

  # Build master regions
  master_regions <- build_master_regions(gr_list = gr_list, conditions = conditions, min_support = min_support)
  if (length(master_regions) == 0) {
    stop("No master regions remain after support filtering.")
  }

  if (is.null(window_size)) {
    peak_set <- master_regions
  } else {
    peak_set <- tile_region_centered(gr = master_regions, chr_sizes = chr_sizes, window_size = window_size)
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