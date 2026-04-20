is_bam_paired <- function(bam_path) {
  param <- ScanBamParam(what = "flag")
  flags <- scanBam(bam_path[1], param = param)[[1]]$flag
  if (length(flags) == 0) return(FALSE)
  any(bitwAnd(head(flags, 1000), 0x1) != 0)
}

# Build Count Matrix Function
# Post: Build a fragment-overlap count matrix from paired-end BAM files over user-specified genomic regions and save it as a Feather file.
# Parameter:
#   bam_path  : Character vector of BAM file paths. Each BAM file becomes one column in the output matrix.
#   regions   : Either
#               - a single integer (bin size in bp), e.g. regions = 5000, which tiles the genome into fixed-size bins; or
#               - a single file path to a BED / TSV / TXT / CSV file containing custom genomic regions.
#   out_dir  : Directory where the output Feather file will be written. Default "./".
#   ref_genome: Reference genome used when regions is numeric. One of "hg38" or "mm10". Ignored when custom regions are provided.
#   sample_name: Optional character. If provided, it is prepended to the output filename.
#   force_chr_coord: When TRUE, region IDs ("pos" column) are always  "CHR_start_end", even if gene_id is available. When FALSE and a non-empty gene_id column exists, gene_id is used as the region identifier.
#   by_single_crf  : When TRUE, after building the pair-level matrix (columns named "CRF1-CRF2"), aggregate counts by individual CRF: each output column is the sum of all pair columns that contain that CRF. Default FALSE.
# Output: Writes a Feather file whose first column is 'pos' and remaining columns are fragment-overlap counts per BAM file. Returns the full output file path (character).

build_count_matrix <- function(bam_path, regions, out_dir = "./", ref_genome = "hg38", sample_name = NULL, force_chr_coord = FALSE, by_single_crf = FALSE) {
  # Load Libraries
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(Rsamtools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(arrow)
    library(matrixStats)
    library(rtracklayer)
  })

  # Set up parallel processing
  BPPARAM <- get_BPPARAM()
  
  # Create folder
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Define Regions
  if (is.numeric(regions)) {
    if (length(regions) != 1) {
      stop("Error: numeric 'regions' must be a single bin size.")
    }
    BINSIZE <- regions
    use_custom_region <- FALSE
  } else if (is.character(regions) && all(file.exists(regions))) {
    if (length(regions) != 1) {
      stop("Error: Only one region file is allowed. Please provide a merged file.")
    }
    region_path <- regions
    ext <- tools::file_ext(region_path)
    if (all(ext == "tsv") || all(ext == "txt")) {
      region_df <- data.table::fread(region_path)
    } else if (all(ext == "csv")) {
      region_df <- read.csv(region_path, row.names = NULL)
    } else if (all(ext == "bed")) {
      region_df <- NULL
    } else {
      stop("Error: Invalid region file format. Only support .csv, .tsv, .txt, .bed")
    }
    use_custom_region <- TRUE
  } else {
    stop("Error: Custom regions must be provided. 'regions' argument is missing or invalid.")
  }

  # Detect chromosome naming style
  bam_path <- sort(bam_path) # sort by name
  bam_header <- scanBamHeader(bam_path[1])
  bam_seqinfo <- Seqinfo(seqnames = names(bam_header[[1]]$targets), seqlengths = bam_header[[1]]$targets)
  bam_style <- seqlevelsStyle(bam_seqinfo)[1]

  # Detect single or paired
  is_paired <- is_bam_paired(bam_path)

  # Generate Count Matrix for each chr
  # Process custom regions (no chromosome loop needed)
  if (use_custom_region) {
    # Handle CSV/GTF file
    if (!is.null(region_df)) {
      region_df <- fix_region_colnames(region_df)
      bin <- GRanges(seqnames = region_df$seqnames, ranges = IRanges(start = region_df$start, end = region_df$end))
      seqlevelsStyle(bin) <- bam_style
      binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
      colnames(binChriDataframe)[1] <- "CHR"
      binChriDataframe$CHR <- as.character(binChriDataframe$CHR)

      if ("gene_id" %in% colnames(region_df)) {
        mcols(bin)$gene_id <- region_df$gene_id
        binChriDataframe$gene_id <- region_df$gene_id
      }
    }
    # Handle BED file
    else {
      bed_data <- read.table(region_path, sep = "\t", stringsAsFactors = FALSE)
      bin <- GRanges(seqnames = bed_data$V1, ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3))
      if (ncol(bed_data) >= 4) {
        mcols(bin)$gene_id <- bed_data$V4
      }
      seqlevelsStyle(bin) <- bam_style
      binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
      colnames(binChriDataframe)[1] <- "CHR"
      binChriDataframe$CHR <- as.character(binChriDataframe$CHR)

      if (!is.null(mcols(bin)) && "gene_id" %in% colnames(mcols(bin))) {
        binChriDataframe$gene_id <- mcols(bin)$gene_id
      }
    }

    # Process BAM files for custom regions
    for (k in seq_along(bam_path)) {
      bam <- bam_path[k]
      overlapCount <- numeric(length(bin))

      if (is_bam_paired(bam)) {
        temp <- readGAlignmentPairs(bam)
        if (length(temp) == 0) {
          bamContent <- GRanges()
        } else {
          locus <- data.frame(
            first_start = start(temp@first),
            first_end   = end(temp@first),
            last_start  = start(temp@last),
            last_end    = end(temp@last)
          )
          frag_start <- rowMin(as.matrix(locus))
          frag_end   <- rowMax(as.matrix(locus))
          bamContent <- makeGRangesFromDataFrame(data.frame(
            seqnames = as.vector(seqnames(temp)), strand = "*",
            start = frag_start, end = frag_end
          ))
        }
      } else {
        temp <- readGAlignments(bam)
        if (length(temp) == 0) {
          bamContent <- GRanges()
        } else {
          bamContent <- GRanges(
            seqnames = seqnames(temp),
            ranges   = IRanges(start = start(temp), end = end(temp)),
            strand   = "*"
          )
        }
      }

      if (length(bamContent) > 0) {
        overlaps <- findOverlaps(bamContent, bin, ignore.strand = TRUE)
        qh <- queryHits(overlaps)
        sh <- subjectHits(overlaps)

        fragment_starts <- start(bamContent)[qh]
        fragment_ends   <- end(bamContent)[qh]
        fragment_lengths <- fragment_ends - fragment_starts + 1

        bin_starts <- start(bin)[sh]
        bin_ends   <- end(bin)[sh]

        overlap_starts  <- pmax(fragment_starts, bin_starts)
        overlap_ends    <- pmin(fragment_ends, bin_ends)
        overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
        proportions     <- overlap_lengths / fragment_lengths

        if (length(sh) > 0) {
          summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
          overlapCount[summed$bin_id] <- summed$x
        }
      }

      bamName <- tools::file_path_sans_ext(basename(bam))
      binChriDataframe[[bamName]] <- overlapCount
    }
    # Format output
    has_gene_id <- "gene_id" %in% colnames(binChriDataframe) &&
      all(!is.na(binChriDataframe$gene_id)) &&
      all(binChriDataframe$gene_id != "")
    if (has_gene_id && !force_chr_coord) {
      pos_vec <- binChriDataframe$gene_id
      drop_cols <- c("CHR", "start", "end", "gene_id")
    } else {
      pos_vec <- paste0(binChriDataframe$CHR, "_", binChriDataframe$start, "_", binChriDataframe$end)
      drop_cols <- c("CHR", "start", "end")
    }
    pos_df <- data.frame(pos = pos_vec)
    tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% drop_cols), drop = FALSE]
    binChriDataframe_full <- cbind(pos_df, tmp_wgc)
  } else {
    # Get reference genome size
    if (ref_genome == "hg38") {
      chrSizes0 <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    } else if (ref_genome == "mm10") {
      chrSizes0 <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
    } else {
      stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
    }

    ref_seqinfo <- Seqinfo(
      seqnames = names(chrSizes0),
      seqlengths = chrSizes0
    )
    seqlevelsStyle(ref_seqinfo) <- bam_style

    chr_list <- GenomeInfoDb::standardChromosomes(ref_seqinfo)
    chr_list <- chr_list[!tolower(chr_list) %in% c("mt", "chrm", "m", "mito")]
    chrSizes <- seqlengths(ref_seqinfo)[chr_list]

    # Process fixed bins with parallel chromosome processing
    binChriDataframe_list <- bplapply(chr_list, function(chr_i) {
      chrSizei <- chrSizes[chr_i]
      bin <- tileGenome(chrSizei, tilewidth = BINSIZE, cut.last.tile.in.chrom = TRUE)
      binChriDataframe <- as.data.frame(bin)[, c("start", "end")]
      chr_df <- data.frame(CHR = names(chrSizei), stringsAsFactors = FALSE)
      binChriDataframe <- cbind(chr_df, binChriDataframe)

      for (k in seq_along(bam_path)) {
        bam <- bam_path[k]
        param <- ScanBamParam(which = GRanges(chr_i, IRanges(1, chrSizei)))
        overlapCount <- numeric(length(bin))

        if (is_paired) {
          temp <- readGAlignmentPairs(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else {
            locus <- data.frame(
              first_start = start(temp@first),
              first_end   = end(temp@first),
              last_start  = start(temp@last),
              last_end    = end(temp@last)
            )
            frag_start <- rowMin(as.matrix(locus))
            frag_end   <- rowMax(as.matrix(locus))
            bamContent <- makeGRangesFromDataFrame(data.frame(
              seqnames = as.vector(seqnames(temp)), strand = "*",
              start = frag_start, end = frag_end
            ))
          }
        } else {
          temp <- readGAlignments(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else {
            bamContent <- GRanges(
              seqnames = seqnames(temp),
              ranges   = IRanges(start = start(temp), end = end(temp)),
              strand   = "*"
            )
          }
        }

        if (length(bamContent) > 0) {
          overlaps <- findOverlaps(bamContent, bin, ignore.strand = TRUE)
          qh <- queryHits(overlaps)
          sh <- subjectHits(overlaps)

          fragment_starts  <- start(bamContent)[qh]
          fragment_ends    <- end(bamContent)[qh]
          fragment_lengths <- fragment_ends - fragment_starts + 1

          bin_starts <- start(bin)[sh]
          bin_ends   <- end(bin)[sh]

          overlap_starts  <- pmax(fragment_starts, bin_starts)
          overlap_ends    <- pmin(fragment_ends, bin_ends)
          overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
          proportions     <- overlap_lengths / fragment_lengths

          if (length(sh) > 0) {
            summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
            overlapCount[summed$bin_id] <- summed$x
          }
        }

        bamName <- tools::file_path_sans_ext(basename(bam))
        binChriDataframe[[bamName]] <- overlapCount
      }

      tmp_pos <- binChriDataframe[, c("CHR", "start", "end")]
      tmp_pos$pos <- paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
      pos_df <- data.frame(pos = tmp_pos$pos)
      tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end"))]
      binChriDataframe_final <- cbind(pos_df, tmp_wgc)

      return(binChriDataframe_final)
    }, BPPARAM = BPPARAM)

    binChriDataframe_full <- as.data.frame(data.table::rbindlist(binChriDataframe_list))
  }


  # Aggregate by single CRF: sum all pair columns that contain each CRF name
  if (by_single_crf) {
    pos_col   <- binChriDataframe_full[, "pos", drop = FALSE]
    count_cols <- binChriDataframe_full[, -1, drop = FALSE]
    pair_names <- colnames(count_cols)
    crf_names  <- unique(unlist(strsplit(pair_names, "-")))
    crf_mat <- sapply(crf_names, function(crf) {
      matched <- pair_names[sapply(strsplit(pair_names, "-"), function(parts) crf %in% parts)]
      rowSums(count_cols[, matched, drop = FALSE])
    })
    binChriDataframe_full <- cbind(pos_col, as.data.frame(crf_mat))
  }

  if (is.numeric(regions)) {
    filename <- paste0("Count_Matrix_", BINSIZE)
  } else if (is.character(regions)) {
    prefix <- basename(tools::file_path_sans_ext(regions))
    filename <- paste0("Count_Matrix_", prefix[1])
  }

  if (is.null(sample_name)) {
    output_filename <- paste0(filename, ".feather")
  } else {
    output_filename <- paste0(sample_name, "_", filename, ".feather")
  }

  # Report
  output_path <- file.path(out_dir, output_filename)
  write_feather(binChriDataframe_full, output_path)
  cat("Successfully saved to: ", output_path, "\n")
  return(output_path)
}
