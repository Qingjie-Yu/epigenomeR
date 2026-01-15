# Extract contiguous non-zero coverage blocks from an Rle object, compute block-level AUC, and estimate local background signal using flanking regions.
extract_signal_blocks <- function(cov, chr, chr_size) {
  # Get run values and lengths from Rle
  all_values <- runValue(cov)
  all_lengths <- runLength(cov)
  
  # Calculate genomic coordinates
  all_ends <- cumsum(all_lengths)
  all_starts <- all_ends - all_lengths + 1
  
  # Filter out zero-coverage regions
  keep <- all_values > 0
  if (!any(keep)) return(NULL)
  
  values <- all_values[keep]
  lengths <- all_lengths[keep]
  starts <- all_starts[keep]
  ends <- all_ends[keep]

  n <- length(values)
  if (n == 0) return(NULL)
  
  # Initialize first block
  blocks <- list()
  current_block <- list(
    chr = chr,
    start = starts[1],
    end = ends[1],
    auc = values[1] * lengths[1]
  )
  
  # Merge consecutive positions
  if (n > 1) {
    for (i in 2:n) {
      # Check if consecutive (end+1 == next_start)
      if (ends[i-1] + 1 == starts[i]) {
        # Merge into current block
        current_block$end <- ends[i]
        current_block$auc <- current_block$auc + (values[i] * lengths[i])
      } else {
        # Save current block and start new one
        blocks[[length(blocks) + 1]] <- current_block
        current_block <- list(
          chr = chr,
          start = starts[i],
          end = ends[i],
          auc = values[i] * lengths[i]
        )
      }
    }
  }
  # Add last block
  blocks[[length(blocks) + 1]] <- current_block

  # Helper 1: overlap indexing
  get_overlap_idx <- function(starts, ends, qstart, qend) {
    if (qstart > qend) return(integer(0))
    i1 <- findInterval(qstart - 1L, ends) + 1L
    i2 <- findInterval(qend, starts)
    if (i1 <= i2) i1:i2 else integer(0)
  }

  # Helper 2: vectorized overlap AUC summation
  sum_overlap_auc <- function(starts, ends, values, idx, qstart, qend) {
    if (!length(idx)) return(0)
    ov_start <- pmax.int(starts[idx], qstart)
    ov_end   <- pmin.int(ends[idx], qend)
    ov_len   <- ov_end - ov_start + 1L
    ov_len[ov_len < 0L] <- 0L
    sum(values[idx] * ov_len)
  }

  # Calculate background for each block
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    block_length <- block$end - block$start + 1
    extend_length <- ceiling(4.5 * block_length)
    
    # Define flanking regions (before boundary checking)
    up_start <- block$start - extend_length
    up_end <- block$start - 1
    down_start <- block$end + 1
    down_end <- block$end + extend_length
    
    # Apply chromosome boundaries
    up_start <- max(1, up_start)
    up_end <- min(chr_size, up_end)
    down_start <- max(1, down_start)  
    down_end <- min(chr_size, down_end)

    bg_auc <- block$auc  # Start with block's own AUC
    
    # Add upstream background (only if valid region exists)
    if (up_end >= up_start) {
      up_len <- up_end - up_start + 1
      idx_up <- get_overlap_idx(all_starts, all_ends, up_start, up_end)
      bg_auc <- bg_auc + sum_overlap_auc(all_starts, all_ends, all_values, idx_up, up_start, up_end)
    } else {
      up_len <- 0
    }

    # Add downstream background (only if valid region exists)
    if (down_end >= down_start) {
      down_len <- down_end - down_start + 1
      idx_dn <- get_overlap_idx(all_starts, all_ends, down_start, down_end)
      bg_auc <- bg_auc + sum_overlap_auc(all_starts, all_ends, all_values, idx_dn, down_start, down_end)
    } else {
      down_len <- 0
    } 

    blocks[[i]]$length <- block_length
    blocks[[i]]$bg_auc <- bg_auc
    blocks[[i]]$bg_length <- up_len + block_length + down_len
  }
  
  # Convert to data.frame
  data.table::rbindlist(blocks)
}


# Extract signal blocks from BAM files
peak_calling  <- function(bam_path, out_dir = "./", ref_genome = "hg38", qvalue_cutoff = 0.05, fc_cutoff = 2) {
  # Load Libraries
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(GenomeInfoDb)
    library(Rsamtools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
  })

  # Create folder
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Set up parallel processing
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  cat(sprintf("Using %d CPU cores\n", n_cores))

  if (n_cores > 1) {
    if (.Platform$OS.type == "unix") {
      BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores)
    } else {
      BPPARAM <- BiocParallel::SnowParam(workers = n_cores)
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }

  # Detect chromosome naming style
  bam_header <- scanBamHeader(bam_path[1])
  bam_seqinfo <- Seqinfo(seqnames = names(bam_header[[1]]$targets), seqlengths = bam_header[[1]]$targets)
  bam_style <- seqlevelsStyle(bam_seqinfo)[1]

  # Get reference genome size
  if (ref_genome == "hg38") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
  } else if (ref_genome == "mm10") {
    genome <- BSgenome.Mmusculus.UCSC.mm10
  } else {
    stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
  }
  seqlevelsStyle(genome) <- bam_style
  chr_names <- GenomeInfoDb::standardChromosomes(genome)
  chr_names <- chr_names[!tolower(chr_names) %in% c("mt", "chrm", "m", "mito")]
  chr_sizes <- seqlengths(genome)[chr_names]

  result_list <- lapply(bam_path, function(bam) {
    header <- scanBamHeader(bam)[[1]]
    chr_info <- header$targets
    
    # Extract signal blocks in parallel
    blocks_list <- bplapply(chr_names, function(chr){
      if (!chr %in% names(chr_info)) {
        return(NULL)
      }
      param <- ScanBamParam(which = GRanges(chr, IRanges(1, chr_info[chr])))
      reads <- readGAlignments(bam, param = param)
      if (length(reads) == 0) return(NULL)
      cov <- coverage(reads)[[chr]]
      if (length(cov) == 0) return(NULL)
      extract_signal_blocks(cov, chr, unname(chr_sizes[chr]))
    }, BPPARAM = BPPARAM)

    blocks_list <- blocks_list[!sapply(blocks_list, is.null)]
    if (length(blocks_list) == 0) return(NULL)
    blocks <- data.table::rbindlist(blocks_list)

    blocks$fc <- blocks$auc * blocks$bg_length / (blocks$bg_auc * blocks$length)
    blocks$p_value <- pbinom(
      blocks$auc - 1, 
      size = blocks$bg_auc, 
      prob = blocks$length / blocks$bg_length, 
      lower.tail = FALSE
    )
    blocks$q_value <- p.adjust(blocks$p_value, method = "BH")
    blocks <- blocks[blocks$q_value < qvalue_cutoff & blocks$fc >= fc_cutoff, ]

    blocks$pValue <- -log10(pmax(blocks$p_value, .Machine$double.xmin))
    blocks$qValue <- -log10(pmax(blocks$q_value, .Machine$double.xmin))
    blocks$score <- pmin(as.integer(round(blocks$qValue * 10)), 1000L)
    blocks$name <- paste0(blocks$chr, ":", blocks$start, "-", blocks$end)
    blocks$strand <- "."
    blocks$chromStart <- blocks$start - 1L

    bed <- blocks[, c("chr", "chromStart", "end", "name", "score", "strand", "fc", "pValue", "qValue")]
    colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue")

    pair <- tools::file_path_sans_ext(basename(bam))
    output_file <- file.path(out_dir, paste0(pair, "_peaks.bed"))
    write.table(bed, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  })

  names(result_list) <- tools::file_path_sans_ext(basename(bam_path))
  result_list <- result_list[!sapply(result_list, is.null)]
  invisible(result_list)
}