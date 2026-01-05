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
#   do_qc     : Whether to perform quality control filtering on BAM files.
#   qc_percent: Numeric in (0, 1); percentile threshold for filtering low-count BAM files. Default 0.25.
#   force_chr_coord: When TRUE, region IDs ("pos" column) are always  "CHR_start_end", even if gene_id is available. When FALSE and a non-empty gene_id column exists, gene_id is used as the region identifier.
# Output: Writes a Feather file whose first column is 'pos' and remaining columns are fragment-overlap counts per BAM file. Returns the full output file path (character).

build_count_matrix <- function(bam_path, regions, out_dir = "./", ref_genome = "hg38", sample_name = NULL, do_qc = FALSE, qc_percent = 0.25, force_chr_coord = FALSE) {
    start_time <- Sys.time()

    # Create folder
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    # initiate packages
    suppressPackageStartupMessages({
        library(R.utils)
        library(GenomicAlignments)
        library(GenomicRanges)
        library(GenomeInfoDb)
        library(Rsamtools)
        library(Biostrings)
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(plyranges)
        library(arrow)
        library(preprocessCore)
        library(tibble)
        library(matrixStats)
        library(BiocParallel)
        library(data.table)
        library(rtracklayer)
    })

    num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    register(MulticoreParam(workers = num_cores))

    # Define Regions
    if (is.numeric(regions)) {
        if (length(regions) != 1) {
            stop("Error: numeric 'regions' must be a single bin size.")
        }
        BINSIZE <- regions
        use_custom_region <- FALSE
    } else if (is.character(regions) && all(file.exists(regions))) {
        if (length(regions)!=1) {
            stop("Error: Only one region file is allowed. Please provide a merged file.")
        }
        region_path <- regions
        ext <- tools::file_ext(region_path)
        if (all(ext == "tsv") || all(ext == "txt")) {   
            region_df <- data.table::fread(region_path)
        } else if (all(ext == "csv")){
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

    # qc
    if (do_qc == TRUE) {
        result <- qc(file_path = bam_path, filtered_percentile = qc_percent, save = FALSE)
        vector_crf <- result$filtered_crf
        bamFiles <- bam_path[tools::file_path_sans_ext(basename(bam_path)) %in% vector_crf]
    } else {
        bamFiles <- bam_path
    }
    if (length(bamFiles) == 0) {
        stop("Error: Bam files missing.")
    }

    # Detect chromosome naming style
    bam_header <- scanBamHeader(bamFiles[1])
    bam_seqinfo <- Seqinfo(seqnames = names(bam_header[[1]]$targets), seqlengths = bam_header[[1]]$targets)
    bam_style <- seqlevelsStyle(bam_seqinfo)[1]

    # Post: Use bplapply for parallel chromosome processing
    #       Generate Count Matrix for each chr
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
        for (k in seq_along(bamFiles)) {
            bamFile <- bamFiles[k]
            temp <- readGAlignmentPairs(bamFile)
            locus <- data.frame(first_start = start(temp@first), 
                                first_end = end(temp@first),
                                last_start = start(temp@last), 
                                last_end = end(temp@last))
            overlapCount <- numeric(length(bin))
            if (nrow(locus) == 0) {
                bamContent <- GRanges()
            } else {
                frag_start <- rowMin(as.matrix(locus))
                frag_end <- rowMax(as.matrix(locus))
                bamContent <- makeGRangesFromDataFrame(data.frame(
                    seqnames = as.vector(seqnames(temp)), strand = "*",
                    start = frag_start, end = frag_end))

                overlaps <- findOverlaps(bamContent, bin, ignore.strand = TRUE)
                qh <- queryHits(overlaps)
                sh <- subjectHits(overlaps)

                fragment_starts <- start(bamContent)[qh]
                fragment_ends <- end(bamContent)[qh]
                fragment_lengths <- fragment_ends - fragment_starts + 1

                bin_starts <- start(bin)[sh]
                bin_ends <- end(bin)[sh]

                overlap_starts <- pmax(fragment_starts, bin_starts)
                overlap_ends <- pmin(fragment_ends, bin_ends)
                overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
                proportions <- overlap_lengths / fragment_lengths
                
                if (length(sh) > 0) {
                    summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
                    overlapCount[summed$bin_id] <- summed$x
                }
            }

            bamName <- tools::file_path_sans_ext(basename(bamFile))
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
            refGenome <- BSgenome.Hsapiens.UCSC.hg38
        } else if (ref_genome == "mm10") {
            refGenome <- BSgenome.Mmusculus.UCSC.mm10
        } else {
            stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
        }

        seqlevelsStyle(refGenome) <- bam_style
        chr_list <- GenomeInfoDb::standardChromosomes(refGenome)
        chr_list <- chr_list[!tolower(chr_list) %in% c("mt", "chrm", "m", "mito")]
        chrSizes <- seqlengths(refGenome)[chr_list]

        # Process fixed bins with parallel chromosome processing
        binChriDataframe_list <- bplapply(chr_list, function(chr_i) {
            chrSizei <- chrSizes[chr_i]
            bin <- tileGenome(chrSizei, tilewidth=BINSIZE, cut.last.tile.in.chrom=TRUE)
            binChriDataframe <- as.data.frame(bin)[, c("start", "end")]
            chr_df <- data.frame(CHR = names(chrSizei), stringsAsFactors = FALSE)
            binChriDataframe <- cbind(chr_df, binChriDataframe)

            for (k in seq_along(bamFiles)) {
                bamFile <- bamFiles[k]
                param <- ScanBamParam(which = GRanges(chr_i, IRanges(1, chrSizei)))
                temp <- readGAlignmentPairs(bamFile, param = param)
                locus <- data.frame(first_start = start(temp@first), 
                                    first_end = end(temp@first),
                                    last_start = start(temp@last), 
                                    last_end = end(temp@last))
                overlapCount <- numeric(length(bin))
                if (nrow(locus) == 0) {
                    bamContent <- GRanges()
                } else {
                    frag_start <- rowMin(as.matrix(locus))
                    frag_end <- rowMax(as.matrix(locus))
                    bamContent <- makeGRangesFromDataFrame(data.frame(
                        seqnames = as.vector(seqnames(temp)), strand = "*",
                        start = frag_start, end = frag_end))

                    overlaps <- findOverlaps(bamContent, bin, ignore.strand = TRUE)
                    qh <- queryHits(overlaps)
                    sh <- subjectHits(overlaps)

                    fragment_starts <- start(bamContent)[qh]
                    fragment_ends <- end(bamContent)[qh]
                    fragment_lengths <- fragment_ends - fragment_starts + 1

                    bin_starts <- start(bin)[sh]
                    bin_ends <- end(bin)[sh]

                    overlap_starts <- pmax(fragment_starts, bin_starts)
                    overlap_ends <- pmin(fragment_ends, bin_ends)
                    overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
                    proportions <- overlap_lengths / fragment_lengths
                    
                    if (length(sh) > 0) {
                        summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
                        overlapCount[summed$bin_id] <- summed$x
                    }
                }

                bamName <- tools::file_path_sans_ext(basename(bamFile))
                binChriDataframe[[bamName]] <- overlapCount
            }

            tmp_pos <- binChriDataframe[, c("CHR", "start", "end")]
            tmp_pos$pos <- paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
            pos_df <- data.frame(pos = tmp_pos$pos)
            tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end"))]
            binChriDataframe_final <- cbind(pos_df, tmp_wgc)

            return(binChriDataframe_final)
        })

        binChriDataframe_full <- as.data.frame(do.call(rbind, binChriDataframe_list))
    }

    
    if (is.numeric(regions)) {
        filename <- paste0("Count_Matrix_", BINSIZE)
    } else if (is.character(regions)) {
        prefix <- basename(tools::file_path_sans_ext(regions))
        filename <- paste0("Count_Matrix_", prefix[1])
    }
    
    if (is.null(sample_name)) {
        output_filename <- paste0(filename,'.feather')
    } else {
        output_filename <- paste0(sample_name, '_', filename,'.feather')
    }

  # Report
    preprocess_time <- Sys.time()
    preprocess_time_taken <- round(preprocess_time - start_time, 2)
    cat("preprocess time taken: ", preprocess_time_taken, "\n")

    output_path <- file.path(out_dir, output_filename)
    write_feather(binChriDataframe_full, output_path)

    saving_time <- Sys.time()
    saving_time_taken <- round(saving_time - preprocess_time, 2)
    cat("saving time original taken: ", saving_time_taken, "\n")
    cat("Successfully saved to: ", output_path, "\n")
    return(output_path)
}