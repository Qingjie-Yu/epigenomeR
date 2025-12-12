# Identify the nearest gene for each genomic interval.
# 
# This function takes a set of query genomic regions and a reference set of gene
# annotations (both as GRanges objects), harmonizes their chromosome naming
# styles, and for each query region determines the nearest gene on the same
# chromosome. If the query midpoint overlaps one or more genes, the shortest
# overlapping gene is selected. Otherwise, the function evaluates the nearest
# upstream and downstream genes and chooses the closest gene body; ties are
# resolved by selecting the shorter gene. The output returns, for each query
# region, the gene name, gene length, and directed distance to the gene's TSS.
#
# Parameters:
#   query_gr: GRanges object of query regions.
#   genes_gr: GRanges object of gene annotations, with gene names stored in names().
#
# Returns:
#   A data.frame with one row per query region containing:
#     - chr: chromosome of the query region
#     - start, end: genomic coordinates of the query region
#     - nearest_gene: name of the nearest gene
#     - gene_length: width of the selected gene
#     - tss_distance: signed distance to the gene TSS 
#                     (negative = upstream, positive = downstream)

get_nearest_gene <- function(query_gr, genes_gr) {
    query_mid <- resize(query_gr, width = 1, fix = "center")
    overlaps <- findOverlaps(query_mid, genes_gr)
    n_query <- length(query_gr)
    nearest_gene <- rep(NA_character_, n_query)
    gene_length <- rep(NA_integer_, n_query)
    tss_distance <- rep(NA_integer_, n_query)

    # Handle overlapping genes - select shortest gene
    if (length(overlaps) > 0) {
        query_hits <- queryHits(overlaps)
        subject_hits <- subjectHits(overlaps)

        overlap_df <- data.frame(
            query_idx = query_hits,
            gene_idx = subject_hits,
            gene_width = width(genes_gr[subject_hits]),
            stringsAsFactors = FALSE
        )

        best_overlaps <- overlap_df %>% 
            group_by(query_idx) %>% 
            slice_min(gene_width, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        for (i in seq_len(nrow(best_overlaps))) {
            q_idx <- best_overlaps$query_idx[i]
            g_idx <- best_overlaps$gene_idx[i]
            best_gene <- genes_gr[g_idx]
            q_mid_pos <- start(query_mid[q_idx])

            nearest_gene[q_idx] <- names(best_gene)
            gene_length[q_idx] <- width(best_gene)

            # TSS distance: negative = upstream, positive = downstream
            if (as.character(strand(best_gene)) == "-") {
                tss_distance[q_idx] <- q_mid_pos - end(best_gene)
            } else {
                tss_distance[q_idx] <- q_mid_pos - start(best_gene)
            }
        }
    }

    no_overlap_idx <- which(is.na(nearest_gene))

    # Handle non-overlapping genes
    if (length(no_overlap_idx) > 0) {
        query_chr <- as.character(seqnames(query_mid))

        for (chr in unique(query_chr[no_overlap_idx])) {
            chr_query_idx <- no_overlap_idx[query_chr[no_overlap_idx] == chr]
            chr_genes <- genes_gr[seqnames(genes_gr) == chr]
            
            if (length(chr_genes) == 0) next
            
            chr_query_mid <- query_mid[chr_query_idx]
            nearest_info <- distanceToNearest(chr_query_mid, chr_genes)
            
            if (length(nearest_info) > 0) {
                # Handle ties by selecting shorter gene
                dist_df <- data.frame(
                    q_local_idx = queryHits(nearest_info),
                    g_local_idx = subjectHits(nearest_info),
                    distance = mcols(nearest_info)$distance,
                    gene_width = width(chr_genes[subjectHits(nearest_info)])
                )
                
                best_matches <- dist_df %>%
                    group_by(q_local_idx) %>%
                    filter(distance == min(distance)) %>%
                    slice_min(gene_width, n = 1, with_ties = FALSE) %>%
                    ungroup()
                
                for (i in seq_len(nrow(best_matches))) {
                    q_local_idx <- best_matches$q_local_idx[i]
                    g_local_idx <- best_matches$g_local_idx[i]
                    q_global_idx <- chr_query_idx[q_local_idx]
                    
                    best_gene <- chr_genes[g_local_idx]
                    q_mid_pos <- start(chr_query_mid[q_local_idx])
                    
                    nearest_gene[q_global_idx] <- names(best_gene)
                    gene_length[q_global_idx] <- width(best_gene)
                    
                    # TSS distance: negative = upstream, positive = downstream
                    if (as.character(strand(best_gene)) == "-") {
                        tss_distance[q_global_idx] <- q_mid_pos - end(best_gene)
                    } else {
                        tss_distance[q_global_idx] <- q_mid_pos - start(best_gene)
                    }
                }
            }
        }
    } 

    result_final <- data.frame(
        chr = as.character(seqnames(query_gr)),
        start = start(query_gr),
        end = end(query_gr),
        nearest_gene = nearest_gene,
        gene_length = gene_length,
        tss_distance = tss_distance,
        stringsAsFactors = FALSE
    )
    result_final <- result_final[!is.na(result_final$nearest_gene), ]
    return(result_final)
}


# Generate matched control regions
#
# This function generates background control regions by randomly selecting genes
# with similar length to the nearest gene of each target peak, then placing
# control regions at the same TSS distance on the selected genes. This approach
# preserves the relationship between peaks and gene structure while avoiding
# biases from using the same genomic neighborhoods.
#
# Parameters:
#   target_gr: GRanges object of target peak regions
#   ref_genome: Reference genome version ("hg38" or "mm10")
#   ref_source: Gene annotation source used to define gene models and TSS
#               coordinates. Supported options are:
#               - "knownGene": Uses the UCSC knownGene annotation obtained via
#                 the Bioconductor package TxDb.Hsapiens.UCSC.hg38.knownGene
#               - "GENCODE": Uses GENCODE gene annotations (GENCODE v49).
#   style: Chromosome naming style (default: "UCSC")
#   n_rep: Number of control regions to generate per target peak
#   regions: Width of control regions in base pairs
#   seed: Random seed for reproducibility
#   length_tolerance: Initial tolerance for gene length matching (±20% default)
#                     Automatically relaxed up to ±100% if insufficient candidates
#
# Returns:
#   A GRanges object containing all control regions

get_matched_control <- function(target_gr, ref_genome = "hg38", ref_source = "knownGene", style = "UCSC", n_rep = 1, regions = 800, seed = 42, length_tolerance = 0.2) {
    set.seed(seed)

    if (!ref_genome %in% names(genome_config)) {
        stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
    }
    
    genome_config <- list(
    hg38 = list(
        bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
        txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        gencode_file = "GENCODE_v49_hg38.rds"
    ),
    mm10 = list(
        bsgenome = "BSgenome.Mmusculus.UCSC.mm10",
        txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        gencode_file = "GENCODE_vM35_mm10.rds"
    )
    )
    config <- genome_config[[ref_genome]]

    library(config$bsgenome, character.only = TRUE, quietly = TRUE)
    bs <- getExportedValue(config$bsgenome, config$bsgenome)

    if (ref_source == "knownGene") {
    library(config$txdb, character.only = TRUE, quietly = TRUE)
    txdb <- getExportedValue(config$txdb, config$txdb)
    genes_gr <- suppressMessages(genes(txdb))
    } else {
    genes_gr <- readRDS(system.file("extdata", config$gencode_file, package = "epigenomeR"))
    genes_gr <- genes_gr[genes_gr$type == "gene"]
    mcols(genes_gr) <- NULL
    }

    seqlevelsStyle(bs) <- style
    chr_list <- GenomeInfoDb::standardChromosomes(bs)
    chr_list <- chr_list[!tolower(chr_list) %in% c("mt", "chrm", "m", "mito")]
    chr_sizes <- seqlengths(bs)[chr_list]

    seqlevelsStyle(genes_gr) <- style
    genes_gr <- genes_gr[seqnames(genes_gr) %in% chr_list]
    seqlevels(genes_gr) <- chr_list
    genes_gr <- sort(genes_gr)

    gene_lengths <- width(genes_gr)
    gene_names <- names(genes_gr)
    gene_chrs <- as.character(seqnames(genes_gr))
    gene_strands <- as.character(strand(genes_gr))
    gene_starts <- start(genes_gr)
    gene_ends <- end(genes_gr)

    target_with_genes <- get_nearest_gene(target_gr, genes_gr)
    n_targets <- nrow(target_with_genes)
    half_width <- floor(regions / 2)

    control_list <- vector("list", n_targets * n_rep)
    control_idx <- 1

    for (i in seq_len(n_targets)) {
        orig_gene <- target_with_genes$nearest_gene[i]
        orig_gene_length <- target_with_genes$gene_length[i]
        tss_distance <- target_with_genes$tss_distance[i]

        tol <- length_tolerance
        lower_bound <- orig_gene_length * (1 - tol)
        upper_bound <- orig_gene_length * (1 + tol)
        candidate_mask <- (gene_lengths >= lower_bound) & 
                          (gene_lengths <= upper_bound) & 
                          (gene_names != orig_gene)

        while (sum(candidate_mask) < n_rep * 2 && tol <= 0.9) {
            tol <- tol + 0.1
            lower_bound <- orig_gene_length * (1 - tol)
            upper_bound <- orig_gene_length * (1 + tol)
            candidate_mask <- (gene_lengths >= lower_bound) & 
                              (gene_lengths <= upper_bound) & 
                              (gene_names != orig_gene)
        }

        if (sum(candidate_mask) < n_rep) {
            candidate_mask <- (gene_names != orig_gene)
        }
        
        candidate_indices <- which(candidate_mask)
        
        if (length(candidate_indices) == 0) {
            warning(sprintf("No candidate genes found for target %d", i))
            next
        }

        n_valid <- 0
        max_attempts <- min(length(candidate_indices), n_rep * 5)
        sampled_indices <- sample(candidate_indices, size = max_attempts, replace = FALSE)
        
        for (idx in sampled_indices) {
            if (n_valid >= n_rep) break

            # Control center = TSS + tss_distance
            # Since tss_distance is already (peak_pos - TSS), this correctly reconstructs position
            if (gene_strands[idx] == "-") {
                control_center <- gene_ends[idx] + tss_distance
            } else {
                control_center <- gene_starts[idx] + tss_distance
            }
            
            control_start <- control_center - half_width
            control_end <- control_center + half_width - 1
            control_chr <- gene_chrs[idx]
            
            if (control_start >= 1 && control_end <= chr_sizes[control_chr]) {
                control_list[[control_idx]] <- GRanges(
                    seqnames = control_chr,
                    ranges = IRanges(start = control_start, end = control_end)
                )
                control_idx <- control_idx + 1
                n_valid <- n_valid + 1
            }
        }

        if (n_valid < n_rep) {
            warning(sprintf("Only found %d/%d valid control regions for target %d", 
                           n_valid, n_rep, i))
        }
    }

    control_list <- control_list[1:(control_idx - 1)]
    if (length(control_list) == 0) {
        warning("No valid control regions found")
        return(GRanges())
    }

    control_gr <- unlist(GRangesList(control_list))
    return(control_gr)
}
