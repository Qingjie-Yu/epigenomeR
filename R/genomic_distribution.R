# Assigns feature annotations to query regions based on overlaps with a library.
# Returns a data.frame with feature counts and percentages.
annotate_GRanges <- function(queries, library, feature_col, mode, feature_order = NULL) {
  n_queries <- length(queries)

  # Find overlaps between queries and library
  overlaps <- GenomicRanges::findOverlaps(queries, library)
  if (length(overlaps) == 0) {
    result <- data.frame(feature = "other", count = n_queries, percentage = 100, stringsAsFactors = FALSE)
    return(result)
  }

  # Extract overlapping regions
  query_hits <- S4Vectors::queryHits(overlaps)
  lib_hits <- S4Vectors::subjectHits(overlaps)
  query_gr <- queries[query_hits]
  lib_gr <- library[lib_hits]
  hit_features <- as.character(GenomicRanges::mcols(lib_gr)[[feature_col]])

  if (mode == "nearest") {
    # Assign each query to closest feature by distance to center
    query_centers <- GenomicRanges::start(query_gr) +
      (GenomicRanges::end(query_gr) - GenomicRanges::start(query_gr)) %/% 2
    dist_to_start <- abs(query_centers - GenomicRanges::start(lib_gr))
    dist_to_end <- abs(query_centers - GenomicRanges::end(lib_gr))
    min_dists <- pmin(dist_to_start, dist_to_end)

    dt <- data.table::data.table(query_idx = query_hits, feature = hit_features, min_dist = min_dists)
    dt_nearest <- dt[dt[, .I[which.min(min_dist)], by = query_idx]$V1]

    features <- rep("other", n_queries)
    features[dt_nearest$query_idx] <- dt_nearest$feature
    feature_counts <- table(features)

    result <- data.frame(
      feature = names(feature_counts),
      count = as.vector(feature_counts),
      percentage = as.vector(feature_counts) / n_queries * 100,
      stringsAsFactors = FALSE
    )
  } else {
    # Weight features by overlap length proportion
    overlap_ranges <- GenomicRanges::pintersect(query_gr, lib_gr)
    overlap_widths <- GenomicRanges::width(overlap_ranges)
    query_widths <- GenomicRanges::width(query_gr)
    overlap_props <- overlap_widths / query_widths

    dt <- data.table::data.table(feature = hit_features, weight = overlap_props)
    feature_weights <- dt[, .(count = sum(weight)), by = feature]
    other_weight <- n_queries - sum(feature_weights$count)

    result <- data.frame(
      feature = c(feature_weights$feature, "other"),
      count = c(feature_weights$count, other_weight),
      percentage = c(feature_weights$count, other_weight) / n_queries * 100,
      stringsAsFactors = FALSE
    )
  }

  # Sort by feature_order if provided
  if (!is.null(feature_order)) {
    full_order <- feature_order[feature_order %in% result$feature]
    result <- result[match(full_order, result$feature), ]
  } else {
    result <- result[order(result$percentage, decreasing = TRUE), ]
  }

  rownames(result) <- NULL
  return(result)
}

# Wrapper function that processes GRanges or GRangesList inputs and returns
# feature annotation percentages in table format (samples × features).
annotate_by_overlap <- function(queries, library, feature_col, mode = c("nearest", "weighted"), feature_order = NULL) {
  mode <- match.arg(mode)

  if (!feature_col %in% colnames(GenomicRanges::mcols(library))) {
    stop("Column '", feature_col, "' not found in library metadata")
  }

  # Prepare feature order
  all_features <- unique(as.character(GenomicRanges::mcols(library)[[feature_col]]))
  if (!is.null(feature_order)) {
    missing_features <- setdiff(feature_order, all_features)
    if (length(missing_features) > 0) {
      feature_order <- feature_order[feature_order %in% all_features]
    }
    if (!"other" %in% feature_order) {
      feature_order <- c(feature_order, "other")
    }
  }

  if (inherits(queries, "GRanges")) {
    result_df <- annotate_GRanges(
      queries = queries,
      library = library,
      feature_col = feature_col,
      mode = mode,
      feature_order = feature_order
    )
    pct_vec <- setNames(result_df$percentage, result_df$feature)

    if (!is.null(feature_order)) {
      full_vec <- setNames(rep(0, length(feature_order)), feature_order)
      full_vec[names(pct_vec)] <- pct_vec
      pct_vec <- full_vec
    }

    result_table <- as.data.frame(t(pct_vec))
    rownames(result_table) <- "query"
    return(result_table)
  }

  if (!inherits(queries, "GRangesList")) {
    stop("queries must be GRanges or GRangesList")
  }

  result_matrix <- t(sapply(queries, function(query_gr) {
    result <- annotate_GRanges(
      queries = query_gr,
      library = library,
      feature_col = feature_col,
      mode = mode,
      feature_order = feature_order
    )

    if (!is.null(feature_order)) {
      pct_vec <- setNames(rep(0, length(feature_order)), feature_order)
      pct_vec[result$feature] <- result$percentage
    } else {
      pct_vec <- setNames(result$percentage, result$feature)
    }
    pct_vec
  }))

  rownames(result_matrix) <- names(queries)
  result_df <- as.data.frame(result_matrix)
  return(result_df)
}



# Repeat Element Annotation Function
# Description:
#   Annotates genomic regions with repeat elements using RepeatMasker data,
#   computes overlap-based repeat class percentages, and generates a summary
#   table along with an optional visualization.
#
# Parameters:
#   query_grl: GRangesList containing query region sets (e.g., biclustering-derived clusters).
#   out_dir: Directory to save output files (repeat tables and plots).
#   ref_genome: Reference genome version. Must be either "hg38" or "mm10".
#   mode: Annotation mode for assigning repeat classes to query regions. Default: "nearest"
#         - "nearest": Assigns each query region to the closest overlapping repeat element
#                      based on distance from query center to repeat boundaries. Each region
#                      receives exactly one repeat class annotation.
#         - "weighted": Assigns repeat classes proportionally based on overlap length. A single
#                       query region can contribute fractionally to multiple repeat classes,
#                       weighted by (overlap_length / query_length). Total weights sum to 100%.
#   plot: Boolean, whether to generate barplot.
#
#   Output files (always saved):
#       - repeat_annotation.tsv: Tab-delimited table of repeat class percentages
#         (rows = samples, columns = repeat categories).
#       - repeat_annotation.pdf: Stacked barplot visualization of repeat distributions.
repeat_distribution <- function(query_grl, out_dir = "./", ref_genome = "hg38", mode = "nearest", plot = TRUE) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
  })

  repeat_order <- c(
    "SINE", "LINE", "LTR", "Retroposon", "RC", "DNA",
    "Satellite", "Simple_repeat", "Low_complexity",
    "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "RNA",
    "Unknown"
  )

  # Load reference library
  if (ref_genome == "hg38") {
    repeat_library_file <- download_rds("RepeatMasker_hg38.rds")
  } else if (ref_genome == "mm10") {
    repeat_library_file <- download_rds("RepeatMasker_mm10.rds")
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  repeat_library <- readRDS(system.file("extdata", repeat_library_file, package = "epigenomeR"))
  message(glue("Using reference genome {ref_genome} with {length(repeat_library)} repeats"))
  style <- seqlevelsStyle(query_grl)[1]
  if (seqlevelsStyle(repeat_library)[1] != style) {
    seqlevelsStyle(repeat_library) <- style
  }

  # Annotate query regions with repeat elements
  repeat_df <- annotate_by_overlap(queries = query_grl, library = repeat_library, feature_col = "repClass", mode = mode, feature_order = repeat_order)
  write.table(repeat_df, file = file.path(out_dir, "repeat_annotation.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Convert to long format for plotting
  repeat_df_with_rownames <- tibble::rownames_to_column(repeat_df, var = "cluster")
  repeat_long <- repeat_df_with_rownames %>% pivot_longer(cols = -cluster, names_to = "Feature", values_to = "Frequency")

  # Set factor levels for proper ordering
  feature_levels <- colnames(repeat_df)
  repeat_long$Feature <- factor(repeat_long$Feature, levels = feature_levels)
  cluster_levels <- rownames(repeat_df)
  repeat_long$cluster <- factor(repeat_long$cluster, levels = rev(cluster_levels))
  
  # Create stacked barplot
  if (plot) {
    repeat_plot <- ggplot(repeat_long, aes(x = cluster, y = Frequency, fill = Feature)) +
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() + 
      scale_fill_brewer(palette = "Dark2") + 
      labs(title = "Repeat Elements", x = NULL, y = "Percentage (%)") +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_text(size = 20)
      ) +
      guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))
    
    ggsave(file.path(out_dir, "repeat_annotation.pdf"), plot = repeat_plot, width = 6, height = 5)
  }
}


# ChromHMM Annotation Function
# Description:
#   Annotates genomic regions with chromatin states using ChromHMM data,
#   computes overlap-based chromatin state percentages, and generates a summary
#   table along with an optional visualization.
#
# Parameters:
#   query_grl: GRangesList containing query region sets (e.g., biclustering-derived clusters).
#   out_dir: Directory to save output files (ChromHMM tables and plots). Default: "./"
#   ref_genome: Reference genome version. Must be either "hg38" or "mm10". Default: "hg38"
#   mode: Annotation mode for assigning chromatin states to query regions. Default: "nearest"
#         - "nearest": Assigns each query region to the closest overlapping chromatin state
#                      based on distance from query center to state boundaries. Each region
#                      receives exactly one chromatin state annotation.
#         - "weighted": Assigns chromatin states proportionally based on overlap length. A single
#                       query region can contribute fractionally to multiple chromatin states,
#                       weighted by (overlap_length / query_length). Total weights sum to 100%.
#   plot: Boolean, whether to generate barplot.
#
# Output files (always saved):
#   - chromhmm_annotation.tsv: Tab-delimited table of chromatin state percentages
#     (rows = samples, columns = chromatin state categories).
#   - chromhmm_annotation.pdf: Stacked barplot visualization of chromatin state distributions.
chromhmm_distribution <- function(query_grl, out_dir = "./", ref_genome = "hg38", mode = "nearest", plot = TRUE) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
  })

  hmm_order <- c(
    "Acet", "EnhWk", "EnhA", "PromF", "TSS",        # active enhancers & promoters
    "TxWk", "TxEx", "Tx",                           # transcription continuum
    "OpenC", "TxEnh",                               # accessible chromatin & enhancer-coupled transcription
    "ReprPCopenC",
    "BivProm", "ZNF",                               # developmental regulation & ZNF clusters
    "ReprPC", "HET",                                # Polycomb repression & constitutive heterochromatin
    "Quies"                                         # quiescent background
  )

  # Load reference library
  if (ref_genome == "hg38") {
    hmm_library_file <- download_rds("ChromHMM_hg38.rds")
  } else if (ref_genome == "mm10") {
    hmm_library_file <- download_rds("ChromHMM_mm10.rds")
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  hmm_library <- readRDS(system.file("extdata", hmm_library_file, package = "epigenomeR"))
  message(glue("Using reference genome {ref_genome} with {length(hmm_library)} ChromHMMs"))
  style <- seqlevelsStyle(query_grl)[1]
  if (seqlevelsStyle(hmm_library)[1] != style) {
    seqlevelsStyle(hmm_library) <- style
  }

  # Annotate query regions with ChromHMM elements
  hmm_df <- annotate_by_overlap(queries = query_grl, library = hmm_library, feature_col = "clean_class", mode = mode, feature_order = hmm_order)
  write.table(hmm_df, file = file.path(out_dir, "chromhmm_annotation.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Convert to long format for plotting
  hmm_df_with_rownames <- tibble::rownames_to_column(hmm_df, var = "cluster")
  hmm_long <- hmm_df_with_rownames %>% pivot_longer(cols = -cluster, names_to = "Feature", values_to = "Frequency")

  # Set factor levels for proper ordering
  feature_levels <- colnames(hmm_df)
  hmm_long$Feature <- factor(hmm_long$Feature, levels = feature_levels)
  cluster_levels <- rownames(hmm_df)
  hmm_long$cluster <- factor(hmm_long$cluster, levels = rev(cluster_levels))
  
  # Create stacked barplot
  if (plot) {
    hmm_plot <- ggplot(hmm_long, aes(x = cluster, y = Frequency, fill = Feature)) +
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() + 
      scale_fill_brewer(palette = "Paired") + 
      labs(title = "ChromHMM States", x = NULL, y = "Percentage (%)") +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_text(size = 20)
      ) +
      guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))
    
    ggsave(file.path(out_dir, "chromhmm_annotation.pdf"), plot = hmm_plot, width = 6, height = 5)
  }
}


# cCRE Element Annotation Function
# Description:
#   Annotates genomic regions with cCRE elements using ENCODE data,
#   computes overlap-based repeat class percentages, and generates a summary
#   table along with an optional visualization.
#
# Parameters:
#   query_grl: GRangesList containing query region sets (e.g., biclustering-derived clusters).
#   out_dir: Directory to save output files (repeat tables and plots).
#   ref_genome: Reference genome version. Must be either "hg38" or "mm10".
#   ref_source: Gene annotation source used to define gene models and TSS
#               coordinates. Supported options are:
#               - "knownGene": Uses the UCSC knownGene annotation obtained via
#                 the Bioconductor package TxDb.Hsapiens.UCSC.hg38.knownGene
#               - "GENCODE": Uses GENCODE gene annotations (GENCODE v49).
#   mode: Annotation mode for assigning repeat classes to query regions. Default: "nearest"
#         - "nearest": Assigns each query region to the closest overlapping repeat element
#                      based on distance from query center to repeat boundaries. Each region
#                      receives exactly one repeat class annotation.
#         - "weighted": Assigns repeat classes proportionally based on overlap length. A single
#                       query region can contribute fractionally to multiple repeat classes,
#                       weighted by (overlap_length / query_length). Total weights sum to 100%.
#   plot: Boolean, whether to generate barplot.
#
#   Output files (always saved):
#       - repeat_annotation.tsv: Tab-delimited table of repeat class percentages
#         (rows = samples, columns = repeat categories).
#       - repeat_annotation.pdf: Stacked barplot visualization of repeat distributions.
ccre_distribution <- function(query_grl, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
  })

  ccre_order <- c(
    "dELS", "pELS", "PLS",
    "5' UTR", "Exon", "Intron", "3' UTR",
    "CA-H3K4me3", "CA-CTCF", "CA-TF", "TF", "CA"
  )

  genome_config <- list(
    hg38 = list(
        ccre_library_file = download_rds("ENCODE_cCRE_v4_hg38.rds"),
        knowngene_file = download_rds("knownGene_hg38.rds"),
        gencode_file = download_rds("GENCODE_v49_hg38.rds")
    ),
    mm10 = list(
        ccre_library_file = download_rds("ENCODE_cCRE_v4_mm10.rds"),
        knowngene_file = download_rds("knownGene_mm10.rds"),
        gencode_file = download_rds("GENCODE_vM35_mm10.rds")
    )
  )
  config <- genome_config[[ref_genome]]

  ccre_library <- readRDS(system.file("extdata", config$ccre_library_file, package = "epigenomeR"))
  
  if (ref_source == "knownGene") {
    gene_library <- readRDS(system.file("extdata", config$knowngene_file, package = "epigenomeR"))
  } else if (ref_source == "GENCODE") {
    gene_library <- readRDS(system.file("extdata", config$gencode_file, package = "epigenomeR"))
  } else {
    stop("ref_source must be 'knownGene' or 'GENCODE'")
  }
  gene_library <- gene_library[gene_library$ccre_class != "Promoter"]

  style <- seqlevelsStyle(query_grl)[1]
  if (seqlevelsStyle(ccre_library)[1] != style) {
    seqlevelsStyle(ccre_library) <- style
  }
  if (seqlevelsStyle(gene_library)[1] != style) {
    seqlevelsStyle(gene_library) <- style
  }

  # Combine libraries with gene_library having priority
  combine_gr <- function(gr1, gr2) {
    overlaps <- findOverlaps(gr2, gr1, type = "any")
    if (length(overlaps) == 0) {
      return(sort(c(gr1, gr2)))
    }
    overlap_idx <- unique(queryHits(overlaps))
    non_overlap_idx <- setdiff(seq_along(gr2), overlap_idx)
    gr2_no_overlap <- gr2[non_overlap_idx]
    gr2_overlap <- gr2[overlap_idx]
    gr2_remain <- setdiff(gr2_overlap, gr1)
    combined <- c(gr1, gr2_remain, gr2_no_overlap)
    return(sort(combined))
  }
  combined_library <- combine_gr(gene_library, ccre_library)

  # Annotate query regions with combined elements
  ccre_df <- annotate_by_overlap(queries = query_grl, library = combined_library, feature_col = "ccre_class", mode = mode, feature_order = ccre_order)
  write.table(ccre_df, file = file.path(out_dir, "ccre_annotation.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Convert to long format for plotting
  ccre_df_with_rownames <- tibble::rownames_to_column(ccre_df, var = "cluster")
  ccre_long <- ccre_df_with_rownames %>% pivot_longer(cols = -cluster, names_to = "Feature", values_to = "Frequency")

  # Set factor levels for proper ordering
  feature_levels <- colnames(ccre_df)
  ccre_long$Feature <- factor(ccre_long$Feature, levels = feature_levels)
  cluster_levels <- rownames(ccre_df)
  ccre_long$cluster <- factor(ccre_long$cluster, levels = rev(cluster_levels))
  
  # Create stacked barplot
  if (plot) {
    ccre_plot <- ggplot(ccre_long, aes(x = cluster, y = Frequency, fill = Feature)) +
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() + 
      scale_fill_brewer(palette = "Set3") + 
      labs(title = "cCRE States", x = NULL, y = "Percentage (%)") +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_text(size = 20)
      ) +
      guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))
    
    ggsave(file.path(out_dir, "ccre_annotation.pdf"), plot = ccre_plot, width = 6, height = 5)
  }
}