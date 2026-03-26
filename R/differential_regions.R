# ── Internal: two-round voom + limma ────────────────────────────────────────
# combined   : region × sample integer matrix, colnames = sample_names
# conditions : character vector, same order as colnames(combined)
# returns    : named list, one entry per comparison + meta (n_before/n_nonzero/n_tested)
run_limma_voom <- function(combined, conditions, min_support = 2, mean_quantile = 0, lfc_threshold = 0.5, fdr_threshold = 0.05) {
  conditions_f <- factor(conditions)
  cond_levels  <- levels(conditions_f)
  design       <- model.matrix(~ 0 + conditions_f)
  colnames(design) <- cond_levels

  comparisons <- combn(cond_levels, 2,
                       function(x) c(x[2], x[1]), simplify = FALSE)

  # pre-filter: at least one condition with >= min_support non-zero replicates
  keep_by_condition <- sapply(cond_levels, function(cond) {
    cols   <- colnames(combined)[conditions == cond]
    nz_cnt <- rowSums(combined[, cols, drop = FALSE] != 0)
    nz_cnt >= min_support
  })
  if (is.vector(keep_by_condition))
    keep_by_condition <- matrix(keep_by_condition, ncol = 1)
  combined_f0 <- combined[rowSums(keep_by_condition) > 0, , drop = FALSE]
  if (nrow(combined_f0) < 2)
    return(list(skipped = "too few rows after non-zero filter"))

  # first voom → mean expression filter
  d0 <- edgeR::calcNormFactors(edgeR::DGEList(combined_f0))
  y0 <- limma::voom(d0, design, plot = FALSE)
  th <- quantile(rowMeans(y0$E), probs = mean_quantile, na.rm = TRUE)
  combined_f1 <- combined_f0[rowMeans(y0$E) >= th, , drop = FALSE]
  if (nrow(combined_f1) < 2)
    return(list(skipped = paste0("too few rows after mean_quantile=", mean_quantile)))

  # second voom → fit
  d1  <- edgeR::calcNormFactors(edgeR::DGEList(combined_f1))
  y1  <- limma::voom(d1, design, plot = FALSE)
  fit <- limma::lmFit(y1, design)

  results <- list(
    n_before  = nrow(combined),
    n_nonzero = nrow(combined_f0),
    n_tested  = nrow(combined_f1)
  )

  for (cmp in comparisons) {
    test <- cmp[1]; ref <- cmp[2]
    contr <- limma::makeContrasts(
      contrasts = paste0("`", test, "`-`", ref, "`"),
      levels    = colnames(coef(fit))
    )
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr))
    tt   <- limma::topTable(fit2, sort.by = "P", n = Inf)
    tt$P.Value[tt$P.Value       == 0] <- 1e-300
    tt$adj.P.Val[tt$adj.P.Val   == 0] <- 1e-300

    cmp_tag <- paste0(test, "_vs_", ref)
    results[[cmp_tag]] <- list(
      tt  = tt,
      sig = tt[tt$adj.P.Val < fdr_threshold & abs(tt$logFC) > lfc_threshold, ,
               drop = FALSE]
    )
  }
  results
}


# ── Internal: write DAR results to disk ─────────────────────────────────────
# limma_result : output of run_limma_voom()
# combined     : original region × sample matrix (for log2 feather output)
# group_name   : pair name or cluster name (used in filenames and summary)
# out_dir      : passed through
write_da_results <- function(limma_result, combined, group_name, cmp_dir_map) {
  if (!is.null(limma_result$skipped)) {
    warning("Group ", group_name, " skipped: ", limma_result$skipped)
    return(invisible(NULL))
  }

  meta_keys <- c("n_before", "n_nonzero", "n_tested", "skipped")
  cmp_tags  <- names(cmp_dir_map)
  sig_paths <- list()

  for (cmp_tag in cmp_tags) {
    cmp_dir <- cmp_dir_map[[cmp_tag]]
    prefix <- file.path(cmp_dir, group_name)

    res <- limma_result[[cmp_tag]]
    if (is.null(res)) next

    tt <- res$tt
    sig <- res$sig
    all_path <- paste0(prefix, "_all.tsv")
    write.table(cbind(pos = rownames(tt), tt), all_path,
                sep = "\t", quote = FALSE, row.names = FALSE)

    # sig stats
    if (nrow(sig) > 0) {
      sig_path <- paste0(prefix, "_sig.tsv")
      write.table(cbind(pos = rownames(sig), sig), sig_path,
                  sep = "\t", quote = FALSE, row.names = FALSE)
      sig_paths[[cmp_tag]] <- sig_path
    }

    summary_tsv <- file.path(cmp_dir, "summary.tsv")
    summary_row <- data.frame(
      group                = group_name,
      n_rows_before        = limma_result$n_before,
      n_rows_after_nonzero = limma_result$n_nonzero,
      n_rows_after_mean    = limma_result$n_tested,
      n_sig                = nrow(sig),
      n_up                 = sum(sig$logFC >  0),
      n_down               = sum(sig$logFC < 0),
      stringsAsFactors     = FALSE
    )
    write.table(summary_row, summary_tsv,
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(summary_tsv),
                append     = file.exists(summary_tsv))

    message(sprintf("%-40s | tested: %d | sig: %d (up: %d, down: %d)",
                    glue::glue("{cmp_tag} [{group_name}]"),
                    limma_result$n_tested,
                    nrow(sig), sum(sig$logFC > 0), sum(sig$logFC <= 0)))
  }
  invisible(sig_paths)
}

# ── Internal: generate summary plot ─────────────────────────────────────
differential_summary_plot <- function(tsv, pdf) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
  })

  df <- data.table::fread(tsv)
  if (!all(c("group", "n_up", "n_down") %in% colnames(df))) {
    warning(glue::glue("Summary file missing required columns (group/n_up/n_down): {tsv}"))
    return(invisible(NULL))
  }
  df$group  <- factor(df$group, levels = unique(df$group))
  df$n_down <- -df$n_down

  plot_df <- rbind(
    data.frame(group = df$group, n = df$n_up,   direction = "Up-regulated"),
    data.frame(group = df$group, n = df$n_down, direction = "Down-regulated")
  )

  y_max   <- max(abs(c(df$n_up, df$n_down)), na.rm = TRUE)
  y_limit <- ceiling(y_max * 1.15)

  p <- ggplot(plot_df, aes(x = group, y = n, fill = direction)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    scale_fill_manual(
      values = c("Up-regulated" = "#E24B4A", "Down-regulated" = "#378ADD"),
      name   = NULL
    ) +
    scale_y_continuous(
      limits = c(-y_limit, y_limit),
      labels = function(x) abs(x)
    ) +
    labs(
      x = NULL,
      y = "Number of DA regions"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      legend.position    = "top",
      legend.key.size    = unit(0.4, "cm"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )

  ggsave(pdf, p,
         width  = max(4, length(levels(df$group)) * 0.6 + 2),
         height = 4,
         dpi    = 150)
  message(glue::glue("Saved summary plot: {pdf}"))
}

# ── Public: conventional DA (one feather per sample, columns = pairs) ───────
differential_regions <- function(cm_path, conditions, sample_names = NULL, out_dir = "./", col_cluster_file_path = NULL, min_support = 2, mean_quantile = 0, lfc_threshold = 0.5, fdr_threshold = 0.05) {
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(edgeR)
    library(limma)
  })
  # input validation
  if (length(cm_path) != length(conditions))
    stop("cm_path and conditions must have the same length.")
  cond_table <- table(conditions)
  if (any(cond_table < 2))
    stop("Each condition must have at least 2 replicates. Failed for: ",
         paste(names(cond_table)[cond_table < 2], collapse = ", "))

  ord        <- order(conditions)
  cm_path    <- cm_path[ord]
  conditions <- conditions[ord]

  # create output dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  cond_levels <- unique(conditions)
  cmp_tags <- combn(
    cond_levels, 2,
    function(x) paste0(x[2], "_vs_", x[1]),
    simplify = TRUE
  ) 
  cmp_dir_map <- setNames(file.path(out_dir, cmp_tags), cmp_tags)
  invisible(lapply(cmp_dir_map, dir.create, recursive = TRUE, showWarnings = FALSE))

  # build sample names
  if (is.null(sample_names)) {
    sample_names <- paste0(conditions,
                           ave(seq_along(conditions), conditions, FUN = seq_along))
  } else {
    if (length(sample_names) != length(conditions))
      stop("Length of sample_names must equal length of conditions.")
    if (anyDuplicated(sample_names))
      stop("sample_names must be unique.")
    sample_names <- sample_names[ord]
  }

  # read feathers
  cm_list <- lapply(cm_path, function(f) {
    df <- arrow::read_feather(f)
    if (!("pos" %in% colnames(df))) stop("Missing 'pos' column in: ", f)
    pos <- df$pos; df$pos <- NULL
    mat <- as.matrix(df); mode(mat) <- "numeric"
    rownames(mat) <- as.character(pos)
    mat
  })
  names(cm_list) <- sample_names

  # harmonize regions and pairs
  harmonize <- function(cm_list, dim, label) {
    get_names <- if (dim == "row") rownames else colnames
    common    <- Reduce(intersect, lapply(cm_list, get_names))
    if (!length(common)) stop("No common ", label, " found.")
    same <- length(unique(vapply(cm_list, function(x)
      paste(sort(get_names(x)), collapse = "\r"), character(1)))) == 1
    if (!same) warning(label, " differ across samples. Using ", length(common), " common ", label, ".")
    message("Common ", label, ": ", length(common))
    if (dim == "row") lapply(cm_list, function(m) m[common, , drop = FALSE])
    else              lapply(cm_list, function(m) m[, common, drop = FALSE])
  }
  cm_list  <- harmonize(cm_list, "row", "regions")
  cm_list  <- harmonize(cm_list, "col", "pairs")
  cm_pairs <- colnames(cm_list[[1]])

  # build group list
  if (is.null(col_cluster_file_path)) {
    group_list  <- as.list(setNames(cm_pairs, cm_pairs))
    use_cluster <- FALSE
  } else {
    ct <- read.table(col_cluster_file_path, header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("pair", "cluster") %in% colnames(ct)))
      stop("col_cluster_file must contain columns: pair and cluster")
    ct <- ct[ct$pair %in% cm_pairs, , drop = FALSE]
    if (nrow(ct) == 0) stop("After filtering by cm_pairs, col_cluster_table is empty.")
    not_found <- setdiff(cm_pairs, ct$pair)
    if (length(not_found) > 0)
      warning(length(not_found), " pair(s) not in col_cluster_file and will be ignored.")
    group_list  <- split(ct$pair, factor(ct$cluster, levels = unique(ct$cluster)))
    names(group_list) <- paste0("cluster", names(group_list))
    use_cluster <- TRUE
  }

  # main loop
  result <- list()
  for (grp_name in names(group_list)) {
    target_pairs <- intersect(group_list[[grp_name]], cm_pairs)
    if (!length(target_pairs)) next

    combined <- do.call(cbind, lapply(sample_names, function(s) {
      if (use_cluster) rowSums(cm_list[[s]][, target_pairs, drop = FALSE])
      else cm_list[[s]][, target_pairs, drop = FALSE]
    }))
    colnames(combined) <- sample_names

    limma_result <- run_limma_voom(
      combined, conditions,
      min_support = min_support,
      mean_quantile = mean_quantile,
      lfc_threshold = lfc_threshold,
      fdr_threshold = fdr_threshold
    )

    grp_result <- write_da_results(limma_result, combined, grp_name, cmp_dir_map)

    # write pair-level counts for sig regions
    if (!is.null(grp_result)) {
      for (cmp_tag in names(grp_result)) {
        sig_regions <- rownames(limma_result[[cmp_tag]]$sig)
        if (!length(sig_regions)) next

        if (use_cluster) {
          # expand back to pair level: columns = sample:pair
          counts_out <- do.call(cbind, lapply(sample_names, function(s) {
            m <- cm_list[[s]][sig_regions, target_pairs, drop = FALSE]
            colnames(m) <- paste0(s, ":", colnames(m))
            m
          }))
        } else {
          # single pair: columns = sample names
          counts_out <- combined[sig_regions, , drop = FALSE]
        }

        counts_path <- file.path(cmp_dir_map[[cmp_tag]], glue("{grp_name}_sig_counts.tsv"))
        write.table(
          cbind(pos = rownames(counts_out), as.data.frame(counts_out)),
          counts_path, sep = "\t", quote = FALSE, row.names = FALSE
        )
        result[[cmp_tag]][[grp_name]] <- counts_path
      }
    }
  }

  # summary plots
  cond_levels <- sort(unique(conditions))
  for (cmp in cmp_tags) {
    summary_tsv <- file.path(cmp_dir_map[[cmp_tag]], "summary.tsv")
    summary_pdf <- file.path(cmp_dir_map[[cmp_tag]], "summary.pdf")
    if (file.exists(summary_tsv))
      differential_summary_plot(summary_tsv, summary_pdf)
  }

  invisible(result)
}


# ── Public: peak-based DA (one feather per sample per pair) ─────────────────
differential_regions_single_peak <- function(bam_path, conditions, regions, pair, sample_names = NULL, out_dir = "./", min_support = 2, lfc_threshold = 0.5, fdr_threshold = 0.05, mean_quantile = 0) {
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(Rsamtools)
    library(matrixStats)
    library(edgeR)
    library(limma)
  })

  # input validation 
  if (length(bam_path) != length(conditions))
    stop("bam_path and conditions must have the same length.")

  cond_table <- table(conditions)
  if (any(cond_table < 2))
    stop("Each condition must have at least 2 replicates. Failed for: ",
         paste(names(cond_table)[cond_table < 2], collapse = ", "))

  if (!file.exists(regions))
    stop("regions file not found: ", regions)

  ord        <- order(conditions)
  bam_path   <- bam_path[ord]
  conditions <- conditions[ord]

  # create output dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  cond_levels <- unique(conditions)
  cmp_tags <- combn(
    cond_levels, 2,
    function(x) paste0(x[2], "_vs_", x[1]),
    simplify = TRUE
  ) 
  cmp_dir_map <- setNames(file.path(out_dir, cmp_tags), cmp_tags)
  invisible(lapply(cmp_dir_map, dir.create, recursive = TRUE, showWarnings = FALSE))

  # build sample names
  if (is.null(sample_names)) {
    sample_names <- paste0(conditions,
                           ave(seq_along(conditions), conditions, FUN = seq_along))
  } else {
    if (length(sample_names) != length(conditions))
      stop("Length of sample_names must equal length of conditions.")
    if (anyDuplicated(sample_names))
      stop("sample_names must be unique.")
    sample_names <- sample_names[ord]
  }

  # Detect chromosome naming style from first BAM 
  bam_header  <- Rsamtools::scanBamHeader(bam_path[1])
  bam_seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames   = names(bam_header[[1]]$targets),
    seqlengths = bam_header[[1]]$targets
  )
  bam_style <- GenomeInfoDb::seqlevelsStyle(bam_seqinfo)[1]

  # Read BED regions 
  bed_data <- read.table(regions, sep = "\t", stringsAsFactors = FALSE)
  bin <- GenomicRanges::GRanges(
    seqnames = bed_data$V1,
    ranges   = IRanges::IRanges(start = bed_data$V2 + 1L, end = bed_data$V3)
  )
  GenomeInfoDb::seqlevelsStyle(bin) <- bam_style
  region_ids <- paste0(bed_data$V1, "_", bed_data$V2, "_", bed_data$V3)

  # Count fragments per BAM 
  combined <- do.call(cbind, lapply(seq_along(bam_path), function(k) {
    bam  <- bam_path[k]
    temp <- GenomicAlignments::readGAlignmentPairs(
      bam,
      param = Rsamtools::ScanBamParam(which = bin)
    )

    overlapCount <- numeric(length(bin))
    if (length(temp) == 0) return(overlapCount)

    locus <- data.frame(
      first_start = start(temp@first),
      first_end   = end(temp@first),
      last_start  = start(temp@last),
      last_end    = end(temp@last)
    )
    frag_start <- matrixStats::rowMins(as.matrix(locus))
    frag_end   <- matrixStats::rowMaxs(as.matrix(locus))
    bamContent <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
      seqnames = as.vector(seqnames(temp)),
      strand   = "*",
      start    = frag_start,
      end      = frag_end
    ))

    overlaps <- GenomicRanges::findOverlaps(bamContent, bin, ignore.strand = TRUE)
    qh <- queryHits(overlaps)
    sh <- subjectHits(overlaps)
    if (!length(sh)) return(overlapCount)

    overlap_lengths <- pmax(0,
      pmin(end(bamContent)[qh],   end(bin)[sh]) -
      pmax(start(bamContent)[qh], start(bin)[sh]) + 1
    )
    proportions <- overlap_lengths / (end(bamContent)[qh] - start(bamContent)[qh] + 1)
    summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
    overlapCount[summed$bin_id] <- summed$x
    overlapCount
  }))
  rownames(combined) <- region_ids
  colnames(combined) <- sample_names

  # limma 
  limma_result <- run_limma_voom(
    combined,
    conditions = conditions,
    min_support = min_support,
    mean_quantile = mean_quantile,
    lfc_threshold = lfc_threshold,
    fdr_threshold = fdr_threshold
  )
  sig_paths <- write_da_results(limma_result, combined, pair, cmp_dir_map)

  invisible(sig_paths)
}