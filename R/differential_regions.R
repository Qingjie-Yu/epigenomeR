differential_regions <- function(cm_path, conditions, sample_names = NULL, col_cluster_file_path = NULL, out_dir = "./", lfc_threshold = 0.5, fdr_threshold = 0.05, mean_quantile = 0.25, pseudocount = 1) {
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(edgeR)
    library(matrixStats)
    library(limma)
  })

  if (length(cm_path) != length(conditions)) {
    stop("cm_path and conditions must have the same length.")
  }

  cond_table <- table(conditions)
  if (any(cond_table < 2)) {
    stop(paste0(
      "Each condition must have at least 2 replicates. Failed for: ",
      paste(names(cond_table)[cond_table < 2], collapse = ", ")
    ))
  }

  ord <- order(conditions)
  cm_path    <- cm_path[ord]
  conditions <- conditions[ord]

  if (is.null(sample_names)) {
    sample_names <- paste0(conditions, ave(seq_along(conditions), conditions, FUN = seq_along))
  } else {
    if (length(sample_names) != length(conditions))
      stop("Length of sample_names must equal length of conditions.")
    if (anyDuplicated(sample_names))
      stop("sample_names must be unique.")
    sample_names <- sample_names[ord]
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Read feather files
  cm_list <- lapply(cm_path, function(f) {
    df <- arrow::read_feather(f)
    if (!("pos" %in% colnames(df))) stop("Missing required column 'pos' in: ", f)
    pos    <- df$pos
    df$pos <- NULL
    mat    <- as.matrix(df)
    mode(mat) <- "numeric"
    rownames(mat) <- as.character(pos)
    mat
  })
  names(cm_list) <- sample_names

  # Harmonize row names
  all_regions <- lapply(cm_list, rownames)
  if (length(unique(sapply(all_regions, paste, collapse = "|"))) > 1) {
    warning("Row names (pos) do not match across all count matrices. Taking intersection.")
    common_regions <- Reduce(intersect, all_regions)
    if (length(common_regions) == 0) stop("No common row names (pos) found.")
    cm_list <- lapply(cm_list, function(mat) mat[common_regions, , drop = FALSE])
  }

  # Harmonize col names
  all_pairs <- lapply(cm_list, colnames)
  if (length(unique(sapply(all_pairs, paste, collapse = "|"))) > 1) {
    warning("Column names (pairs) do not match across all count matrices. Taking intersection.")
    common_pairs <- Reduce(intersect, all_pairs)
    if (length(common_pairs) == 0) stop("No common column names found.")
    cm_list <- lapply(cm_list, function(mat) mat[, common_pairs, drop = FALSE])
  }

  cm_regions <- rownames(cm_list[[1]])
  cm_pairs   <- colnames(cm_list[[1]])

  # Column cluster table
  if (is.null(col_cluster_file_path)) {
    col_cluster_table <- data.frame(
      pair    = cm_pairs,
      cluster = seq_along(cm_pairs),
      stringsAsFactors = FALSE
    )
  } else {
    col_cluster_table <- read.table(col_cluster_file_path, header = TRUE,
                                     sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("pair", "cluster") %in% colnames(col_cluster_table)))
      stop("col_cluster_file must contain columns: pair and cluster")
    col_cluster_table <- col_cluster_table[col_cluster_table$pair %in% cm_pairs, , drop = FALSE]
    if (nrow(col_cluster_table) == 0)
      stop("After filtering by cm_pairs, col_cluster_table is empty.")
    col_cluster_table$cluster <- as.integer(col_cluster_table$cluster)
    if (any(is.na(col_cluster_table$cluster)))
      stop("col_cluster_file: 'cluster' column contains non-integer / NA values.")
  }

  cm_not_in_cluster <- setdiff(cm_pairs, col_cluster_table$pair)
  if (length(cm_not_in_cluster) > 0)
    warning("Count matrices contain ", length(cm_not_in_cluster),
            " column(s) not found in col_cluster_file (they will be ignored).")

  # Design matrix
  conditions  <- factor(conditions)
  design      <- model.matrix(~ 0 + conditions)
  colnames(design) <- levels(conditions)

  # All pairwise comparisons
  cond_levels <- levels(conditions)
  n_cond      <- length(cond_levels)
  comparisons <- list()
  idx <- 1
  for (i in 1:(n_cond - 1)) {
    for (j in (i + 1):n_cond) {
      comparisons[[idx]] <- c(cond_levels[j], cond_levels[i])
      idx <- idx + 1
    }
  }

  # Main loop
  cluster_list  <- sort(unique(col_cluster_table$cluster))

  for (cl in cluster_list) {

    target_pairs <- intersect(col_cluster_table$pair[col_cluster_table$cluster == cl], cm_pairs)
    if (length(target_pairs) == 0) next

    # Aggregate counts across pairs in this cluster
    combined <- do.call(cbind,
      lapply(sample_names, function(s)
        rowSums(cm_list[[s]][, target_pairs, drop = FALSE])
      )
    )
    colnames(combined) <- sample_names
    rownames(combined) <- cm_regions

    # Pre-filter: at least one condition with >=2 non-zero AND >50% non-zero replicates
    keep_by_condition <- sapply(cond_levels, function(cond) {
      cols   <- sample_names[conditions == cond]
      nz_cnt <- rowSums(combined[, cols, drop = FALSE] != 0)
      (nz_cnt >= 2) & (nz_cnt > length(cols) * 0.5)
    })
    if (is.vector(keep_by_condition))
      keep_by_condition <- matrix(keep_by_condition, ncol = 1)
    keep_rows  <- rowSums(keep_by_condition) > 0
    combined_f0 <- combined[keep_rows, , drop = FALSE]
    if (nrow(combined_f0) < 2) {
      warning(glue("Cluster {cl}: too few rows after non-zero filter; skipped.")) 
      next
    }

    # First voom for mean-expression filter
    d0 <- DGEList(combined_f0)
    d0 <- calcNormFactors(d0)
    y0 <- voom(d0, design, plot = FALSE)
    mean_log2_cpm <- rowMeans(y0$E)

    th          <- as.numeric(stats::quantile(mean_log2_cpm, probs = mean_quantile, na.rm = TRUE))
    combined_f1 <- combined_f0[mean_log2_cpm >= th, , drop = FALSE]
    if (nrow(combined_f1) < 2) {
      warning(glue("Cluster {cl}: too few rows after mean_quantile={mean_quantile}; skipped."))
      next
    }

    # Second voom on filtered matrix
    d1  <- DGEList(combined_f1)
    d1  <- calcNormFactors(d1)
    y1  <- voom(d1, design, plot = FALSE)
    fit <- lmFit(y1, design)

    # Pairwise comparisons
    for (cmp in comparisons) {
      test <- cmp[1]
      ref  <- cmp[2]

      contr <- makeContrasts(contr = paste0("`", test, "`-`", ref, "`"),
                              levels = colnames(coef(fit)))
      fit2  <- contrasts.fit(fit, contr)
      fit2  <- eBayes(fit2)

      tt <- topTable(fit2, sort.by = "P", n = Inf)
      tt$P.Value[tt$P.Value   == 0] <- 1e-300
      tt$adj.P.Val[tt$adj.P.Val == 0] <- 1e-300

      sig     <- tt[(tt$adj.P.Val < fdr_threshold) & (abs(tt$logFC) > lfc_threshold), , drop = FALSE]
      cmp_tag <- glue("{test}_vs_{ref}")
      prefix  <- glue("{cmp_tag}_cluster{cl}")
      # Save full table
      tt_out     <- rownames_to_column(tt, var = "pos")
      feather_all <- file.path(out_dir, glue("{prefix}_limma_all.feather"))
      write_feather(tt_out, feather_all)

      # Save sig table
      sig_tsv <- file.path(out_dir, glue("{prefix}_limma_sig.tsv"))
      if (nrow(sig) > 0) {
        write.table(sig, sig_tsv, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
      } else {
        writeLines("No significant regions under current thresholds.", sig_tsv)
      }

      # Save per-sample log2 matrix for sig regions
      if (nrow(sig) > 0) {
        sig_regions <- rownames(sig)
        for (s in sample_names) {
          m <- cm_list[[s]][sig_regions, target_pairs, drop = FALSE]
          m <- log2(m + pseudocount)
          m <- rownames_to_column(as.data.frame(m), var = "pos")
          write_feather(m, file.path(out_dir, glue("{prefix}_{s}_log2.feather")))
        }
      }

      # Summary
      summary_tsv <- file.path(out_dir, glue("{cmp_tag}_limma_summary.tsv"))
      file_exists <- file.exists(summary_tsv)
      summary_row <- data.frame(
        cluster              = cl,
        n_rows_before        = nrow(combined),
        n_rows_after_nonzero = nrow(combined_f0),
        n_rows_after_mean    = nrow(combined_f1),
        n_sig                = nrow(sig),
        n_up                 = nrow(sig[sig$logFC >  0, ]),
        n_down               = nrow(sig[sig$logFC <= 0, ]),
        stringsAsFactors     = FALSE
      )
      write.table(summary_row, summary_tsv,
        sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = !file_exists, append = file_exists
      )
      message(glue("Cluster {cl} | {cmp_tag} | tested: {nrow(combined_f1)} | sig: {nrow(sig)} (up: {nrow(sig[sig$logFC > 0,])}, down: {nrow(sig[sig$logFC <= 0,])})"))
    }
  }
}