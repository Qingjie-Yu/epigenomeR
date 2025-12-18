draw_heatmap <- function(data, out_path, col_fun, name, apply_cluster = FALSE, aspect_ratio = 2) {
  if (apply_cluster) {
    h_temp <- Heatmap(data, name = name, col = col_fun, cluster_columns = TRUE, show_heatmap_legend = FALSE)
    ht_drawn <- draw(h_temp)
    col_order <- column_order(ht_drawn)
    ordered_colnames <- rev(colnames(data)[col_order])
    plot_data <- data[, ordered_colnames, drop = FALSE]
  } else {
    plot_data <- data
  }
  
  cell_width <- 8 
  col_fontsize <- 10
  cell_height <- cell_width / aspect_ratio
  row_fontsize <- col_fontsize / aspect_ratio
  legend_height <- nrow(plot_data) * unit(cell_height, "mm") / 3

  h <- Heatmap(
    plot_data, 
    name = name, 
    col = col_fun,
    cluster_columns = FALSE, 
    cluster_rows = FALSE,
    row_names_side = "left", 
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_side = "bottom", 
    column_names_gp = gpar(fontsize = col_fontsize),
    column_names_rot = 45,
    width = ncol(plot_data) * unit(cell_width, "mm"), 
    height = nrow(plot_data) * unit(cell_height, "mm"), 
    show_heatmap_legend = FALSE
  )

  lgd <- Legend(
    col_fun = col_fun,
    title = name,
    direction = "vertical",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 6),
    legend_height = legend_height
  )

  pdf_width <- convertWidth(
    unit(ncol(plot_data) * cell_width, "mm") + unit(50, "mm"), 
    "inches", 
    valueOnly = TRUE
  )
  pdf_height <- convertHeight(
    unit(nrow(plot_data) * cell_height, "mm") + unit(20, "mm"), 
    "inches", 
    valueOnly = TRUE
  )
  
  pdf(out_path, width = pdf_width, height = pdf_height)
  draw(h, heatmap_legend_side = "right", padding = unit(c(4, 4, 4, 4), "mm"))
  draw(lgd,  x = unit(1, "npc") - unit(30, "mm"), y = unit(0.5, "npc"), just = c("left", "center"))
  dev.off()  
  cat(glue("Saved heatmap: {out_path}"), "\n")
}

# Peak Enrichment Heatmap
# Post: Generate heatmaps of TFBS enrichment (log2 odds ratio) from enrichment result tsv files.
# Parameter: tsv_path: vector of paths to TFBS enrichment result tsv files
#            label: vector of sample/cluster labels corresponding to tsv_path
#            out_dir: output directory
#            top_n: Number of top enriched TFBS (based on coefficient of variation) to display in second heatmap
#            selected_tfs: Vector of transcription factor names or substrings to filter rows
#            apply_cluster: Whether to cluster columns and reorder them
# Output: Heatmap PDF (all filtered TFBS): TFBS_enrichment_all.pdf
#         Heatmap PDF (top n TFBS by CV): TFBS_heatmap_top<n>.pdf (if top_n specified)
#         log2 odds ratio matrix: TFBS_odds_ratio_log2.csv
#         FDR matrix: TFBS_FDR.csv

TFBS_enrichment_heatmap <- function(tsv_path, label, out_dir, top_n = NULL, selected_tfs = NULL, apply_cluster = FALSE) {
  # load library
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(grid)
    library(glue)
    library(circlize)
  })

  if (length(tsv_path) != length(label)) {
    stop("tsv_path and label must have the same length")
  }
  
  if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
  }

  # reorder
  sorted_indices <- tryCatch({
    order(as.numeric(label))
  }, warning = function(w) {
    order(label)
  }, error = function(e) {
    order(label)
  })
  label <- label[sorted_indices]
  tsv_path <- tsv_path[sorted_indices]

  # read first tsv file
  first_file <- read.table(tsv_path[[1]], sep = "\t", header = TRUE, row.names = 1)
  tfbs_ids <- rownames(first_file)
  n_tfbs <- length(tfbs_ids)
  n_samples <- length(label)

  # build matrices
  odds_ratio_mat <- matrix(NA, nrow = n_tfbs, ncol = n_samples, dimnames = list(tfbs_ids, label))
  FDR_mat <- matrix(NA, nrow = n_tfbs, ncol = n_samples, dimnames = list(tfbs_ids, label))
  target_hit_mat <- matrix(NA, nrow = n_tfbs, ncol = n_samples, dimnames = list(tfbs_ids, label))

  for (i in seq_along(tsv_path)) {
    file_path <- tsv_path[[i]]
    sample_label <- label[i]
    if (!file.exists(file_path)) {
      warning(glue("File {file_path} does not exist. Skipping."))
      next
    }
    df <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1)

    if (!all(rownames(df) == tfbs_ids)) {
      warning(glue("TFBS IDs mismatch in {sample_label}, attempting to match by row names"))
      common_ids <- intersect(tfbs_ids, rownames(df))
      odds_ratio_mat[common_ids, i] <- df[common_ids, "odds_ratio"]
      FDR_mat[common_ids, i] <- df[common_ids, "FDR"]
      target_hit_mat[common_ids, i] <- df[common_ids, "target_hit"]
    } else {
      odds_ratio_mat[, i] <- df$odds_ratio
      FDR_mat[, i] <- df$FDR
      target_hit_mat[, i] <- df$target_hit
    }
  }

  # remove row/columns with all NAs
  valid_rows <- rowSums(!is.na(odds_ratio_mat)) > 0
  valid_cols <- colSums(!is.na(odds_ratio_mat)) > 0
  odds_ratio_mat <- odds_ratio_mat[valid_rows, valid_cols, drop = FALSE]
  FDR_mat <- FDR_mat[valid_rows, valid_cols, drop = FALSE]
  target_hit_mat <- target_hit_mat[valid_rows, valid_cols, drop = FALSE]

  # filter
  target_hit_thresh <- quantile(target_hit_mat, 0.1, na.rm = TRUE)
  keep_tfbs <- rowSums(odds_ratio_mat >= 2, na.rm = TRUE) >= 1 &
    rowSums(FDR_mat <= 0.05, na.rm = TRUE) >= 1 &
    rowSums(target_hit_mat >= target_hit_thresh, na.rm = TRUE) >= 1 
  odds_ratio_mat <- odds_ratio_mat[keep_tfbs, , drop = FALSE]
  FDR_mat <- FDR_mat[keep_tfbs, , drop = FALSE]
  cat(glue("Filtered to {sum(keep_tfbs)} TFBS from {length(keep_tfbs)} total"), "\n")

  # apply selected_tfs filter if specified
  if (!is.null(selected_tfs)) {
    selected_pattern <- paste(selected_tfs, collapse = "|")
    selected_rows <- grepl(selected_pattern, rownames(odds_ratio_mat), ignore.case = TRUE)
    if (sum(selected_rows) == 0) {
      warning("No TFBS matched the selected_tfs patterns.")
    } else {
      odds_ratio_mat <- odds_ratio_mat[selected_rows, , drop = FALSE]
      FDR_mat <- FDR_mat[selected_rows, , drop = FALSE]
      cat(glue("Filtered to {sum(selected_rows)} TFBS after applying selected_tfs filter"), "\n")
    }
  }

  # stop check
  if (nrow(odds_ratio_mat) == 0) {
    stop("No TFBS left after filtering. Cannot generate heatmap.")
  }
  cat(glue("Final matrix: {nrow(odds_ratio_mat)} TFBS × {ncol(odds_ratio_mat)} samples"), "\n")

  # first heatmap (all)
  odds_ratio_log2 <- log2(odds_ratio_mat)
  max_abs_val <- max(abs(odds_ratio_log2), na.rm = TRUE)
  col_fun <- colorRamp2(c(-max_abs_val, 0, max_abs_val), c("#3155C3", "white", "#AF0525"))
  draw_heatmap(data = odds_ratio_log2, out_path = file.path(out_dir, "TFBS_enrichment_all.pdf"), col_fun = col_fun, name = "log2(Odds Ratio)", apply_cluster = apply_cluster) 

  # second heatmap (top n)
  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    cv <- apply(odds_ratio_mat, 1, sd, na.rm = TRUE) / abs(rowMeans(odds_ratio_mat, na.rm = TRUE))
    top_tfbs <- names(head(sort(cv, decreasing =  TRUE), top_n))
    odds_ratio_log2_topn <- odds_ratio_log2[top_tfbs, , drop = FALSE]
    cat(glue("Selected top {top_n} TFBS based on coefficient of variation"), "\n")
    draw_heatmap(data = odds_ratio_log2_topn, out_path = file.path(out_dir, paste0("TFBS_heatmap_top", top_n, ".pdf")), col_fun = col_fun, name = "log2(Odds Ratio)", apply_cluster = apply_cluster, aspect_ratio = 1)
  }

  # save .csv
  write.csv(odds_ratio_log2, file.path(out_dir, "TFBS_odds_ratio_log2.csv"))
  write.csv(FDR_mat, file.path(out_dir, "TFBS_FDR.csv"))
  cat(glue("Saved odds ratio matrix: {file.path(out_dir, 'TFBS_odds_ratio_log2.csv')}"), "\n")
  cat(glue("Saved FDR matrix: {file.path(out_dir, 'TFBS_FDR.csv')}"), "\n")
}
