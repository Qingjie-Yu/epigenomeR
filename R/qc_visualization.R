#' Generate CRF-by-CRF Matrix from Pairwise Data
#'
#' Creates a symmetric matrix from CRF pairs, filling only the lower triangle.
#' Pairs that don't pass QC (not in filtered_df) are set to NA.
#'
#' @param all_df Data frame with columns `file` and `read_count`
#' @param filtered_df Data frame with columns `file` and `read_count`
#' @param group_csv Optional: file path to CSV with CRF grouping info
#' @param crf_col Column name in group_csv for CRF identifiers (default: "crf")
#' @param category_col Column name in group_csv for categories (default: "category")
#' @param by Delimiter used in pair names (default: "-")
#'
#' @return List with:
#'   - mat: Numeric matrix (lower triangle filled, upper = 0)
#'   - category: Category labels for each CRF (NULL if no group_csv)
create_read_count_matrix <- function(all_df, filtered_df, group_csv = NULL, crf_col = "crf", category_col = "category", by = "-") {
  # Extract CRFs
  if (!is.null(group_csv)) {
    if (!file.exists(group_csv)) {
      stop(glue::glue("File not found: {group_csv}"))
    }
    group_df <- read.csv(group_csv, stringsAsFactors = FALSE)
    missing_cols <- setdiff(c(crf_col, category_col), colnames(group_df))
    if (length(missing_cols) > 0) {
      stop(glue::glue("Missing required column(s): {paste(missing_cols, collapse = ', ')}"))
    }

    # Sort by category
    group_df <- group_df[order(group_df[[category_col]]), ]
    crfs <- group_df[[crf_col]]
    category_vector <- group_df[[category_col]]
  } else {
    # Extract unique CRFs from all pairs
    pairs <- all_df$pair
    valid_pairs <- pairs[grepl(by, pairs, fixed = TRUE)]

    if (length(valid_pairs) == 0) {
      stop(glue::glue("No valid pairs found with delimiter '{by}'"))
    }
    crfs <- unique(unlist(strsplit(valid_pairs, by, fixed = TRUE)))
    crfs <- crfs[grepl("^[A-Za-z]", crfs)] # keep only CRFs starting with letter

    if (length(crfs) == 0) {
      stop("No valid CRF names extracted from pairs")
    }
    category_vector <- NULL
  }

  # Create and fill matrix
  n_crfs <- length(crfs)
  mat <- matrix(0,
    nrow = n_crfs, ncol = n_crfs,
    dimnames = list(crfs, crfs)
  )

  pair_split <- strsplit(as.character(all_df$pair), by, fixed = TRUE)
  crf1_vec <- sapply(pair_split, `[`, 1)
  crf2_vec <- sapply(pair_split, `[`, 2)
  values <- all_df$read_count

  crf_to_idx <- setNames(seq_along(crfs), crfs)
  valid_mask <- (crf1_vec %in% crfs) & (crf2_vec %in% crfs)
  crf1_vec <- crf1_vec[valid_mask]
  crf2_vec <- crf2_vec[valid_mask]
  values <- values[valid_mask]

  idx1 <- crf_to_idx[crf1_vec]
  idx2 <- crf_to_idx[crf2_vec]

  row_idx <- pmax(idx1, idx2)
  col_idx <- pmin(idx1, idx2)
  mat[cbind(row_idx, col_idx)] <- values

  # Mark non-QC-passed pairs as NA
  pass_pairs <- unique(as.character(filtered_df$pair))
  pass_set <- new.env(hash = TRUE, parent = emptyenv())
  for (pair in pass_pairs) {
    pass_set[[pair]] <- TRUE
    pair_split <- strsplit(pair, by, fixed = TRUE)[[1]]
    reversed <- paste(rev(pair_split), collapse = by)
    pass_set[[reversed]] <- TRUE
  }

  lower_tri_idx <- which(lower.tri(mat, diag = FALSE), arr.ind = TRUE)
  pair_names <- paste(crfs[lower_tri_idx[, 1]],
    crfs[lower_tri_idx[, 2]],
    sep = by
  )
  is_fail <- !sapply(pair_names, function(p) exists(p, envir = pass_set))
  mat[lower_tri_idx[is_fail, , drop = FALSE]] <- NA

  return(list(mat = mat, category = category_vector))
}

#' Plot CRF Pair Read Count Heatmap
#'
#' Generates a heatmap (PDF) showing log2-transformed read counts for CRF pairs.
#' Only the lower triangle is filled; filtered pairs are shown in grey.
#'
#' @param mat Numeric matrix from create_read_count_matrix()
#' @param category Category vector for each CRF (NULL if no grouping)
#' @param out_dir Output directory for PDF (default: "./")
#' @param filename Output PDF filename (default: "qc_heatmap.pdf")
#' @param width PDF width in inches (default: 12)
#' @param height PDF height in inches (default: 12)
plot_read_count_heatmap <- function(mat, category = NULL, out_dir = "./", pdf_name = "qc_heatmap.pdf", width = 12, height = 12) {
  # Load libraries
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(latex2exp)
    library(circlize)
    library(grid)
  })

  # Validate input
  if (!is.matrix(mat)) {
    stop("mat must be a matrix")
  }
  if (nrow(mat) != ncol(mat)) {
    stop("mat must be a square matrix")
  }

  # Check and create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  output_pdf <- file.path(out_dir, pdf_name)

  # Log2 transformation
  log_mat <- log2(mat + 1)
  log_min <- min(log_mat, na.rm = TRUE)
  log_max <- max(log_mat, na.rm = TRUE)
  log_avg <- (log_min + log_max) / 2

  # Create NA legend
  na_legend <- Legend(
    labels = "Filtered Pairs",
    legend_gp = gpar(fill = "grey90"),
    title = ""
  )

  ht_opt$ROW_ANNO_PADDING <- unit(0.4, "cm")
  ht_opt$COLUMN_ANNO_PADDING <- unit(0.4, "cm")


  if (!is.null(category)) {
    # With category grouping
    unique_categories <- unique(category)
    n_categories <- length(unique_categories)

    # 12-color palette
    palette_12 <- c(
      "#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2",
      "#BEB8DC", "#E7DAD2", "#999999", "#BC3C29",
      "#0072B5", "#E18727", "#20854E", "#7876B1"
    )
    group_colors <- rep_len(palette_12, n_categories)

    # Create annotations
    row_annotation <- rowAnnotation(
      category = anno_block(
        gp = gpar(fill = group_colors),
        labels = unique_categories,
        labels_gp = gpar(col = "white", fontsize = 11)
      )
    )

    col_annotation <- HeatmapAnnotation(
      category = anno_block(
        gp = gpar(fill = group_colors),
        labels = unique_categories,
        labels_gp = gpar(col = "white", fontsize = 11)
      )
    )

    # Create heatmap with grouping
    pdf(output_pdf, width = width, height = height)
    ht <- Heatmap(
      log_mat,
      rect_gp = gpar(type = "none"),
      na_col = "grey90",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      # Only draw lower triangle
      cell_fun = function(j, i, x, y, w, h, fill) {
        if (i >= j) {
          grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }
      },
      # Labels
      row_labels = TeX(rownames(mat)),
      column_labels = TeX(colnames(mat)),
      row_names_side = "left",
      column_title = "CRF Pair Read Count",
      # Color scheme
      col = colorRamp2(
        c(log_min, log_avg, log_max),
        c("blue", "white", "red")
      ),
      # Legend
      heatmap_legend_param = list(
        title = "log2(read count + 1)",
        color_bar = "continuous",
        title_gp = gpar(fontsize = 15)
      ),
      # Splits by category
      column_split = category,
      row_split = category,
      row_title = rep("", n_categories),

      # Annotations
      left_annotation = row_annotation,
      bottom_annotation = col_annotation
    )

    draw(ht, annotation_legend_list = list(na_legend), merge_legend = TRUE)
    dev.off()
  } else {
    # Without category grouping
    pdf(output_pdf, width = width, height = height)
    ht <- Heatmap(
      log_mat,
      rect_gp = gpar(type = "none"),
      na_col = "grey90",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      # Only draw lower triangle
      cell_fun = function(j, i, x, y, w, h, fill) {
        if (i >= j) {
          grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }
      },
      # Labels
      row_labels = TeX(rownames(mat)),
      column_labels = TeX(colnames(mat)),
      row_names_side = "left",
      column_title = "CRF Pair Read Count",
      # Color scheme
      col = colorRamp2(
        c(log_min, log_avg, log_max),
        c("blue", "white", "red")
      ),
      # Legend
      heatmap_legend_param = list(
        title = "log2(read count + 1)",
        color_bar = "continuous",
        title_gp = gpar(fontsize = 15)
      )
    )

    draw(ht, annotation_legend_list = list(na_legend), merge_legend = TRUE)
    dev.off()
  }

  cat("Heatmap saved to:", output_pdf, "\n")
}


#' Visualize QC results generated by qc_count.sh
#'
#' It converts the pairwise table into a symmetric CRF-by-CRF matrix (lower triangle filled),
#' sets non-passing pairs to NA (shown as grey in the heatmap), applies log2(read_count + 1)
#' transformation, and saves a PDF heatmap under out_dir.
#'
#' Requirements for input TSVs:
#' - Must be tab-delimited with header
#' - Must contain columns: `pair` and `read_count`
#' - Pair names must be splittable by `split_pair_by` (default "-")
#'
#' @param all_read_count_path Path to all_qc.tsv generated by qc_count.sh
#' @param filtered_read_count_path Path to filtered_qc.tsv generated by qc_count.sh
#' @param out_dir Output directory for the heatmap PDF (default: "./")
#' @param split_pair_by Delimiter used in pair names to split CRFs (default: "-")
#' @param group_csv Optional CSV path defining CRF grouping (columns: crf_col, category_col)
#' @param crf_col Column name in group_csv for CRF identifiers (default: "crf")
#' @param category_col Column name in group_csv for CRF categories (default: "category")
#'
#' @return Invisibly returns a list containing `mat` and `category` used for plotting.
qc_visualization <- function(all_read_count_path, filtered_read_count_path, out_dir = "./", split_pair_by = "-", group_csv = NULL, crf_col = "crf", category_col = "category") {

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  all_df <- read.table(all_read_count_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  filtered_df <- read.table(filtered_read_count_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # sanity check
  req_cols <- c("pair", "read_count")
  if (length(setdiff(req_cols, colnames(all_df))) > 0) stop("all_df must contain: pair, read_count")
  if (!("pair" %in% colnames(filtered_df))) stop("filtered_df must contain: pair")

  cat("The following pairs pass the quality control:\n")
  cat(paste(filtered_df$pair, collapse = ", "), "\n")

  result <- create_read_count_matrix(all_df = all_df, filtered_df = filtered_df, group_csv = group_csv, crf_col = crf_col, category_col = category_col, by = split_pair_by
  )
  plot_read_count_heatmap(mat = result$mat, category = result$category, out_dir = out_dir)
}