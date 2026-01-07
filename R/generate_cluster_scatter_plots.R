# differential -3
# Post: Generate scatter plots comparing log2 fold changes between genomic bins and RNA-seq expression for genes within specified distance, performing Fisher's exact test to assess concordance between epigenetic and transcriptomic changes across clusters.
# Parameter: rnaseq_file_path: Path to CSV file containing RNA-seq differential expression results with gene IDs as row names
#            promoter_file_path: Path to TSV file containing promoter annotations with genomic coordinates and gene information
#            load_dir: Directory containing feather files with limma differential analysis results for each cluster
#            sig_result_dir: Output directory where scatter plots and summary tables will be saved
#            cluster_idx_list: Vector of cluster IDs to process (default: all cluster)
#            nearest_dist_cutoff: Maximum distance in bp for mapping genomic bins to promoters (default: 5000)
#            l2fc_thres: Log2 fold change threshold for defining significant changes (default: 0.5)
#            ylim_top: Upper limit for y-axis (RNA-seq log2FC) in scatter plots (default: 7)
#            ylim_bottom: Lower limit for y-axis (RNA-seq log2FC) in scatter plots (default: -4)
#            xlim_top: Upper limit for x-axis (genomic bin logFC) in scatter plots (default: 3)
#            xlim_bottom: Lower limit for x-axis (genomic bin logFC) in scatter plots (default: -3)
# Output: Saves scatter plots as PDF files and TSV summary tables showing quadrant analysis of concordant/discordant changes, with Fisher's exact test results for each cluster
generate_cluster_scatter_plots <- function(rnaseq_file_path, promoter_file_path, load_dir, sig_result_dir, cluster_idx_list = NULL, nearest_dist_cutoff = 5000, l2fc_thres = 0.5, ylim_top = 7, ylim_bottom = -4, xlim_top = 3, xlim_bottom = -3) {
  # load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(GenomicRanges)
    library(glue)
    library(ggpointdensity)
    library(viridis)
    library(ComplexHeatmap)
    library(ggplot2)
    library(ggpubr)
    library(ggrastr)
    library(data.table)
    library(dplyr)
    library(circlize)
  })

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  # load RNA-seq results
  rnaseq_raw <- read.table(rnaseq_file_path, sep = ",", header = TRUE, check.names = FALSE, row.names = 1)
  rnaseq_raw <- rnaseq_raw[!is.na(rnaseq_raw$log2FoldChange), ]
  rnaseq_raw$log2FoldChange_TPM <- rnaseq_raw$log2_T - rnaseq_raw$log2_V

  # load promoter annotations
  promoter_table <- read.table(promoter_file_path, sep = "\t", header = TRUE)
  promoter_table$seqnames <- paste0("chr", promoter_table$seqnames)
  promoter_grange <- makeGRangesFromDataFrame(
    promoter_table,
    seqnames.field = "seqnames", start.field = "start", end.field = "end",
    strand.field = "strand", keep.extra.columns = TRUE
  )

  if (is.null(cluster_idx_list)) {
    cluster_idx_list <- sort(unique(col_cluster_full$label))
    message(glue("cluster_idx_list not specified. Using all unique labels from column cluster file: {paste(cluster_idx_list, collapse=', ')}"))
  }

  for (cluster_id in cluster_idx_list) {
    message(glue("Processing cluster {cluster_id}..."))
    file_name <- glue("result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-0.25_column_cluster-{cluster_id}_limma.feather")
    file_path <- file.path(load_dir, file_name)
    pc <- read_feather(file_path)

    colnames(pc)[1] <- "pos"
    pc <- pc[, c("pos", "adj.P.Val", "logFC")]
    pc <- cbind(pc, do.call(rbind, strsplit(pc$pos, "[_]")))
    df_wgc <- pc[, c("logFC", "adj.P.Val", "1", "2", "3", "pos")]
    colnames(df_wgc)[ncol(df_wgc)] <- "genomic_bin_pos"
    df_wgc$mid <- round((as.numeric(df_wgc[, "3"]) - as.numeric(df_wgc[, "2"])) %/% 2) + as.numeric(df_wgc[, "2"])

    gr_wgc <- makeGRangesFromDataFrame(df_wgc, seqnames.field = "1", start.field = "mid", end.field = "mid", keep.extra.columns = TRUE)

    # map to nearest promoters
    overlap_res <- distanceToNearest(gr_wgc, promoter_grange)
    overlap_dist <- mcols(overlap_res)$distance
    overlap_res <- overlap_res[overlap_dist <= nearest_dist_cutoff]
    overlap_df <- cbind(
      mcols(gr_wgc[queryHits(overlap_res), ]),
      mcols(promoter_grange[subjectHits(overlap_res), ])
    )
    overlap_df <- cbind(
      overlap_df,
      rnaseq_raw[
        overlap_df$gene_id,
        c(
          "log2FoldChange_TPM", "padj",
          "NaB.1.RNA.seq_S3", "NaB.2.RNA.seq_S4",
          "veh.1.RNA.seq_S1", "veh.2.RNA.seq_S2"
        )
      ]
    )
    overlap_df <- overlap_df[which(overlap_df$padj <= 0.05), ]
    overlap_df <-
      overlap_df[, c(
        "logFC",
        "adj.P.Val",
        "genomic_bin_pos",
        "pos",
        "gene_id",
        "gene_name",
        "log2FoldChange_TPM",
        "padj",
        "NaB.1.RNA.seq_S3",
        "NaB.2.RNA.seq_S4",
        "veh.1.RNA.seq_S1",
        "veh.2.RNA.seq_S2"
      )]
    df_to_draw <- overlap_df
    df_to_draw <- df_to_draw[order(df_to_draw$logFC), ]
    df_summarized <- as.data.frame(df_to_draw) %>%
      group_by(gene_id) %>%
      summarise(
        gene_name = first(gene_name),
        pos = first(pos),
        genomic_bin_pos = paste(unique(genomic_bin_pos), collapse = ","),
        adj.P.Val = min(adj.P.Val, na.rm = TRUE),
        logFC = mean(logFC, na.rm = TRUE),
        padj = first(padj),
        log2FoldChange_TPM = first(log2FoldChange_TPM),
        .groups = "drop"
      )

    # quadrants for overlapping significant bins
    df_to_draw <- df_summarized %>%
      mutate(
        quadrant = case_when(
          logFC >= l2fc_thres & log2FoldChange_TPM >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 1",
          logFC < -l2fc_thres & log2FoldChange_TPM >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 2",
          logFC < -l2fc_thres & log2FoldChange_TPM < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 3",
          logFC >= l2fc_thres & log2FoldChange_TPM < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 4"
        )
      )

    df_to_draw_clean <- df_to_draw[!is.na(df_to_draw$quadrant), ]
    write.table(df_to_draw_clean,
      file.path(sig_result_dir, glue("result_merge_bins_column_cluster-{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    # fisher test for Q1+Q4 vs Q2+Q3
    quadrant_counts <- table(df_to_draw_clean$quadrant)
    test_tb <- matrix(
      c(
        quadrant_counts["Quadrant 2"], quadrant_counts["Quadrant 1"],
        quadrant_counts["Quadrant 3"], quadrant_counts["Quadrant 4"]
      ),
      nrow = 2, byrow = TRUE
    )
    test_tb[is.na(test_tb)] <- 0
    test_tb <- test_tb + 1
    fisher_test_result <- fisher.test(test_tb)
    annotation_text <- glue("Cluster {cluster_id}: OR={round(1/fisher_test_result$estimate,2)}, P={signif(fisher_test_result$p.value,3)}")

    # scatter plot
    df_to_draw$quadrant[is.na(df_to_draw$quadrant)] <- "NONE"
    p <- ggplot(df_to_draw, aes(x = logFC, y = log2FoldChange_TPM)) +
      geom_point(size = 1.5, color = "#e9ecef") +
      geom_pointdensity(data = df_to_draw[df_to_draw$quadrant != "NONE", ], aes(logFC, log2FoldChange_TPM), size = 1.5) +
      scale_color_viridis() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#adb5bd") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#adb5bd") +
      geom_smooth(data = df_to_draw[df_to_draw$quadrant != "NONE", ], aes(logFC, log2FoldChange_TPM), method = "lm", color = "red", se = TRUE) +
      coord_cartesian(xlim = c(xlim_bottom, xlim_top), ylim = c(ylim_bottom, ylim_top)) +
      ggtitle(annotation_text) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
      )

    scatter_plot_filename <- glue("scatter_{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.pdf")
    ggsave(file.path(sig_result_dir, scatter_plot_filename), plot = p, width = 8 / 1.5, height = 8 / 1.5, bg = "transparent")
  }

  message("Finished all clusters.")
}
