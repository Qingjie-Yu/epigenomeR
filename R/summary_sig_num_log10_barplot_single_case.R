# differential -5
# Post: Generate log10-scaled barplot visualization of significant region counts showing upregulated and downregulated regions across column clusters with original count labels and pseudo-log transformation for better visualization of wide dynamic ranges.
# Parameter: result_file_dir: Directory containing TSV file with summary statistics (columns 'column_cluster', 'up', 'down' from summary_sig_num function)
#            sig_result_dir: Output directory where the barplot PDF will be saved
#            mean_per_thres: Mean expression percentile threshold used in summary file name (default: "0.25")
#            fdr_thres: FDR threshold used in summary file name (default: "0.25")
#            l2fc_thres: Log2 fold change threshold used in summary file name (default: 0.5)
# Output: Saves PDF barplot with log10-transformed y-axis showing upregulated regions (red, positive) and downregulated regions (blue, negative) for each column cluster, with original count values as text labels and pseudo-log scale for enhanced readability
summary_sig_num_log10_barplot_single_case <- function(result_file_dir, sig_result_dir, mean_per_thres = "0.25", fdr_thres = "0.25", l2fc_thres = 0.5) {
  # load libraries
  suppressPackageStartupMessages({
    library(glue)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  result_file_path <- glue("{result_file_dir}/summary_sig_num_{mean_per_thres}_FDR-{fdr_thres}_log2FC-{l2fc_thres}.tsv")

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  result_table <- read.table(result_file_path, header = TRUE, sep = "\t")

  # ---- SAVE ORIGINAL COUNTS BEFORE LOG TRANSFORM ----
  result_table$orig_up <- result_table$up
  result_table$orig_down <- result_table$down

  # ---- LOG10 TRANSFORM WITH PSEUDOCOUNT 1 ----
  result_table$up <- log10(result_table$up + 1)
  result_table$down <- log10(result_table$down + 1)

  df_long <- result_table %>%
    select(column_cluster, up, down, orig_up, orig_down) %>%
    pivot_longer(
      cols = c(up, down),
      names_to = "status",
      values_to = "count"
    ) %>%
    mutate(
      # map to plotting length
      count_plot = if_else(status == "down", -count, count),
      # carry original count for labels
      original_count = if_else(status == "down", orig_down, orig_up),
      status = factor(status, levels = c("up", "down"))
    )

  y_top <- max(df_long$count_plot)
  y_bot <- min(df_long$count_plot)

  result_figure <- ggplot(
    df_long,
    aes(
      x = factor(column_cluster, levels = unique(column_cluster)),
      y = count_plot,
      fill = status
    )
  ) +
    # bars
    geom_col(width = 0.7) +

    # text labels: use ORIGINAL counts
    geom_text(
      data = filter(df_long, status == "up"),
      aes(label = original_count),
      vjust = -0.3, size = 4
    ) +
    geom_text(
      data = filter(df_long, status == "down"),
      aes(label = original_count),
      vjust = 1.3, size = 4
    ) +

    # manual colors
    scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +

    # zero line
    geom_hline(yintercept = 0, colour = "black") +

    # log10 scale on y-axis with manual breaks for pseudo_log
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(
        -log10(1000), -log10(100), -log10(10),
        log10(10), log10(100), log10(1000)
      ),
      labels = c("1000", "100", "10", "10", "100", "1000"),
      sec.axis = sec_axis(
        transform = ~.,
        breaks = c(-1.5, 1.5),
        labels = c("Down Regulated", "Up Regulated")
      )
    ) +

    # no legend, no y axis title
    labs(x = "Column Clusters", y = NULL) +
    guides(fill = "none") +

    # minimal theme, axes only
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(
        angle  = 0,
        size   = 14,
        colour = "black",
        hjust  = 0.5,
        vjust  = 0.5
      ),
      axis.text.y.right = element_text(
        angle  = 90,
        size   = 14,
        colour = "black",
        hjust  = 0.5,
        vjust  = 0.5
      ),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  result_figure_filename <- paste0(basename(tools::file_path_sans_ext(result_file_path)), "_log.pdf")
  result_figure_dir_filename <- file.path(sig_result_dir, result_figure_filename)
  ggsave(result_figure_dir_filename, result_figure, bg = "transparent", width = 8, height = 5)
}
