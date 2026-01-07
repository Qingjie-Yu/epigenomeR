# differential -4
# Post: Summarize the number of significant differentially accessible regions by reading pre-computed differential analysis results
# Parameter: load_dir: Directory containing differential analysis results from limma_column_cluster_differential_regions
#            sig_result_dir: Output directory where summary tables will be saved
#            cluster_idx_list: Vector of column cluster labels to summarize
#            l2fc_thres: Log2 fold change threshold for significance calling (default: 0.5)
#            mean_per_thres_list: Vector of mean expression percentile thresholds (default: c(0.25))
#            fdr_thres_list: Vector of FDR thresholds for significance testing (default: c(0.25))
# Output: Saves TSV summary tables showing counts of significant, upregulated, and downregulated regions for each column cluster
summary_sig_num <- function(load_dir, sig_result_dir, cluster_idx_list, l2fc_thres = 0.5, mean_per_thres_list = c(0.25), fdr_thres_list = c(0.25)) {
  suppressPackageStartupMessages({
    library(arrow)
    library(glue)
    library(dplyr)
  })

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  for (fdr_thres in fdr_thres_list) {
    for (mean_per_thres in mean_per_thres_list) {
      fdr_sig_summary_table <- data.frame(
        column_cluster = cluster_idx_list,
        sig = 0,
        up = 0,
        down = 0
      )

      for (i in seq_along(cluster_idx_list)) {
        col_label <- cluster_idx_list[i]

        # Read pre-computed results from limma_column_cluster_differential_regions
        result_file <- glue("{load_dir}/result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_limma.feather")

        if (file.exists(result_file)) {
          top.table <- read_feather(result_file)

          # Calculate statistics based on thresholds
          sig_regions <- top.table %>% filter(adj.P.Val < fdr_thres, abs(logFC) > l2fc_thres)
          up_regions <- top.table %>% filter(adj.P.Val < fdr_thres, logFC > 0)
          down_regions <- top.table %>% filter(adj.P.Val < fdr_thres, logFC <= 0)

          fdr_sig_summary_table[i, "sig"] <- nrow(sig_regions)
          fdr_sig_summary_table[i, "up"] <- nrow(up_regions)
          fdr_sig_summary_table[i, "down"] <- nrow(down_regions)

          message(glue("sig thres: {fdr_thres}, filter: {mean_per_thres}, column cluster: {col_label}, sig regions: {nrow(sig_regions)}"))
        } else {
          warning(glue("File not found: {result_file}"))
        }
      }

      out_file <- glue("summary_sig_num_{mean_per_thres}_FDR-{fdr_thres}_log2FC-{l2fc_thres}.tsv")
      write.table(fdr_sig_summary_table, file.path(sig_result_dir, out_file), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}
