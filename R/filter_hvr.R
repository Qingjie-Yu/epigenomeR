# Filter highly variable regions based on mean-variance analysis
# Description: Filters genomic regions to identify the most informative features
#              based on log2-transformed mean expression and hypervariance values
# Parameters:
#   transformed_cm_path: Path to transformed count matrix file (.feather format with "pos" column)
#   mav_stats_path: Path to mean-variance statistics file (.feather) from MAV screening
#   out_dir: Output directory path (default: "./")
#   n_bins: Number of bins to divide log2(mean) range for stratified sampling (default: 100)
#   keep_percent: Fraction of total regions to keep (0-1), distributed across bins (default: 0.01)
#   log2mean_quantile_thres: Percentile threshold for log2(mean) to retain highly expressed regions (default: 0.99)
#   plot: Whether to generate diagnostic scatter plots (default: FALSE)
# Output:
#   - Filtered count matrix (.feather): Contains only selected informative regions
#   - Optional diagnostic plots (.png): Shows selection process across expression levels

filter_hvr <- function(transformed_cm_path, mav_stats_path, out_dir = "./", n_bins = 100, keep_percent = 0.01, log2mean_quantile_thres = 0.99, plot = FALSE, seed = 42) {
  # Load libraries
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(cowplot)
    library(arrow)
  })

  # Create folder
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Read mean-variance statistics (only needed columns)
  mav_stats <- arrow::read_feather(mav_stats_path,
    col_select = c("pos", "log2_mean", "log2_var", "hypervar")
  )

  # Calculate mean expression threshold
  mean_threshold <- quantile(mav_stats$log2_mean, log2mean_quantile_thres)
  message(
    "Log2(mean) threshold (", log2mean_quantile_thres * 100,
    "th percentile): ", round(mean_threshold, 3)
  )

  # Pre-filter: only keep overdispersed regions (hypervar > 1)
  mav_filtered <- mav_stats %>% filter(hypervar > 1)

  # Stratified sampling: bin regions by log2(mean)
  mav_filtered$bin <- cut(mav_filtered$log2_mean, breaks = n_bins, labels = FALSE, include.lowest = TRUE)

  # Calculate regions to keep per bin
  n_per_bin <- floor(nrow(mav_stats) * keep_percent / n_bins)
  message(
    "Selecting ~", n_per_bin, " regions per bin (",
    keep_percent * 100, "% of total)"
  )

  # Select top hypervariant regions from each bin
  set.seed(42)
  selected_regions <- mav_filtered %>%
    dplyr::group_by(bin) %>%
    dplyr::slice_max(order_by = hypervar, n = n_per_bin, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-bin)

  # Final filter: keep only highly expressed regions
  final_regions <- selected_regions %>%
    dplyr::filter(log2_mean >= mean_threshold)

  message(
    "Final selected regions: ", nrow(final_regions), " (",
    round(nrow(final_regions) / nrow(mav_stats) * 100, 2), "% of total)"
  )

  # Load count matrix and extract selected regions
  count_matrix <- arrow::read_feather(transformed_cm_path)
  selected_matrix <- count_matrix %>%
    dplyr::filter(pos %in% final_regions$pos)

  # Save filtered matrix
  input_name <- tools::file_path_sans_ext(basename(transformed_cm_path))
  feather_path <- file.path(out_dir, paste0(input_name, "_filtered_regions.feather"))
  arrow::write_feather(selected_matrix, feather_path)
  message("Saved filtered matrix to: ", feather_path)

  # Generate diagnostic plots if requested
  if (plot == TRUE) {
    plot_path <- file.path(out_dir, paste0(input_name, "_filtered_regions.png"))

    # Plot 1: Mean-variance relationship
    p1 <- ggplot(mav_stats, aes(x = log2_mean, y = log2_var)) +
      geom_point(alpha = 0.3, color = "#D4E9D9", size = 1) +
      geom_point(
        data = final_regions,
        aes(x = log2_mean, y = log2_var),
        size = 2, color = "#F6E0D3", alpha = 0.8
      ) +
      labs(x = "Log2(mean)", y = "Log2(variance)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    # Plot 2: Hypervariance selection process
    p2 <- ggplot() +
      geom_point(
        data = mav_stats,
        aes(x = log2_mean, y = hypervar),
        alpha = 0.3, color = "blue", size = 0.5
      ) +
      geom_point(
        data = selected_regions,
        aes(x = log2_mean, y = hypervar),
        alpha = 0.3, color = "red", size = 1
      ) +
      geom_point(
        data = final_regions,
        aes(x = log2_mean, y = hypervar),
        size = 3, color = "green", alpha = 0.8
      ) +
      labs(x = "Log2(mean)", y = "Hypervariance") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    # Combine and save
    combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
    ggsave(plot_path, plot = combined_plot, width = 8, height = 5, dpi = 300)
    message("Saved plot to: ", plot_path)
  }

  return(feather_path)
}
