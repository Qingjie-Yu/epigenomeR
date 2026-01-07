# Filter highly variable regions
# Description:
#   This is a convenience wrapper function that combines mean-variance screening
#   and highly variable region (HVR) selection into a single step. It automatically:
#   1. Performs MAV (Mean-Variance) screening to calculate hypervariance metrics
#   2. Filters regions based on hypervariance and expression levels
#   3. Returns a filtered count matrix containing only informative regions
#
# Parameters:
#   transformed_cm_path: Path to transformed count matrix file (.feather format with "pos" column)
#   out_dir: Output directory path (default: "./")
#   keep_percent:
#       Fraction of total regions to keep (0–1), distributed equally across expression
#       bins. This parameter controls the FIRST stage of filtering through stratified
#       binning. The algorithm divides all regions into 100 equal-width bins based on
#       their log2(mean) expression, then ranks regions within each bin by hypervariance
#       from highest to lowest. From each bin, the top N regions are selected, where
#       N = (total_regions × keep_percent) / 100.
#   log2mean_quantile_thres:
#       Percentile threshold (0–1) for log2(mean) expression that controls the SECOND
#       stage of filtering through an expression cutoff. After the first stage of
#       stratified binning, this parameter applies a final quality filter. The algorithm
#       calculates a percentile-based threshold from the ORIGINAL (unfiltered) log2(mean)
#       distribution of all regions, then retains only those regions from Stage 1 that
#       have log2(mean) greater than or equal to this threshold.
#       plot: Whether to generate diagnostic scatter plots (default: FALSE)
#
# Output Files:
#   1. <input>_mav_screen.feather - Mean-variance statistics (intermediate)
#   2. <input>_filtered_regions.feather - Final filtered count matrix (main output)
#   3. <input>_mav_screen.png - MAV diagnostic plot (if plot=TRUE)
#   4. <input>_filtered_regions.png - HVR selection plot (if plot=TRUE)
# Returns:
#   Character string: File path to the filtered count matrix (.feather format)

detect_hvr <- function(transformed_cm_path, out_dir = "./", keep_percent = 0.01, log2mean_quantile_thres = 0.99, plot = FALSE) {
  mav_path <- mav_screen(transformed_cm_path = transformed_cm_path, out_dir = out_dir, plot = plot)
  filtered_cm_path <- filter_hvr(transformed_cm_path = transformed_cm_path, mav_stats_path = mav_path, out_dir = out_dir, keep_percent = keep_percent, log2mean_quantile_thres = log2mean_quantile_thres, plot = plot)
  return(filtered_cm_path)
}
