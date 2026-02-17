# MAV Screen Function
# Purpose: Model mean-variance relationship and identify overdispersed features.
#
# Parameters:
#   transformed_cm_path: Single input .feather file containing transformed count matrix (first column must be "pos").
#   out_dir: Output directory (default = "./")
#   fitting_model: Model to fit mean-variance trend - "gam" (default) or "loess".
#   k: For "gam" - number of basis functions (default = 40 if NULL); ignored for "loess".
#   span: For "loess" - smoothing span (default = 0.5 if NULL); ignored for "gam".
#   seed: Random seed for reproducibility in plot sampling (default = 42).
#   font_size: Font size for plot labels and axes (default = 10).
#   nrow_sample_per: Proportion of rows to randomly sample for density plot (default = 0.2, i.e., 20%).
#   plot: Whether to generate and save diagnostic plots (default = FALSE).
#
# Output:
#   Annotated dataframe saved as a .feather file and (optional) two diagnostic plots saved as .PNG.
#   Returns the path of the .feather file.

mav_screen <- function(transformed_cm_path, out_dir = "./", fitting_model = "gam", k = NULL, span = NULL, seed = 42, font_size = 10, nrow_sample_per = 0.2, plot = FALSE) {
  # Load library
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(matrixStats)
    library(mgcv)
    library(ggplot2)
    library(cowplot)
    library(ggpointdensity)
    library(tools)
  })

  # Create folder
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Load data
  if (!file.exists(transformed_cm_path)) {
    stop("Input file does not exist: ", transformed_cm_path)
  }

  df <- arrow::read_feather(transformed_cm_path)
  pos <- df$pos
  df$pos <- NULL
  mat <- as.matrix(df)
  mode(mat) <- "numeric"

  region_mean <- rowMeans(mat)
  region_var <- rowVars(mat)
  data_fit <- data.frame(mean = region_mean, var = region_var)

  # Set regression model
  message("Fitting ", fitting_model, " model...")
  if (fitting_model == "gam") {
    if (is.null(k)) k <- 40
    fit_model <- gam(log2(var) ~ s(log2(mean), k = k), data = data_fit)
  } else {
    if (is.null(span)) span <- 0.5
    fit_model <- loess(log2(var) ~ log2(mean), data = data_fit, span = span)
  }

  # Initial fitting
  var_expected <- 2^(fit_model$fitted)
  sd_expected <- sqrt(var_expected)
  mat_normalized <- (mat - region_mean) / sd_expected

  norm_mean <- rowMeans(mat_normalized)
  norm_var <- rowVars(mat_normalized)
  hypervar <- rowSums(mat_normalized^2) / (ncol(mat) - 1)

  # Save filtered count matrix
  result <- data.frame(
    pos = as.character(pos),
    mean = region_mean,
    var = region_var,
    var_expect = var_expected,
    norm_mean = norm_mean,
    norm_var = norm_var,
    hypervar = hypervar,
    log2_mean = log2(region_mean),
    log2_var = log2(region_var),
    log2_var_expect = log2(var_expected),
    stringsAsFactors = FALSE
  )

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  input_name <- basename(tools::file_path_sans_ext(transformed_cm_path))
  feather_path <- file.path(out_dir, paste0(input_name, "_mav_screen.feather"))
  write_feather(result, feather_path)
  message("Saved results to: ", feather_path)

  # Plot
  if (plot == TRUE) {
    set.seed(seed)
    message("Generating diagnostic plots...")
    main_plot_path <- file.path(out_dir, paste0(input_name, "_mean_variance.png"))
    density_plot_path <- file.path(out_dir, paste0(input_name, "_fit_density.png"))

    # Plot A: Mean-variance relationship with fitted curve
    p1 <- ggplot(result, aes(log2_mean, log2_var)) +
      geom_point(alpha = 1 / 20) +
      geom_point(aes(y = log2_var_expect), color = "red", size = 0.1) +
      labs(x = "log2(mean)", y = "log2(variance)") +
      theme_bw() +
      theme(
        axis.text = element_text(size = font_size),
        axis.title = element_text(size = font_size)
      )

    # Plot B: Hypervariance distribution across mean expression
    p2 <- ggplot(result, aes(log2_mean, hypervar)) +
      geom_point(alpha = 1 / 20) +
      labs(x = "log2(mean)", y = "Hypervariance") +
      theme_bw() +
      theme(
        axis.text = element_text(size = font_size),
        axis.title = element_text(size = font_size)
      )

    # Combine and save main diagnostic plots
    combined_plot <- plot_grid(p1, p2, labels = c("A", "B"))
    ggsave(main_plot_path, plot = combined_plot, width = 8, height = 5)
    message("  Saved main plot to: ", main_plot_path)

    # Generate density visualization with sampled data points
    nrow_sample <- min(nrow(result), ceiling(nrow(result) * nrow_sample_per))
    result_sample <- result[sample(nrow(result), nrow_sample), ]
    message("  Sampling ", nrow_sample, " regions for density plot")

    p3 <- ggplot(result_sample, aes(log2_mean, log2_var)) +
      geom_pointdensity() +
      scale_color_viridis_c() +
      geom_point(data = result, aes(y = log2_var_expect), color = "red", size = 0.1) +
      labs(x = "log2(mean)", y = "log2(variance)") +
      theme_bw() +
      theme(
        axis.text = element_text(size = font_size),
        axis.title = element_text(size = font_size)
      )

    ggsave(density_plot_path, plot = p3, width = 4, height = 5)
    message("  Saved density plot to: ", density_plot_path)
  }

  return(feather_path)
}
