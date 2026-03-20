peak_differential_analysis <- function(peak_dirs, bam_dirs, conditions, sample_names = NULL, out_dir = "./", pattern = "_peaks\\.narrowPeak$", window_size = NULL, min_support = 2, lfc_threshold = 0.5, fdr_threshold = 0.05, mean_quantile = 0.25) {
  # input validation
  if (length(peak_dirs) != length(conditions)) {
    stop("peak_dirs and conditions must have the same length.")
  }
  if (length(bam_dirs) != length(conditions)) {
    stop("bam_dirs and conditions must have the same length.")
  }

  cond_table <- table(conditions)
  if (any(cond_table < 2)) {
    stop("Each condition must have at least 2 replicates. Failed for: ",
         paste(names(cond_table)[cond_table < 2], collapse = ", "))
  }

  # scan directories
  scan_dir <- function(d, pat) {
    paths <- list.files(d, pattern = pat, full.names = TRUE)
    setNames(paths, stringr::str_replace(basename(paths), pat, ""))
  }
  peak_maps <- lapply(peak_dirs, scan_dir, pat = pattern)
  bam_maps  <- lapply(bam_dirs,  scan_dir, pat = "\\.bam$")

  common_peak_pairs <- Reduce(intersect, lapply(peak_maps, names))
  common_bam_pairs <- Reduce(intersect, lapply(bam_maps,  names))
  common_pairs <- intersect(common_peak_pairs, common_bam_pairs)
  if (length(common_pairs) == 0)
    stop("No common pairs found across all peak_dirs and bam_dirs.")
  message(length(common_pairs), " common pairs found.")

  # output directories
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # get cores
  BPPARAM <- get_BPPARAM()

  # main loop
  results <- bplapply(common_pairs, function(pair) {
    message("\n── Processing pair: ", pair)

    peak_path <- vapply(peak_maps, function(m) m[[pair]], character(1))
    bam_path  <- vapply(bam_maps,  function(m) m[[pair]], character(1))

    # step 1: build consensus peak set
    regions_bed <- tryCatch(
      build_peak_set(
        peak_path   = peak_path,
        conditions  = conditions,
        pair        = pair,
        out_dir     = out_dir,
        window_size = window_size,
        min_support = min_support
      ),
      error = function(e) {
        warning("build_peak_set failed for '", pair, "': ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(regions_bed)) return(NULL)

    # step 2: differential analysis
    tryCatch(
      {
        sig_paths <- differential_regions_single_peak(
          bam_path      = bam_path,
          conditions    = conditions,
          regions       = regions_bed,
          pair          = pair,
          sample_names  = sample_names,
          out_dir       = out_dir,
          min_support   = min_support,
          lfc_threshold = lfc_threshold,
          fdr_threshold = fdr_threshold,
          mean_quantile = mean_quantile
        )
        list(regions_bed = regions_bed, sig_paths = sig_paths)
      },
      error = function(e) {
        warning("differential_regions_single_peak failed for '", pair, "': ", conditionMessage(e))
        NULL
      }
    )
  }, BPPARAM = BPPARAM)
  names(results) <- common_pairs
  invisible(results)
}