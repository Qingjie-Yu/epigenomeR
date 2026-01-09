# Post: Identify two valley positions in fragment length distribution using kernel density estimation for nucleosome positioning analysis
# Parameter: frags: Numeric vector of fragment lengths from sequencing data
#            dens_reso: Resolution for kernel density estimation (default: 2^15 = 32768 points for smooth curve)
#            dens_kernel: Kernel type for density estimation (default: "gaussian")
#            valley1_range: Search range for first valley separating NFR from mononucleosome (default: c(73, 221) bp)
#            valley2_range: Search range for second valley separating mononucleosome from dinucleosome (default: c(222, 368) bp)
# Output: Named numeric vector with valley1 and valley2 positions in base pairs (rounded to integers); returns NA if valleys not found
find_two_valleys <- function(frags, dens_reso = 2^15, dens_kernel = "gaussian", valley1_range = c(73, 221), valley2_range = c(222, 368)) {
  # Input validation
  if (!is.numeric(frags) || length(frags) == 0) {
    stop("Error: 'frags' must be a non-empty numeric vector")
  }
  if (dens_reso <= 0 || !is.numeric(dens_reso)) {
    stop("Error: 'dens_reso' must be a positive number")
  }

  # Compute kernel density estimation to create smooth distribution curve
  density_result <- density(frags, kernel = dens_kernel, n = dens_reso)
  hist_x <- density_result$x # Fragment length positions
  hist_y <- density_result$y # Density values

  # Extract Valley 1 region (NFR/mononucleosome boundary)
  valley1_mask <- hist_x >= valley1_range[1] & hist_x <= valley1_range[2]
  valley1_x_region <- hist_x[valley1_mask]
  valley1_y_region <- hist_y[valley1_mask]

  # Extract Valley 2 region (mononucleosome/dinucleosome boundary)
  valley2_mask <- hist_x >= valley2_range[1] & hist_x <= valley2_range[2]
  valley2_x_region <- hist_x[valley2_mask]
  valley2_y_region <- hist_y[valley2_mask]

  # Find minimum (valley) in Valley 1 region
  if (length(valley1_y_region) == 0) {
    warning("No data points found in valley1_range [", valley1_range[1], ", ", valley1_range[2], "]")
    valley1_x <- NA
  } else {
    valley1_min_idx <- which.min(valley1_y_region)
    valley1_x <- valley1_x_region[valley1_min_idx]
  }

  # Find minimum (valley) in Valley 2 region
  if (length(valley2_y_region) == 0) {
    warning("No data points found in valley2_range [", valley2_range[1], ", ", valley2_range[2], "]")
    valley2_x <- NA
  } else {
    valley2_min_idx <- which.min(valley2_y_region)
    valley2_x <- valley2_x_region[valley2_min_idx]
  }

  return(c(valley1 = round(valley1_x), valley2 = round(valley2_x)))
}


# Post: Generate fragment length histogram and optionally compute per-pair fragment decomposition statistics based on detected valleys
# Parameter: frags_list: Named list where each element contains fragment lengths for one pair (names are pair IDs)
#            out_dir: Directory path to save output files (histogram PDF and decomposition TSV)
#            BPPARAM: BiocParallel parameter object for parallel processing
#            detect_valley: Logical. If TRUE, detect valleys and compute fragment decomposition statistics (default: FALSE)
#            dens_reso: Resolution for kernel density estimation (default: 2^15)
#            dens_kernel: Kernel type for density estimation (default: "gaussian")
#            valley1_range: Search range for first valley (default: c(73, 221) bp)
#            valley2_range: Search range for second valley (default: c(222, 368) bp)
# Output: Saves histogram PDF to out_dir; if detect_valley=TRUE, also saves per-pair fragment decomposition TSV with counts and percentages
frag_hist <- function(frags_list, out_dir, detect_valley = FALSE, dens_reso = 2^15, dens_kernel = "gaussian", valley1_range = c(73, 221), valley2_range = c(222, 368), BPPARAM = BiocParallel::SerialParam()) {
  # Create output directory if needed
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Combine all fragments for global histogram
  frags_combined <- unlist(frags_list, use.names = FALSE)

  cat("\n=== Fragment Length Analysis ===\n")
  cat(sprintf("Total pairs: %d\n", length(frags_list)))
  cat(sprintf("Total fragments: %s\n", format(length(frags_combined), big.mark = ",")))
  cat(sprintf("Fragment length range: %d - %d bp\n", min(frags_combined), max(frags_combined)))

  # Initialize valley detection results
  valley1 <- NULL
  valley2 <- NULL
  valleys_detected <- FALSE

  # Detect valleys if requested
  if (detect_valley) {
    cat("\nDetecting valley positions from combined distribution...\n")
    tryCatch(
      {
        valley <- find_two_valleys(
          frags = frags_combined,
          dens_reso = dens_reso,
          dens_kernel = dens_kernel,
          valley1_range = valley1_range,
          valley2_range = valley2_range
        )
        valley1 <- valley["valley1"]
        valley2 <- valley["valley2"]

        # Check if both valleys were successfully detected
        if (!is.na(valley1) && !is.na(valley2)) {
          cat(sprintf("  Valley 1 (NFR/mononucleosome):   %d bp\n", valley1))
          cat(sprintf("  Valley 2 (mono/dinucleosome):    %d bp\n", valley2))
          valleys_detected <- TRUE
        } else {
          warning("Valley detection failed - one or both valleys not found")
          detect_valley <- FALSE
        }
      },
      error = function(e) {
        warning("Valley detection error: ", e$message)
        detect_valley <- FALSE
      }
    )
  }

  # Generate histogram plot
  cat("\nGenerating histogram...\n")
  hist_path <- file.path(out_dir, "fragment_distribution.pdf")
  pdf(hist_path, width = 10, height = 6)

  hist(
    frags_combined,
    breaks = 160,
    main = "Fragment Length Distribution",
    xlab = "Fragment Length (bp)",
    ylab = "Frequency"
  )

  # Add valley markers if detected
  if (valleys_detected) {
    abline(v = valley1, col = "red", lwd = 2, lty = 2)
    abline(v = valley2, col = "blue", lwd = 2, lty = 3)
  }

  dev.off()
  cat(sprintf("Histogram saved to: %s\n", hist_path))

  # Generate per-pair fragment decomposition if valleys detected
  if (valleys_detected) {
    cat("\nComputing per-pair fragment decomposition...\n")

    decomp_list <- BiocParallel::bplapply(names(frags_list), function(pair_name) {
      pair_frags <- frags_list[[pair_name]]

      # Classify fragments into three categories based on valleys
      subnucleo_count <- sum(pair_frags < valley1)
      monomer_count <- sum(pair_frags >= valley1 & pair_frags < valley2)
      dimer_plus_count <- sum(pair_frags >= valley2)
      total_count <- length(pair_frags)

      # Calculate percentages
      subnucleo_pct <- subnucleo_count / total_count * 100
      monomer_pct <- monomer_count / total_count * 100
      dimer_plus_pct <- dimer_plus_count / total_count * 100

      data.frame(
        pair = pair_name,
        total_count = total_count,
        subnucleo_count = subnucleo_count,
        subnucleo_pct = subnucleo_pct,
        monomer_count = monomer_count,
        monomer_pct = monomer_pct,
        dimer_plus_count = dimer_plus_count,
        dimer_plus_pct = dimer_plus_pct,
        stringsAsFactors = FALSE
      )
    }, BPPARAM = BPPARAM)
    decomp_df <- data.table::rbindlist(decomp_list)

    # Save decomposition data
    decomp_path <- file.path(out_dir, "fragment_decomposition.tsv")
    write.table(decomp_df, file = decomp_path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("\nDecomposition data saved to: %s\n", decomp_path))
  }
}


# Post: Complete pipeline for fragment length analysis including file reading, quality assessment, valley detection, and decomposition statistics
# Parameter: file_path: Character vector of file paths (BAM or BED format) to analyze
#            out_dir: Directory path to save all output files and plots
#            detect_valley: Logical. If TRUE, perform valley detection and fragment decomposition (default: FALSE)
#            dens_reso: Resolution for kernel density estimation (default: 2^15)
#            dens_kernel: Kernel type for density estimation (default: "gaussian")
#            valley1_range: Search range for first valley (default: c(73, 221) bp)
#            valley2_range: Search range for second valley (default: c(222, 368) bp)
# Output: Generates fragment length histogram PDF and optionally per-pair decomposition TSV; all files saved to out_dir
frag_decomposition <- function(file_path, out_dir = "./", detect_valley = FALSE, dens_reso = 2^15, dens_kernel = "gaussian", valley1_range = c(73, 221), valley2_range = c(222, 368)) {
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(BiocParallel)
    library(data.table)
  })

  # Input validation
  if (!all(file.exists(file_path))) {
    missing_files <- file_path[!file.exists(file_path)]
    stop("Error: The following files do not exist:\n  ", paste(missing_files, collapse = "\n  "))
  }

  # Set up parallel processing
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  cat(sprintf("Using %d CPU cores\n", n_cores))

  if (n_cores > 1) {
    if (.Platform$OS.type == "unix") {
      BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores)
    } else {
      BPPARAM <- BiocParallel::SnowParam(workers = n_cores)
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }

  # Determine file format and extract fragment lengths
  file_exts <- tools::file_ext(file_path)

  cat("\nReading fragment lengths from files...\n")
  if (all(file_exts == "bam")) {
    cat("  File format: BAM (paired-end alignment)\n")
    frags_list <- BiocParallel::bplapply(file_path, function(bam_path) {
      ga <- GenomicAlignments::readGAlignmentPairs(bam_path)
      len <- GenomicRanges::width(GenomicRanges::granges(ga))
      as.integer(len)
    }, BPPARAM = BPPARAM)
  } else if (all(file_exts == "bed")) {
    cat("  File format: BED (genomic intervals)\n")
    frags_list <- BiocParallel::bplapply(file_path, function(bed_path) {
      dt <- data.table::fread(
        bed_path,
        header = FALSE,
        sep = "\t",
        select = 2:3, # Start and end columns
        colClasses = c("NULL", "integer", "integer"),
        data.table = FALSE,
        comment.char = "#"
      )
      len <- dt[[2]] - dt[[1]] # End - start
      as.integer(len)
    }, BPPARAM = BPPARAM)
  } else {
    stop("Error: All files must be either BAM or BED format. Mixed formats are not supported.")
  }
  # Assign pair names to list
  pair_names <- tools::file_path_sans_ext(basename(file_path))
  names(frags_list) <- pair_names

  # Generate histogram and perform analysis
  frag_hist(
    frags_list = frags_list,
    out_dir = out_dir,
    detect_valley = detect_valley,
    dens_reso = dens_reso,
    dens_kernel = dens_kernel,
    valley1_range = valley1_range,
    valley2_range = valley2_range,
    BPPARAM = BPPARAM
  )
}
