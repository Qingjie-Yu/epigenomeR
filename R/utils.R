get_BPPARAM <- function(verbose = TRUE, backend = c("snow", "multicore", "serial")) {
  backend <- match.arg(backend)

  slurm_cores <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
  )

  if (is.na(slurm_cores) || slurm_cores < 1) {
    detected <- parallel::detectCores()
    n_cores <- if (is.na(detected) || detected <= 1) 1L else detected - 1L
  } else {
    n_cores <- slurm_cores
  }

  n_cores <- max(1L, n_cores)

  if (backend == "snow") {
    if (n_cores > 1) {
      bp <- BiocParallel::SnowParam(workers = n_cores, type = "SOCK")
    } else {
      bp <- BiocParallel::SerialParam()
    }

  } else if (backend == "multicore") {
    if (.Platform$OS.type != "unix") {
      stop("MulticoreParam is only supported on unix-like systems.")
    }
    if (n_cores > 1) {
      bp <- BiocParallel::MulticoreParam(workers = n_cores)
    } else {
      bp <- BiocParallel::SerialParam()
    }

  } else {
    bp <- BiocParallel::SerialParam()
  }

  if (verbose) {
    message(sprintf("Using %d CPU core(s) for BiocParallel.", n_cores))
    message(sprintf("Backend: %s", class(bp)[1]))
  }

  return(bp)
}