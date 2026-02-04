
get_BPPARAM <- function(verbose = TRUE){
  slurm_cores <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
  )

  if (is.na(slurm_cores) || slurm_cores < 1){
    detected <- parallel::detectCores()
    n_cores <- if (is.na(detected) || detected <= 1) 1L else detected - 1L
  } else {
    n_cores <- slurm_cores
  }

  n_cores <- max(1L, n_cores)

  if (n_cores > 1){
    if (.Platform$OS.type == "unix") {
      bp <- BiocParallel::MulticoreParam(workers = n_cores)
    } else {
      bp <- BiocParallel::SnowParam(workers = n_cores)
    }
  } else {
    bp <- BiocParallel::SerialParam()
  }

  if (verbose){
    message(sprintf("Using %d CPU core(s) for BiocParallel.", n_cores))
  }

  bp
}