# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = FALSE) {
  cache_dir <- rappdirs::user_cache_dir("epigenomeR")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  zip_name <- sub("\\.rds$", ".zip", rds_name)
  zip_path <- file.path(cache_dir, zip_name)
  rds_path <- file.path(cache_dir, rds_name)

  base_url <- sprintf(
    "https://github.com/Qingjie-Yu/epigenomeR/releases/download/%s",
    release_tag
  )
  zip_url <- paste0(base_url, "/", zip_name)

  if (!force && file.exists(rds_path)) {
    return(normalizePath(rds_path))
  }
  utils::unzip(zip_path, exdir = cache_dir)
  normalizePath(rds_path)
}