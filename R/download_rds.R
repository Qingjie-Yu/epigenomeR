# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = FALSE) {
  sha256_list <- c(
    "ChromHMM_hg38.rds" = "7c94b0691e15b318a42cf6316ae9486fbf9ef1325e6baeb3cfa4fbe9fcbc975a",
    "ChromHMM_mm10.rds" = "e44aa16d04cdcf5399a751695c030170b3870c970146d7696dea98afd59e01f5",
    "ENCODE_cCRE_v4_hg38.rds" = "540ea690b6df3fc0116cdf98155d459164ccd3d2e5f68f0f134464b91a471c80",
    "ENCODE_cCRE_v4_mm10.rds" = "ab0cfb2270776f7bdf7236887eed4d498e4c0d77037286ea7d27338af3b3de41",
    "GENCODE_v49_hg38_processed.rds" = "856c8c96b050b8cf3b05c87ff3351b673f1a51c2768507b05b7137f79a749a2d",
    "GENCODE_vM23_mm10_processed.rds" = "b1bc36df63f043fb1cda950d6af247cd79537b1f38d106eb4894846a429569e5","knownGene_hg38_processed.rds" = "7248c3d519f953f96d6c2ca16643ae02d2d9f922fe4f2bcefa1071cbff5b628f",
    "knownGene_mm10_processed.rds" = "61fd67732af331069f7a8151870460b23339b2717a807fa691c4795c78f831a7",
    "RepeatMasker_hg38_processed.rds" = "1ba1ae0c85cf871be32868ef4b21438ea630de44614101cf480a47e942765ee5",
    "RepeatMasker_mm10_processed.rds" = "87830e07f1abf7fd69daae2f45e027b547f223946e357ea167c094ee40476c0a",
    "TFBS_lib_hg38.rds" = "51d5191c9946fa2896e016cd2bd1e6b232021715d2676979b70b1d1c7ba68d0c"
  )

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

  if (force || !file.exists(rds_path)) {
    message("Downloading: ", zip_name)
    utils::download.file(zip_url, zip_path, mode = "wb", quiet = TRUE)
    utils::unzip(zip_path, exdir = cache_dir)
  }
  sha256_local <- digest::digest(rds_path, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
    stop("sha256 mismatch for ", rds_name, "\nExpected: ", sha256_list[[rds_name]], "\nObserved: ", sha256_local)
  }
  
  normalizePath(rds_path)
}