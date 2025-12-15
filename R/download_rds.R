# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = FALSE) {
  sha256_list <- c(
    "ChromHMM_hg38.rds" = "1dc13278981bab4eded9e35423404c9d8f77140f108b275a92bffaad6593b992",
    "GENCODE_v49_hg38.rds" = "59160b125527d5b4e9621b23bee9d323b056b93148921d563669cbed6b858ff6",
    "ENCODE_cCRE_v4_hg38.rds" = "d24441998400770947e0c8c760251c408cb870b329f4445f3201f500e962a9ad",
    "ENCODE_cCRE_v4_mm10.rds" = "42f3e5af98768b505d976538c006a4cbf6e3099f9e3ca3490ef5c9b9b15922a1",
    "ChromHMM_mm10.rds" = "7d7f4007093bfaae0512671e487048e33c22f7962f350fd5553e5d4dad7a69dd",
    "GENCODE_v49_hg38_processed.rds" = "8d08eb97c599075de5aca1c2617146397eb90fee290a1ede7e3b34f41ad466f5",
    "GENCODE_v49_hg38_single_tx_by_evidence.rds" = "a202ee1cdde557b922cd4a02fd0175514555103ab6169059277dd8cc1efce318",
    "GENCODE_vM23_mm10.rds" = "5ae83defc0160ef9354c00d4e5ddeb5a100702ca9c56113f26b3d5445a899e9f",
    "GENCODE_vM23_mm10_processed.rds" = "6f366f7a0d79d48cee07ce4bccf4bea5b933272ddbbd133f0b548387f8879a7c",
    "GENCODE_vM23_mm10_single_tx_by_evidence.rds" = "1b69a9e161b6d6dd99504c97cec945ac697e320b170db39947aa7d4c860ca7c0",
    "knownGene_hg38_processed.rds" = "443c55840283bb34814c5e3b6b6ecd6bde9efd2e475fb27ad5d050af827dce19",
    "knownGene_mm10_processed.rds" = "268a4fe8dcc197d75ea1b4e72252cf8b08437d5a86e47488776ed230aeaf26fe",
    "TFBS_lib_hg38.rds" = "4ee4b749cd599b4d1090dbf938c6a5081965bcae956d8340a417472ab76c6981"
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
    sha256_local <- digest::digest(zip_path, algo = "sha256", serialize = FALSE)
    if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
      stop("sha256 mismatch for ", rds_name, "\nExpected: ", sha256_list[[rds_name]], "\nObserved: ", sha256_local)
    }
    utils::unzip(zip_path, exdir = cache_dir)
  }
  
  normalizePath(rds_path)
}