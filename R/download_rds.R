# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = TRUE) {
  sha256_list <- c(
    "ChromHMM_hg38.rds" = "6c3fe5c004b5c58b3b6f210612f00db0807ec17508507677dfc8d5bbdfeb7ebb",
    "ChromHMM_mm10.rds" = "fbd8c1ec7407667b590246584c3b6b80292c0db000546d5bea5aaa362e5d4199",
    "ENCODE_cCRE_v4_hg38.rds" = "22175f9af598f9168737f7de0eb8e9bf32dc8fc23c365357f744c04e4e484fac",
    "ENCODE_cCRE_v4_mm10.rds" = "8602da7c586d49aa43b5f4fc70ed81a1e6efba8a1c09517844207f4131823388",
    "GENCODE_v49_hg38_processed.rds" = "f324adaec87b318b6ef0e7398c8609a4e9db7862f78f6c2fd4b15cb98f9e5274",
    "GENCODE_vM23_mm10_processed.rds" = "83534c517d87b6489b5e595465b7ca5681a3e0e41acb55ce7b1e225cf6ce95af",
    "GENCODE_v49_hg38_single_tx_by_evidence.rds" = "b1dc38f708bae8320c6ef3ddf2c680dd4b47f69457490402f8328b1017f82e2d",
    "GENCODE_vM23_mm10_single_tx_by_evidence.rds" = "7c304118d43deb01a173e973c2d112e1be61f2350213f775d4603f3bab69aabd",
    "knownGene_hg38_processed.rds" = "ac82fc657e866f5b34acd01e179559df47c96e58b6890efdf21bcb62554e4fd7",
    "knownGene_mm10_processed.rds" = "c72f10ec5e21942df0a1af2e75b173d9f5a211b6dbf4dbd7208157acc5f25b58",
    "RepeatMasker_hg38_processed.rds" = "3993d786df82edef2a13eb70b10399e2c0725bdedc598a07260b7cbb23da291a",
    "RepeatMasker_mm10_processed.rds" = "8b9068c8d165f9a8651fb0f8e38684375f6e6b63224fc09ffab945f9e0f8e23d",
    "TFBS_lib_hg38.rds" = "0f7f1d2407ad429a9617e11cf71654c7e2af4b2ca43a92ca8537f3cb99540d69"
  )

  # Validate rds_name
  if (!rds_name %in% names(sha256_list)) {
    stop("Unknown file: ", rds_name)
  }

  cache_dir <- rappdirs::user_cache_dir("multiEpiCore")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  zip_name <- sub("\\.rds$", ".zip", rds_name)
  zip_path <- file.path(cache_dir, zip_name)
  rds_path <- file.path(cache_dir, rds_name)

  base_url <- sprintf(
    "https://github.com/Qingjie-Yu/multiEpiCore/releases/download/%s",
    release_tag
  )
  zip_url <- paste0(base_url, "/", zip_name)

  needs_download <- force || !file.exists(rds_path)
  if (!needs_download) {
    sha256_local <- digest::digest(file = rds_path, algo = "sha256")
    if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
      message("Checksum mismatch for ", rds_name, ". Re-downloading...")
      needs_download <- TRUE
      unlink(c(rds_path, zip_path), force = TRUE)
    }
  }

  if (needs_download) {
    max_try <- 3
    success <- FALSE 
    for (i in seq_len(max_try)) {
      tryCatch({
        message("Downloading: ", zip_name, " (attempt ", i, "/", max_try, ")")
        utils::download.file(zip_url, zip_path, mode = "wb", quiet = TRUE)
        utils::unzip(zip_path, exdir = cache_dir)
        if (!file.exists(rds_path)) {
          stop("Unzipped file not found: ", rds_path)
        }
        sha256_local <- digest::digest(file = rds_path, algo = "sha256")
        if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
          stop("sha256 mismatch for ", rds_name, "\nExpected: ", sha256_list[[rds_name]], "\nObserved: ", sha256_local)
        }
        success <- TRUE
        unlink(zip_path, force = TRUE) 
        break
      }, error = function(e) {
        message("Attempt ", i, " failed: ", e$message)
        unlink(c(rds_path, zip_path), force = TRUE)
      })

      if (!success && i < max_try) {
        Sys.sleep(2^i)
      }
    }
    if (!success) {
      stop("Failed to download ", rds_name, " after ", max_try, " attempts.")
    }
  }

  normalizePath(rds_path)
}