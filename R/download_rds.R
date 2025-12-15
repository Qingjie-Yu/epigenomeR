# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = FALSE) {
  sha256_list <- c(
    "ChromHMM_hg38.rds" = "6ce5972399d2b712ab1f76722ef98e8573dc6e2366afd4dc10a670b016973edf",
    "GENCODE_v49_hg38.rds" = "279993a6862ae4487a4f5caba5560ca2bb298ad74c40a5ff5c23b5d23534f354",
    "ENCODE_cCRE_v4_hg38.rds" = "21e99022ed54b2c8edb106606a107fb9ed7eab0a480a5512ca50e28e4e81916d",
    "ENCODE_cCRE_v4_mm10.rds" = "2d43916cd954bd86806c9287ce2ce0f573231e214dee531e6efd5d2d72a89c08",
    "ChromHMM_mm10.rds" = "d7bcf63a4e95e9f3c68c36358f583ae4b2492b6ebda8a9ff8834f60aebf69a0b",
    "GENCODE_v49_hg38_processed.rds" = "3b946ab7869b0bd48f902ab6c95582c1a7b4d60b19e0d304d557a6d7fcef5cf4",
    "GENCODE_v49_hg38_single_tx_by_evidence.rds" = "fddab927b927c2874085645df12aa2a30d6ed9c9e78bb7dc714144a695641b3a",
    "GENCODE_vM23_mm10.rds" = "8163c36daf7be227620c93692500058331ef21ebb86a57acb09d246745c8252d",
    "GENCODE_vM23_mm10_processed.rds" = "1ae21ca2f4289389a63ce4331686b6ca3d93e4340188650db80352bf96466585",
    "GENCODE_vM23_mm10_single_tx_by_evidence.rds" = "f3e8429e790420a877dd9df72f70dd0985ad4ed41a429076414fcb64ff9928d7","knownGene_hg38_processed.rds" = "6946011ac979e0fe09d2a27427099419b6b0897b15e5d3d69be77acc98aec895",
    "knownGene_mm10_processed.rds" = "bfcb4bc634eced7451a5226cfbdf9fb6f0dc4041d03638037a978b7eebfc8299",
    "TFBS_lib_hg38.rds" = "59b2a61d9a78caef07e5b7a82fac03aa24266e1b5a5281ec16f15ba4711bbddc"
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
      stop("sha256 mismatch for ", zip_name, "\nExpected: ", sha256_list[[rds_name]], "\nObserved: ", sha256_local)
    }
    utils::unzip(zip_path, exdir = cache_dir)
  }
  
  normalizePath(rds_path)
}