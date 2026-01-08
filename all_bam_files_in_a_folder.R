all_bam_files_in_a_folder <- function(bam_dir, want = c("bam_files_path", "bam_files_name")) {
  bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", recursive = TRUE, full.names = TRUE)
  want <- match.arg(want) 
  
  if (want == "bam_files_path") { 
    return(bam_files) 
  } else if (want == "bam_files_name") { 
    bam_names <- basename(bam_files)
    bam_names_noext <- tools::file_path_sans_ext(bam_names)
    return(bam_names_noext)  
  }
}
