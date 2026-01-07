all_files_in_a_folder <- function(dir, ext = ".bam", want = c("files_path", "files_name"), filter_IgG = TRUE) {
  files <- list.files(path = dir, pattern = paste0("\\", ext, "$"), recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("unknown|is", files, ignore.case = TRUE)]

  if (filter_IgG == TRUE) {
    files <- files[!grepl("IgG", files, ignore.case = TRUE)]
  }
  want <- match.arg(want)

  if (want == "files_path") {
    return(files)
  } else if (want == "files_name") {
    names <- basename(files)
    names_noext <- tools::file_path_sans_ext(names)
    return(names_noext)
  }
}
