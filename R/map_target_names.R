map_target_names <- function(target_pair_list, target_pair_mapping_df_path = NULL, from = "targets", to = "shorthand") {
  if (is.null(target_pair_mapping_df_path)) {
    data("target_pair_short_hand", package = "multiEpiCore", envir = environment())
    target_pair_mapping_df <- target_pair_short_hand
  } else {
    target_pair_mapping_df <- read.table(target_pair_mapping_df_path, header = TRUE, check.names = FALSE)
  }

  cur_names <- target_pair_mapping_df[[from]]
  new_names <- target_pair_mapping_df[[to]]
  result <- target_pair_list
  for (i in seq_along(cur_names)) {
    cur_name <- cur_names[i]
    new_name <- new_names[i]
    result <- gsub(cur_name, new_name, result, fixed = TRUE)
  }
  return(result)
}
