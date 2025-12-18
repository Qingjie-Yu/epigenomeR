filter_target_pairs <- function(percentage_cutoff = 0.25, target_pairs=NULL, frag_len_num_file=NULL) {
  if (is.null(target_pair_mapping_df_path)) {
    data("frag_len_num", package = "multiEpiCore", envir = environment())
  } else {
    frag_len_num <- read.table(frag_len_num_file, sep = "\t", header = TRUE, row.names = 1)
  }

  frag_len_num$total_frag <- rowSums(frag_len_num)
  percentile_res <- quantile(frag_len_num$total_frag, type=3, probs = c(percentage_cutoff))
  # print(percentile_res)
  frag_len_num_cutoff <- unname(percentile_res)
  frag_len_num_filtered <- frag_len_num[frag_len_num$total_frag >= frag_len_num_cutoff, ]
  filtered_target_pairs <- rownames(frag_len_num_filtered)
  if (!is.null(target_pairs)) {
    filtered_target_pairs <- filtered_target_pairs[filtered_target_pairs %in% target_pairs]
  }
  return(filtered_target_pairs)
}
