# Merge Count Matrices
# Post: Merge multiple count matrices (feather files) by summing values across files,
#       with options for strict consistency checking or flexible merging with zero-filling.
# Parameters:
#   cm_path: A vector of feather file paths to be merged.
#   out_dir: Output file directory for the merged matrix. Default "./"
#   check_consistency: Logical. If TRUE, only keep rows (positions) and columns (samples)
#                      that exist in ALL input files. If FALSE, merge all rows and columns,
#                      filling missing values with 0. Default is TRUE.
# Output: A data frame containing the merged count matrix, with values summed across files.
#         The merged matrix is also saved as a feather file.

merge_count_matrices <- function(cm_path, out_dir = "./", check_consistency = TRUE) {
  library(arrow)

  if (length(cm_path) == 0) {
    stop("No feather files provided")
  }

  all_dfs <- lapply(seq_along(cm_path), function(i) {
    read_feather(cm_path[i])
  })

  if (check_consistency) {
    # Find common positions
    all_pos <- lapply(all_dfs, function(df) df$pos)
    common_pos <- Reduce(intersect, all_pos)
    message("Common positions: ", length(common_pos))
    if (length(common_pos) == 0) {
      stop("Error: No common positions found!")
    }

    # Find common columns
    all_colnames <- lapply(all_dfs, function(df) setdiff(names(df), "pos"))
    common_cols <- Reduce(intersect, all_colnames)
    message("Common columns: ", length(common_cols))
    if (length(common_cols) == 0) {
      stop("Error: No common sample columns found!")
    }

    # Filter and align all dataframes
    aligned_dfs <- lapply(all_dfs, function(df) {
      df_filtered <- df[df$pos %in% common_pos, c("pos", common_cols), drop = FALSE]
      df_filtered <- df_filtered[match(common_pos, df_filtered$pos), ]
      return(df_filtered)
    })

    # Sum all files
    merged_df <- aligned_dfs[[1]]
    if (length(aligned_dfs) > 1) {
      for (i in 2:length(aligned_dfs)) {
        merged_df[, common_cols] <- merged_df[, common_cols] + aligned_dfs[[i]][, common_cols]
      }
    }
  } else {
    # Get all unique positions and columns
    all_pos <- unique(unlist(lapply(all_dfs, function(df) df$pos)))
    all_cols <- unique(unlist(lapply(all_dfs, function(df) setdiff(names(df), "pos"))))
    message("Total positions: ", length(all_pos))
    message("Total columns: ", length(all_cols))

    # Initialize result matrix with zeros
    merged_df <- data.frame(pos = all_pos, stringsAsFactors = FALSE)
    for (col in all_cols) {
      merged_df[[col]] <- 0
    }

    # Add each file using merge
    for (i in seq_along(all_dfs)) {
      current_df <- all_dfs[[i]]
      sample_cols <- setdiff(names(current_df), "pos")
      row_indices <- match(current_df$pos, merged_df$pos)
      for (col in sample_cols) {
        merged_df[row_indices, col] <- merged_df[row_indices, col] + current_df[[col]]
      }
    }
  }

  # Save
  output_path <- file.path(out_dir, "Count_Matrix_merged.feather")
  message("\nSaving to: ", output_path)
  write_feather(merged_df, output_path)
  return(output_path)
}
