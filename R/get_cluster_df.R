# Define a helper function to create a small dataframe for one group.
get_cluster_df <- function(df, prefix) {
  # Get the column names from the given dataframe.
  original_names <- colnames(df)

  # Create a dataframe with:
  # - "feature": prefixed column names, e.g. "V1:colName"
  # - "label": the prefix (group label)
  # - "nonprefix": the original column names (to help with ordering)
  data.frame(
    feature = paste0(prefix, ":", original_names),
    label = prefix,
    nonprefix = original_names,
    stringsAsFactors = FALSE
  )
}
