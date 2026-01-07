# Post: Apply transformation operations to a numeric matrix / data frame, in sequence.
# Parameter: df: A numeric matrix or data frame to be transformed.
#            transformations: A character vector of transformation steps to apply in order.
#            Acceptable values (case sensitive):
#                 - "log2": apply log2(x)
#                 - "t": transpose the matrix
#                 - "Z": Z-score normalization (scale each column)
#                 - "cluster": placeholder (does nothing here)
# Output: A transformed data frame or matrix after applying the specified transformations in order.
transform_df <- function(df, transformations = c("log2")) {
  result <- df
  for (t in transformations) {
    if (t == "log2") {
      result <- log2(result)
    } else if (t == "t") {
      result <- t(result)
    } else if (t == "Z") {
      result <- as.data.frame(scale(result))
    } else {
      warning(paste0("Unrecognized transformation: ", t, "; Only accept log2, t, Z"))
      next
    }
  }
  return(result)
}
