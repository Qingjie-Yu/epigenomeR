safe_minmax <- function(df) {
    df_scaled <- df
    for (col in colnames(df)) {
        col_min <- min(df[[col]], na.rm = TRUE)
        col_max <- max(df[[col]], na.rm = TRUE)
        if (col_max != col_min) {
        df_scaled[[col]] <- (df[[col]] - col_min) / (col_max - col_min)
        } else {
        df_scaled[[col]] <- 0.0
        }
    }
    return(df_scaled)
}
