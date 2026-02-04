map_target_names_pdp <- function(target_pair_list, target_pair_mapping_df, from_col = "targets", to_col = "shorthand") {
    cur_names <- target_pair_mapping_df[[from_col]]
    new_names <- target_pair_mapping_df[[to_col]]
    
    cat("\n")
    result <- target_pair_list
    for (i in seq_along(cur_names)) {
        cur_name <- cur_names[i]
        new_name <- new_names[i]
        result <- sapply(result, function(target_pair) {
            gsub(cur_name, new_name, target_pair, fixed = TRUE)
        }, USE.NAMES = FALSE)
    }
    
    return(result)
}
