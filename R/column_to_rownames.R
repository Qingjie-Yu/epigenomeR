column_to_rownames <- function(wgc, var = "pos") {
    pos_list_full <- wgc[[var]]
    rownames(wgc) <- wgc[[var]]
    wgc <- wgc[, !names(wgc) %in% var, drop = FALSE]
    return(wgc)
}
