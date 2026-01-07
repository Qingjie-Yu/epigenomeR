list_to_dataframe <- function(dataList) {
  if (is.null(names(dataList))) {
    return(do.call("rbind", dataList))
  }

  cn <- lapply(dataList, colnames) %>%
    unlist() %>%
    unique()
  cn <- c(".id", cn)
  dataList2 <- lapply(seq_along(dataList), function(i) {
    data <- dataList[[i]]
    data$.id <- names(dataList)[i]
    idx <- !cn %in% colnames(data)
    if (sum(idx) > 0) {
      for (i in cn[idx]) {
        data[, i] <- NA
      }
    }
    return(data[, cn])
  })
  res <- do.call("rbind", dataList2)
  res$.id <- factor(res$.id, levels = rev(names(dataList)))
  return(res)
}
