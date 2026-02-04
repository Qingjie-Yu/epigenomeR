calculate_effect_size <- function(x, y, method = "slope") {  
    if (method == "slope") {
        result <- (y[length(y)] - y[1]) / (x[length(x)] - x[1])
        
    } else if (method == "auc") {
        sorted_indices <- order(x)
        x_sorted <- x[sorted_indices]
        y_sorted <- y[sorted_indices]
        result <- pracma::trapz(x_sorted, y_sorted)
    } else if (method == "end") {
        result <- y[length(y)]
    }
    return(result)
}
