
get_average_lines <- function(pdp_results, target_pair) {
    all_arrays <- list()
    
    # grid length
    first_seed <- names(pdp_results)[1]
    grid_length <- length(pdp_results[[first_seed]][["coding_all"]][[target_pair]][["grid_values"]][[1]])
    
    for (random_seed in names(pdp_results)) {
        individual_lines <- pdp_results[[random_seed]][["coding_all"]][[target_pair]][["individual"]]
        
        n_samples <- length(individual_lines) / grid_length
        
        if (length(individual_lines) %% grid_length == 0) {
            reshaped_matrix <- matrix(individual_lines, 
                                    nrow = n_samples, 
                                    ncol = grid_length, 
                                    byrow = TRUE)
            all_arrays <- append(all_arrays, list(reshaped_matrix))
        } else {
            warning(paste("Cannot reshape data for seed", random_seed))
        }
    }
    
    if (length(all_arrays) > 0) {
        stacked <- do.call(rbind, all_arrays)
        return(stacked)
    } else {
        stop("Could not process individual_lines data")
    }
}
