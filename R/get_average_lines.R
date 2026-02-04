get_average_lines <- function(pdp_results, target_pair) {
    all_arrays <- list()
    
    seed_names <- names(pdp_results)
    seed_names_sorted <- seed_names[order(as.numeric(seed_names))]
    
    cat("Processing:", seed_names_sorted, "\n")  
    for (random_seed in seed_names_sorted) {
        individual_lines <- pdp_results[[random_seed]][["coding_all"]][[target_pair]][["individual"]]
        grid_length <- length(pdp_results[[random_seed]][["coding_all"]][[target_pair]][["grid_values"]][[1]])
        
        n_samples <- length(individual_lines) / grid_length
        
        if (length(individual_lines) %% grid_length == 0) {
            reshaped_matrix <- matrix(individual_lines, 
                                    nrow = n_samples, 
                                    ncol = grid_length, 
                                    byrow = FALSE)
            
            all_arrays <- append(all_arrays, list(reshaped_matrix))
        }
    }
    
    stacked <- do.call(rbind, all_arrays)
    return(stacked)
}
