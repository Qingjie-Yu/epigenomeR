remove_model_in_param_grid <- function(best_params) {
    results <- list()
    for (param in names(best_params)) {
        new_param <- gsub("model__", "", param, fixed = TRUE)
        results[[new_param]] <- best_params[[param]]
    }
    return(results)
}
