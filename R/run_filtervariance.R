#' Filter Features by Low Variance
#'
#' @description
#' Removes features with low variance across samples. Features with minimal
#' variation provide little information for downstream analysis and can be
#' safely removed to reduce dimensionality.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param percentile Numeric. Percentile threshold (0-100) for variance filtering.
#'   Features in the bottom X percentile of variance are removed. For example,
#'   10 removes the 10% of features with lowest variance. Default: 10.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_filtervariance" containing:
#'   \item{data}{Matrix or data frame of filtered data (same class as input `x`)}
#'   \item{features_removed}{Character vector of removed feature names}
#'   \item{features_kept}{Character vector of retained feature names}
#'   \item{variance_values}{Numeric vector of variance values for all original features}
#'   \item{variance_threshold}{Numeric. The variance value corresponding to the percentile cutoff}
#'   \item{n_features_before}{Integer. Number of features before filtering}
#'   \item{n_features_after}{Integer. Number of features after filtering}
#'   \item{n_features_removed}{Integer. Number of features removed}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Why Filter by Variance:**
#' 
#' Features with very low variance across samples:
#' - Provide minimal discriminatory power for classification
#' - Contribute little to multivariate models (PCA, PLS-DA)
#' - May represent technical noise or detection artifacts
#' - Increase computational burden without adding information
#' 
#' **Choosing the Percentile:**
#' 
#' - **Conservative** (5-10%): Remove only the least variable features
#' - **Moderate** (10-20%): Standard approach for most studies
#' - **Aggressive** (20-30%): When high dimensionality is a concern
#' 
#' The default (10%) removes the bottom 10% of features by variance, which
#' is a commonly used threshold in metabolomics.
#' 
#' **Important Notes:**
#' 
#' 1. Variance filtering should be applied **after** normalization and
#'    transformation, as these steps affect feature variance.
#'    
#' 2. The variance calculation uses all samples (including QCs if present).
#'    Consider removing QC samples before variance filtering if they have
#'    artificially low variance.
#'    
#' 3. For scaled data (auto or Pareto), this filter identifies features with
#'    inherently low biological variation across the study groups.
#'    
#' 4. Very high percentiles (>30%) risk removing biologically relevant features
#'    with naturally stable expression levels.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control, Batch Correction,
#' Exploration & Data Cleaning (Version 2.1) Zendono. \doi{10.5281/zenodo.16824822}.
#' Retrieved from \url{https://github.com/broadhurstdavid/QC-MXP}.
#'
#' @importFrom matrixStats colVars
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with varying levels of variance
#' set.seed(123)
#' n_samples <- 100
#' n_features <- 50
#' 
#' # Create features with different variance levels
#' x <- matrix(nrow = n_samples, ncol = n_features)
#' for (i in 1:n_features) {
#'   # Variance increases with feature index
#'   sd_val <- i / 2
#'   x[, i] <- rnorm(n_samples, mean = 100, sd = sd_val)
#' }
#' colnames(x) <- paste0("Feature", 1:n_features)
#' 
#' # Remove bottom 10% (5 features with lowest variance)
#' result1 <- run_filtervariance(x, percentile = 10)
#' 
#' # More aggressive filtering (bottom 20%)
#' result2 <- run_filtervariance(x, percentile = 20)
#' }
run_filtervariance <- function(
    x,
    percentile = 10,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  if (!is.numeric(percentile) || length(percentile) != 1 || 
      percentile < 0 || percentile > 100) {
    stop("'percentile' must be a numeric value between 0 and 100")
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  feature_names <- colnames(x_matrix)
  n_features_before <- ncol(x_matrix)
  
  msg(sprintf("Filtering features by variance (%gth percentile cutoff)...", percentile))
  
  # Calculate variance for each feature
  variances <- matrixStats::colVars(x_matrix, na.rm = TRUE)
  names(variances) <- feature_names
  
  # Handle features with zero or invalid variance
  valid_var_indices <- is.finite(variances) & variances > 0
  
  if (sum(valid_var_indices) == 0) {
    warning("No features with valid positive variance. Returning original data.")
    result <- list(
      data = if (was_matrix) x_matrix else as.data.frame(x_matrix),
      features_removed = character(0),
      features_kept = feature_names,
      variance_values = variances,
      variance_threshold = NA,
      n_features_before = n_features_before,
      n_features_after = n_features_before,
      n_features_removed = 0,
      parameters = list(percentile = percentile)
    )
    class(result) <- c("run_filtervariance", "list")
    return(result)
  }
  
  # Calculate variance threshold at specified percentile
  var_threshold <- quantile(variances[valid_var_indices], percentile / 100, na.rm = TRUE)
  
  # Keep features with variance > threshold OR invalid variance (to be safe)
  features_to_keep <- (variances > var_threshold) | !valid_var_indices
  
  x_filtered <- x_matrix[, features_to_keep, drop = FALSE]
  features_kept <- feature_names[features_to_keep]
  features_removed <- feature_names[!features_to_keep]
  
  n_features_after <- length(features_kept)
  n_removed <- n_features_before - n_features_after
  
  msg(sprintf("Removed %d features (%.1f%%) with variance <= %.6f",
              n_removed, (n_removed / n_features_before) * 100, var_threshold))
  
  # Return in original class
  if (!was_matrix) {
    x_filtered <- as.data.frame(x_filtered)
  }
  
  result <- list(
    data = x_filtered,
    features_removed = features_removed,
    features_kept = features_kept,
    variance_values = variances,
    variance_threshold = var_threshold,
    n_features_before = n_features_before,
    n_features_after = n_features_after,
    n_features_removed = n_removed,
    parameters = list(percentile = percentile)
  )
  
  class(result) <- c("run_filtervariance", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_filtervariance <- function(x, ...) {
  cat("=== Variance Filtering Results ===\n")
  cat(sprintf("Features before: %d\n", x$n_features_before))
  cat(sprintf("Features after:  %d\n", x$n_features_after))
  cat(sprintf("Features removed: %d (%.1f%%)\n", 
              x$n_features_removed,
              (x$n_features_removed / x$n_features_before) * 100))
  cat(sprintf("Percentile cutoff: %g%%\n", x$parameters$percentile))
  cat(sprintf("Variance threshold: %.6f\n", x$variance_threshold))
  invisible(x)
}