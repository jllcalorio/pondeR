#' Scale Metabolomics Data
#'
#' @description
#' Applies various scaling methods to standardize features and prepare data
#' for multivariate analysis.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param method Character. Scaling method to apply. Options:
#'   \itemize{
#'     \item \code{"mean"}: Mean-centering only (centers distribution at zero)
#'     \item \code{"auto"}: Auto-scaling (mean-centering + unit variance scaling)
#'     \item \code{"pareto"}: Pareto-scaling (mean-centering + sqrt(SD) scaling)
#'   }
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_scale" containing:
#'   \item{data}{Matrix or data frame of scaled data (same class as input `x`)}
#'   \item{scaling_factors}{Numeric vector of scaling factors applied to each feature}
#'   \item{center_values}{Numeric vector of centering values (means) for each feature}
#'   \item{method_used}{Character describing scaling method applied}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Scaling Methods:**
#' 
#' - **mean**: Mean-centering only. Each feature is centered to have mean = 0
#'   but retains original variance. Useful when features are already on similar
#'   scales and you want to preserve relative variance information.
#' 
#' - **auto**: Auto-scaling (also called standardization or unit variance scaling).
#'   Each feature is mean-centered and divided by its standard deviation,
#'   resulting in mean = 0 and SD = 1 for all features. This gives equal weight
#'   to all features regardless of their original variance.
#'   
#'   **Recommended for:**
#'   - Principal Component Analysis (PCA)
#'   - Hierarchical clustering
#'   - Correlation analysis
#'   - t-tests and ANOVA
#'   - Most univariate and multivariate methods
#'   
#'   Auto-scaling can amplify noise in low-variance features.
#' 
#' - **pareto**: Pareto-scaling. Each feature is mean-centered and divided by
#'   the square root of its standard deviation. This provides a compromise between
#'   no scaling and auto-scaling, reducing the influence of large values while
#'   preserving some of the original variance structure.
#'   
#'   **Recommended for:**
#'   - PLS-DA (Partial Least Squares Discriminant Analysis)
#'   - OPLS-DA (Orthogonal PLS-DA)
#'   - sPLS-DA (Sparse PLS-DA)
#'   - Other PLS-type methods
#'   
#'   Pareto-scaling is generally preferred for discriminant analysis as it
#'   balances the importance of large and small features better than auto-scaling.
#' 
#' **Choosing a Scaling Method:**
#' 
#' The choice depends on your downstream analysis:
#' - For exploratory data analysis (PCA, clustering): **auto-scaling**
#' - For classification/discrimination (PLS-DA): **Pareto-scaling**
#' - When feature magnitudes carry meaning: **mean-centering**
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., & van der Werf, M.J. (2006).
#' Centering, scaling, and transformations: improving the biological information content of
#' metabolomics data. BMC Genomics, 7, 142. \doi{10.1186/1471-2164-7-142}
#' 
#' Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). The New S Language.
#' Wadsworth & Brooks/Cole. \doi{10.1201/9781351074988}
#'
#' @importFrom matrixStats colSds
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 50, mean = 100, sd = 50),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' 
#' # Auto-scaling for PCA
#' result1 <- run_scale(x, method = "auto")
#' 
#' # Pareto-scaling for PLS-DA
#' result2 <- run_scale(x, method = "pareto")
#' 
#' # Mean-centering only
#' result3 <- run_scale(x, method = "mean")
#' }
run_scale <- function(
    x,
    method = "auto",
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  
  method <- tolower(method)
  valid_methods <- c("mean", "auto", "pareto")
  if (!method %in% valid_methods) {
    stop("Unknown scaling method '", method, "'. ",
         "Supported: ", paste(valid_methods, collapse = ", "))
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  
  msg(sprintf("Applying '%s' scaling...", method))
  
  # Calculate column means (always needed)
  col_means <- colMeans(x_matrix, na.rm = TRUE)
  
  x_scaled <- switch(method,
                     
                     "mean" = {
                       msg("Applying mean-centering...")
                       scaling_factors <- rep(1, ncol(x_matrix))
                       scale(x_matrix, center = TRUE, scale = FALSE)
                     },
                     
                     "auto" = {
                       msg("Applying auto-scaling (mean-centering + unit variance)...")
                      #  col_sds <- matrixStats::colSds(x_matrix, na.rm = TRUE)
                      #  col_sds[col_sds == 0] <- 1  # Avoid division by zero
                      #  scaling_factors <- col_sds
                       scale(x_matrix, center = TRUE, scale = TRUE)
                     },
                     
                     "pareto" = {
                       msg("Applying Pareto-scaling (mean-centering + sqrt(SD))...")
                       col_sds <- matrixStats::colSds(x_matrix, na.rm = TRUE)
                       col_sds[col_sds == 0] <- 1  # Avoid division by zero
                       scaling_factors <- sqrt(col_sds)
                       scale(x_matrix, center = TRUE, scale = sqrt(col_sds))
                     }
  )
  
  msg("Scaling complete.")
  
  # Remove scale attributes to return clean data
  attributes(x_scaled) <- list(dim = dim(x_scaled), dimnames = dimnames(x_scaled))
  
  # Return in original class
  if (!was_matrix) {
    x_scaled <- as.data.frame(x_scaled)
  }
  
  result <- list(
    data = x_scaled,
    # scaling_factors = scaling_factors,
    center_values = col_means,
    method_used = method,
    parameters = list(method = method)
  )
  
  class(result) <- c("run_scale", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_scale <- function(x, ...) {
  cat("=== Scaling Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Data dimensions: %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  invisible(x)
}