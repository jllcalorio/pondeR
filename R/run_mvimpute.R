#' Impute Missing Values in Metabolomics Data
#'
#' @description
#' Replaces missing values (NAs) using either a simple deterministic method
#' (fraction of minimum) or a statistical imputation method (quantile regression).
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   Missing values should be represented as NA.
#' @param method Numeric or character. If numeric, represents the fraction of the minimum
#'   positive value per feature to use for imputation (e.g., 0.2 = 1/5th of minimum).
#'   If character, must be one of the available imputation methods. Currently supported:
#'   "quantileregression" (uses QRILC algorithm). Default: 0.2 (1/5th of minimum).
#' @param tune_sigma Numeric. Only used when `method = "quantileregression"`.
#'   Controls the standard deviation of the left-censored missing value distribution.
#'   Smaller values (e.g., 0.5) assume tighter distribution; larger values (e.g., 2)
#'   allow more spread. Default: 1.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_mvimpute" containing:
#'   \item{data}{Matrix or data frame with imputed values (same class as input `x`)}
#'   \item{n_missing_before}{Integer. Number of missing values before imputation}
#'   \item{n_missing_after}{Integer. Number of missing values after imputation (should be 0)}
#'   \item{imputed_summary}{Data frame with per-feature imputation statistics}
#'   \item{parameters}{List of parameters used}
#'   \item{method_used}{Character describing the imputation method applied}
#'
#' @details
#' **Missing Value Imputation Strategies:**
#' 
#' 1. **Deterministic (fraction of minimum)**: When `method` is numeric, each missing
#'    value in a feature is replaced by the minimum positive value in that feature
#'    multiplied by the specified fraction. This assumes missing values represent
#'    low abundance (Missing Not At Random - MNAR) and is common in metabolomics
#'    for values below detection limit.
#'    
#'    - **Recommended range**: 0.01 to 0.5 (1% to 50% of minimum)
#'    - Values > 1 will trigger a warning as they exceed the observed minimum
#'    - Default (0.2) replaces NAs with 1/5th of the feature's minimum
#'
#' 2. **Quantile Regression (QRILC)**: When `method = "quantileregression"`, uses
#'    the Quantile Regression Imputation of Left-Censored data algorithm. This
#'    method assumes missing values are Missing Not At Random (MNAR) due to low
#'    abundance and models them from a left-censored distribution.
#'    
#'    - The `tune_sigma` parameter controls the spread of imputed values
#'    - Requires the `imputeLCMD` package
#'    - More sophisticated than simple fraction-of-minimum
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
#' Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data.
#' Scientific Reports, 8(1), 663. \doi{10.1038/s41598-017-19120-0}
#' 
#' Lazar, C., Gatto, L., Ferro, M., Bruley, C., & Burger, T. (2016).
#' Accounting for the Multiple Natures of Missing Values in Label-Free Quantitative
#' Proteomics Data Sets to Compare Imputation Strategies.
#' Journal of Proteome Research, 15(4), 1116-1125. \doi{10.1021/acs.jproteome.5b00981}
#'
#' @importFrom imputeLCMD impute.QRILC
#' @importFrom matrixStats colMins
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with missing values
#' set.seed(123)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
#' x[sample(length(x), 200)] <- NA
#' colnames(x) <- paste0("Feature", 1:50)
#' 
#' # Simple imputation (1/5th of minimum)
#' result1 <- run_mvimpute(x, method = 0.2)
#' 
#' # Quantile regression imputation
#' result2 <- run_mvimpute(x, method = "quantileregression", tune_sigma = 1)
#' }
run_mvimpute <- function(
    x,
    method = 0.2,
    tune_sigma = 1,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  
  n_missing_before <- sum(is.na(x_matrix))
  
  if (n_missing_before == 0) {
    msg("No missing values detected. Returning original data.")
    result <- list(
      data = if (was_matrix) x_matrix else as.data.frame(x_matrix),
      n_missing_before = 0,
      n_missing_after = 0,
      imputed_summary = data.frame(
        Feature = colnames(x_matrix),
        N_Missing = 0,
        Imputation_Value = NA
      ),
      parameters = list(method = method, tune_sigma = tune_sigma),
      method_used = "None (no missing values)"
    )
    class(result) <- c("run_mvimpute", "list")
    return(result)
  }
  
  msg(sprintf("Starting missing value imputation (%d missing values)...", n_missing_before))
  
  # Imputation strategy
  if (is.numeric(method)) {
    
    if (length(method) != 1) {
      stop("'method' when numeric must be a single value")
    }
    if (method <= 0) {
      stop("'method' when numeric must be positive")
    }
    if (method > 1) {
      warning("'method' > 1 means imputed values will exceed the observed minimum. ",
              "This is unusual for MNAR imputation. Consider using a value < 1.")
    }
    
    msg(sprintf("Using deterministic imputation: %.2f x minimum value per feature", method))
    
  # # Calculate minimum per column (ignoring NAs)
  # min_vals <- matrixStats::colMins(x_matrix, na.rm = TRUE)
  
  # # Handle features that were all NA (min will be Inf)
  # min_vals[!is.finite(min_vals)] <- 1e-9
  
  # # Impute: NA -> fraction of minimum
  # imputation_vals <- min_vals * method
  # na_indices <- is.na(x_matrix)
  # x_matrix[na_indices] <- imputation_vals[col(x_matrix)][na_indices]
  
  # method_description <- sprintf("Deterministic: %.2f x min", method)
  
  # imputed_summary <- data.frame(
  #   Feature = colnames(x_matrix),
  #   N_Missing = colSums(na_indices),
  #   Min_Value = min_vals,
  #   Imputation_Value = imputation_vals
  # )

  # Calculate minimum per column (ignoring NAs)
  min_vals <- matrixStats::colMins(x_matrix, na.rm = TRUE)

  # Handle features that were all NA (min will be Inf)
  min_vals[!is.finite(min_vals)] <- 1e-9

  # Apply the fraction
  min_vals <- min_vals * method

  # Impute: NA -> fraction of minimum (direct indexing)
  na_indices <- is.na(x_matrix)
  x_matrix[na_indices] <- min_vals[col(x_matrix)][na_indices]

  method_description <- sprintf("Deterministic: %.2f x min", method)

  imputed_summary <- data.frame(
    Feature = colnames(x_matrix),
    N_Missing = colSums(na_indices),
    Min_Value = min_vals / method,  # Store original min
    Imputation_Value = min_vals
  )
    
  } else if (is.character(method)) {
    
    method <- tolower(method)
    
    if (method == "quantileregression") {
      
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' is required for quantile regression imputation. ",
             "Install it with: BiocManager::install('imputeLCMD')")
      }
      
      if (!is.numeric(tune_sigma) || length(tune_sigma) != 1 || tune_sigma <= 0) {
        stop("'tune_sigma' must be a positive numeric value")
      }
      
      msg(sprintf("Using quantile regression imputation (QRILC) with tune_sigma = %.2f", tune_sigma))
      
      tryCatch({
        # QRILC expects features in rows
        x_imputed <- imputeLCMD::impute.QRILC(t(x_matrix), tune.sigma = tune_sigma)
        x_matrix <- t(x_imputed[[1]])
        
        method_description <- sprintf("QRILC (tune_sigma = %.2f)", tune_sigma)
        
        imputed_summary <- data.frame(
          Feature = colnames(x_matrix),
          N_Missing = colSums(is.na(as.matrix(x))),
          Method = "QRILC"
        )
        
      }, error = function(e) {
        stop("Quantile regression imputation failed: ", e$message,
             "\nConsider using a deterministic method instead.")
      })
      
    } else {
      stop("Unknown imputation method '", method, "'. ",
           "Supported methods: 'quantileregression' or a numeric fraction.")
    }
    
  } else {
    stop("'method' must be either numeric (fraction of minimum) or character (method name)")
  }
  
  n_missing_after <- sum(is.na(x_matrix))
  
  msg(sprintf("Imputation complete. Missing values: %d -> %d", 
              n_missing_before, n_missing_after))
  
  # Return in original class
  if (!was_matrix) {
    x_matrix <- as.data.frame(x_matrix)
  }
  
  result <- list(
    data = x_matrix,
    n_missing_before = n_missing_before,
    n_missing_after = n_missing_after,
    imputed_summary = imputed_summary,
    parameters = list(
      method = method,
      tune_sigma = if (is.character(method) && method == "quantileregression") tune_sigma else NA
    ),
    method_used = method_description
  )
  
  class(result) <- c("run_mvimpute", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_mvimpute <- function(x, ...) {
  cat("=== Missing Value Imputation Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Missing values before: %d\n", x$n_missing_before))
  cat(sprintf("Missing values after:  %d\n", x$n_missing_after))
  invisible(x)
}