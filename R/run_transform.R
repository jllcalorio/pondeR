#' Transform Metabolomics Data
#'
#' @description
#' Applies various transformation methods to stabilize variance and reduce
#' heteroscedasticity in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param method Character. Transformation method to apply. Options:
#'   \itemize{
#'     \item \code{"log2"}: Log base 2 transformation
#'     \item \code{"log10"}: Log base 10 transformation
#'     \item \code{"sqrt"}: Square root transformation
#'     \item \code{"cbrt"}: Cube root transformation
#'     \item \code{"vsn"}: Variance Stabilizing Normalization (requires 'vsn' package)
#'     \item \code{"glog"}: Generalized logarithm transformation (requires 'pmp' package)
#'   }
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Required for glog transformation to identify QC samples. Default: NULL.
#' @param groups Character. Name of column in `metadata` containing group labels.
#'   Required for glog transformation. Default: "Group".
#' @param qc_types Character vector. Group labels identifying QC samples for glog.
#'   Default: c("QC", "SQC", "EQC").
#' @param num_cores Integer or "max". Number of CPU cores to use for VSN transformation.
#'   Use "max" for automatic detection (uses all available - 2). Default: 1.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_transform" containing:
#'   \item{data}{Matrix or data frame of transformed data (same class as input `x`)}
#'   \item{method_used}{Character describing transformation method applied}
#'   \item{shift_applied}{Numeric value added to data before transformation (if applicable)}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Transformation Methods:**
#' 
#' - **log2 / log10**: Logarithmic transformations compress large values and expand
#'   small values. Useful for data spanning multiple orders of magnitude. If data
#'   contains zeros or negative values, a shift is applied first.
#' 
#' - **sqrt**: Square root transformation. Moderate variance stabilization. Works
#'   with zero values but requires non-negative data.
#' 
#' - **cbrt**: Cube root transformation. Similar to sqrt but can handle negative
#'   values. Provides gentler compression than log.
#' 
#' - **vsn**: Variance Stabilizing Normalization. Model-based method from microarray
#'   analysis that calibrates variance across intensity range. Requires 'vsn' package.
#'   Can be computationally intensive; use `num_cores` for parallel processing.
#' 
#' - **glog**: Generalized logarithm transformation. Hybrid log transformation that
#'   handles low values better than standard log. Requires 'pmp' package and metadata
#'   with QC sample labels.
#' 
#' **When to Use Each Method:**
#' 
#' - **log2/log10**: Most common for metabolomics. Use when variance increases
#'   with mean (heteroscedasticity).
#' 
#' - **sqrt**: When variance is proportional to mean.
#' 
#' - **vsn**: When you need sophisticated variance modeling. Recommended for
#'   large datasets with clear intensity-dependent variance.
#' 
#' - **glog**: Alternative to log when you have many low-intensity values.
#'   Good for LC-MS data with detection limits.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Huber, W., von Heydebreck, A., Sueltmann, H., Poustka, A., & Vingron, M. (2002).
#' Variance stabilization applied to microarray data calibration and to the
#' quantification of differential expression. Bioinformatics, 18(Suppl 1), S96-S104.
#' \doi{10.1093/bioinformatics/18.suppl_1.s96}
#' 
#' Parsons, H.M., Ludwig, C., Günther, U.L., & Viant, M.R. (2007).
#' Improved classification accuracy in 1- and 2-dimensional NMR metabolomics data
#' using the variance stabilising generalised logarithm transformation.
#' BMC Bioinformatics, 8, 234. \doi{10.1186/1471-2105-8-234}
#'
#' @importFrom vsn vsn2 predict
#' @importFrom pmp glog_transformation
#' @importFrom BiocParallel register SnowParam SerialParam
#' @importFrom parallelly availableCores
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100, sd = 50)),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' 
#' # Log transformation
#' result1 <- run_transform(x, method = "log2")
#' 
#' # VSN with parallel processing
#' result2 <- run_transform(x, method = "vsn", num_cores = 4)
#' 
#' # Glog transformation (requires metadata)
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
#' )
#' result3 <- run_transform(x, method = "glog", metadata = metadata)
#' }
run_transform <- function(
    x,
    method = "log2",
    metadata = NULL,
    groups = "Group",
    qc_types = c("QC", "SQC", "EQC"),
    num_cores = 1,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  
  method <- tolower(method)
  valid_methods <- c("log2", "log10", "sqrt", "cbrt", "vsn", "glog")
  if (!method %in% valid_methods) {
    stop("Unknown transformation method '", method, "'. ",
         "Supported: ", paste(valid_methods, collapse = ", "))
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  shift_value <- 0
  
  msg(sprintf("Applying '%s' transformation...", method))
  
  # Handle num_cores
  if (method == "vsn") {
    max_available_cores <- parallelly::availableCores(omit = 2)
    
    if (num_cores == "max") {
      num_cores <- max_available_cores
    } else if (is.numeric(num_cores)) {
      num_cores <- as.integer(num_cores)
      if (num_cores < 1 || num_cores > max_available_cores) {
        stop(sprintf("'num_cores' must be between 1 and %d (available cores)", 
                     max_available_cores))
      }
    } else {
      stop("'num_cores' must either be 'max' or a positive integer")
    }
  }
  
  x_transformed <- switch(method,
                          
                          "log2" = {
                            min_val <- min(x_matrix, na.rm = TRUE)
                            if (min_val <= 0) {
                              shift_value <- abs(min_val) + 1
                              x_matrix <- x_matrix + shift_value
                              msg(sprintf("Added shift of %.2f to handle non-positive values", shift_value))
                            }
                            msg("Applying log2 transformation...")
                            log2(x_matrix)
                          },
                          
                          "log10" = {
                            min_val <- min(x_matrix, na.rm = TRUE)
                            if (min_val <= 0) {
                              shift_value <- abs(min_val) + 1
                              x_matrix <- x_matrix + shift_value
                              msg(sprintf("Added shift of %.2f to handle non-positive values", shift_value))
                            }
                            msg("Applying log10 transformation...")
                            log10(x_matrix)
                          },
                          
                          "sqrt" = {
                            min_val <- min(x_matrix, na.rm = TRUE)
                            if (min_val < 0) {
                              shift_value <- abs(min_val)
                              x_matrix <- x_matrix + shift_value
                              msg(sprintf("Added shift of %.2f to handle negative values", shift_value))
                            }
                            msg("Applying square root transformation...")
                            sqrt(x_matrix)
                          },
                          
                          "cbrt" = {
                            min_val <- min(x_matrix, na.rm = TRUE)
                            if (min_val < 0) {
                              shift_value <- abs(min_val) + 1
                              x_matrix <- x_matrix + shift_value
                              msg(sprintf("Added shift of %.2f to handle negative values", shift_value))
                            }
                            msg("Applying cube root transformation...")
                            x_matrix^(1/3)
                          },
                          
                          "vsn" = {
                            if (!requireNamespace("vsn", quietly = TRUE)) {
                              stop("Package 'vsn' required for VSN transformation. ",
                                   "Install it with: BiocManager::install('vsn')")
                            }
                            
                            if (!requireNamespace("BiocParallel", quietly = TRUE)) {
                              stop("Package 'BiocParallel' required for VSN. ",
                                   "Install it with: BiocManager::install('BiocParallel')")
                            }
                            
                            tryCatch({
                              # Setup parallel processing
                              BPPARAM_to_use <- BiocParallel::SnowParam(workers = num_cores)
                              BiocParallel::register(BPPARAM_to_use)
                              
                              msg(sprintf("Applying VSN transformation (using %d cores)...", num_cores))
                              vsn_fit <- vsn::vsn2(x_matrix)
                              result <- vsn::predict(vsn_fit, newdata = x_matrix)
                              
                              # Reset to serial processing
                              BiocParallel::register(BiocParallel::SerialParam())
                              
                              result
                              
                            }, error = function(e) {
                              # Always reset to serial on error
                              BiocParallel::register(BiocParallel::SerialParam())
                              
                              warning("VSN transformation failed: ", e$message,
                                      ". Using log10 instead.")
                              min_val <- min(x_matrix, na.rm = TRUE)
                              if (min_val <= 0) {
                                shift_value <<- abs(min_val) + 1
                                x_matrix <<- x_matrix + shift_value
                              }
                              log10(x_matrix)
                            })
                          },
                          
                          "glog" = {
                            if (!requireNamespace("pmp", quietly = TRUE)) {
                              stop("Package 'pmp' required for glog transformation. ",
                                   "Install it with: BiocManager::install('pmp')")
                            }
                            
                            if (is.null(metadata)) {
                              stop("'metadata' required for glog transformation (to identify QC samples)")
                            }
                            if (nrow(metadata) != nrow(x_matrix)) {
                              stop("'metadata' must have same number of rows as 'x'")
                            }
                            if (!groups %in% colnames(metadata)) {
                              stop(sprintf("'%s' column not found in metadata", groups))
                            }
                            
                            # Create Group_ column
                            metadata$Group_ <- ifelse(metadata[[groups]] %in% qc_types,
                                                      "QC", metadata[[groups]])
                            
                            tryCatch({
                              msg("Applying glog transformation...")
                              glog_result <- pmp::glog_transformation(
                                t(x_matrix),
                                classes = metadata$Group_,
                                qc_label = "QC"
                              )
                              t(glog_result)
                            }, error = function(e) {
                              warning("glog transformation failed: ", e$message,
                                      ". Using log10 instead.")
                              min_val <- min(x_matrix, na.rm = TRUE)
                              if (min_val <= 0) {
                                shift_value <<- abs(min_val) + 1
                                x_matrix <<- x_matrix + shift_value
                              }
                              log10(x_matrix)
                            })
                          }
  )
  
  msg("Transformation complete.")
  
  # Return in original class
  if (!was_matrix) {
    x_transformed <- as.data.frame(x_transformed)
  }
  
  result <- list(
    data = x_transformed,
    method_used = method,
    shift_applied = shift_value,
    parameters = list(
      method = method,
      num_cores = if (method == "vsn") num_cores else NA,
      groups = groups,
      qc_types = qc_types
    )
  )
  
  class(result) <- c("run_transform", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_transform <- function(x, ...) {
  cat("=== Transformation Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  if (x$shift_applied > 0) {
    cat(sprintf("Shift applied: %.4f\n", x$shift_applied))
  }
  cat(sprintf("Data dimensions: %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  invisible(x)
}