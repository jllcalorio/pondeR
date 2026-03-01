#' Normalize Metabolomics Data
#'
#' @description
#' Applies various normalization methods to account for dilution effects and
#' sample-to-sample variation in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   May contain normalization factors, group labels, and QC sample identifiers.
#' @param method Character or numeric vector. Normalization method to apply. Options:
#'   \itemize{
#'     \item \code{"sum"}: Total sum normalization
#'     \item \code{"median"}: Median normalization
#'     \item \code{"specific_factor"}: Use values from a metadata column (specify via `factor_col`)
#'     \item \code{"pqn_global"}: Probabilistic Quotient Normalization using global median
#'     \item \code{"pqn_reference"}: PQN using a specific reference sample (specify via `ref_sample`)
#'     \item \code{"pqn_group"}: PQN using pooled QC samples as reference
#'     \item \code{"quantile"}: Quantile normalization
#'     \item Numeric vector: Custom normalization factors (length must equal nrow(x) excluding QCs)
#'   }
#' @param factor_col Character. Name of column in `metadata` containing normalization
#'   factors when `method = "specific_factor"`. Default: "Normalization".
#' @param ref_sample Character. Sample name to use as reference for `method = "pqn_reference"`.
#'   Must match a value in the sample identifier column of metadata. Default: NULL.
#' @param group_sample Character. Group label for `method = "pqn_group"`. Uses samples
#'   with this group label as reference. Default: "QC".
#' @param qc_normalize Character. How to normalize QC samples when using specific factors.
#'   Options: "mean", "median", "none". Default: "median".
#' @param groups Character. Name of column in `metadata` containing group labels.
#'   Required for PQN and QC identification. Default: "Group".
#' @param qc_types Character vector. Group labels identifying QC samples. Default: c("QC", "SQC", "EQC").
#' @param reference_method Character. Method for computing reference in quantile normalization.
#'   Options: "mean", "median". Default: "mean".
#' @param sample_id_col Character. Name of column in `metadata` containing sample identifiers
#'   for matching `ref_sample`. Default: "Sample".
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_normalize" containing:
#'   \item{data}{Matrix or data frame of normalized data (same class as input `x`)}
#'   \item{normalization_factors}{Numeric vector of factors used for normalization}
#'   \item{method_used}{Character describing normalization method applied}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Normalization Methods:**
#' 
#' - **sum**: Divides each sample by its total sum. Assumes similar total metabolite
#'   abundance across samples.
#'   
#' - **median**: Divides each sample by its median value. More robust to outliers than sum.
#' 
#' - **specific_factor**: Uses externally measured factors (e.g., osmolality for urine,
#'   protein concentration for plasma) from a metadata column. QC samples are normalized
#'   separately using mean/median of biological sample factors.
#' 
#' - **pqn_global**: Probabilistic Quotient Normalization using the global median spectrum
#'   as reference. Robust method that accounts for dilution effects.
#' 
#' - **pqn_reference**: PQN using a specific reference sample. Useful when one sample
#'   represents the "ideal" composition.
#' 
#' - **pqn_group**: PQN using pooled QC samples as reference. Recommended for
#'   metabolomics data with regular QC injections.
#' 
#' - **quantile**: Forces all samples to have identical distributions. Very strong
#'   normalization that may remove biological variation - use with caution.
#' 
#' **Custom Normalization Factors:**
#' 
#' You can provide a numeric vector as `method` with custom normalization factors.
#' The vector length must equal the number of biological samples (excluding QCs).
#' QC samples will be normalized using the mean or median of these factors.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006).
#' Probabilistic Quotient Normalization as Robust Method to Account for Dilution
#' of Complex Biological Mixtures. Analytical Chemistry, 78(13), 4281-4290.
#' \doi{10.1021/ac051632c}
#' 
#' Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and signal
#' batch correction for metabolomics datasets. \doi{10.18129/B9.bioc.pmp},
#' R package version 1.20.0, \url{https://bioconductor.org/packages/pmp}.
#'
#' @importFrom matrixStats rowMedians colMedians rowSums2
#' @importFrom pmp pqn_normalisation
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100, sd = 20)),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
#'   Normalization = runif(100, 0.8, 1.2)
#' )
#' 
#' # Sum normalization
#' result1 <- run_normalize(x, metadata, method = "sum")
#' 
#' # PQN with QC reference
#' result2 <- run_normalize(x, metadata, method = "pqn_group")
#' 
#' # Specific factors from metadata
#' result3 <- run_normalize(x, metadata, method = "specific_factor",
#'                          factor_col = "Normalization")
#' }
run_normalize <- function(
    x,
    metadata,
    method = "sum",
    factor_col = "Normalization",
    ref_sample = NULL,
    group_sample = "QC",
    qc_normalize = "median",
    groups = "Group",
    qc_types = c("QC", "SQC", "EQC"),
    reference_method = "mean",
    sample_id_col = "Sample",
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  if (nrow(x) != nrow(metadata)) {
    stop("Number of rows in 'x' must equal number of rows in 'metadata'")
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  
  # Create Group_ column
  if (groups %in% colnames(metadata)) {
    metadata$Group_ <- ifelse(metadata[[groups]] %in% qc_types, "QC", metadata[[groups]])
  } else if (is.character(method) && grepl("pqn", method)) {
    stop(sprintf("'%s' column required in metadata for PQN methods", groups))
  } else {
    metadata$Group_ <- "Sample"  # Assume all samples if no group info
  }
  
  qc_indices <- metadata$Group_ == "QC"
  non_qc_indices <- !qc_indices
  
  msg(sprintf("Applying '%s' normalization...", 
              if (is.numeric(method)) "custom factors" else method))
  
  # Handle custom numeric factors
  if (is.numeric(method)) {
    
    n_non_qc <- sum(non_qc_indices)
    
    if (length(method) != n_non_qc) {
      stop(sprintf("Custom normalization factors must have length %d (number of non-QC samples), got %d",
                   n_non_qc, length(method)))
    }
    
    norm_factors <- method
    
    # Normalize biological samples
    x_matrix[non_qc_indices, ] <- x_matrix[non_qc_indices, ] / norm_factors
    
    # Normalize QC samples
    if (qc_normalize == "mean") {
      qc_factor <- mean(norm_factors, na.rm = TRUE)
    } else if (qc_normalize == "median") {
      qc_factor <- median(norm_factors, na.rm = TRUE)
    } else {
      qc_factor <- 1
    }
    
    x_matrix[qc_indices, ] <- x_matrix[qc_indices, ] / qc_factor
    
    method_description <- "Custom normalization factors"
    all_factors <- rep(NA, nrow(x_matrix))
    all_factors[non_qc_indices] <- norm_factors
    all_factors[qc_indices] <- qc_factor
    
  } else {
    
    method <- tolower(method)
    
    x_matrix <- switch(method,
                       
                       "sum" = {
                         msg("Normalizing by total sum...")
                         row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
                         all_factors <- row_sums
                         x_matrix / row_sums
                       },
                       
                       "median" = {
                         msg("Normalizing by median...")
                         row_medians <- matrixStats::rowMedians(x_matrix, na.rm = TRUE)
                         all_factors <- row_medians
                         x_matrix / row_medians
                       },
                       
                       "specific_factor" = {
                         if (!factor_col %in% colnames(metadata)) {
                           stop(sprintf("Column '%s' not found in metadata", factor_col))
                         }
                         
                         factor_values <- as.numeric(metadata[[factor_col]][non_qc_indices])
                         
                         if (all(is.na(factor_values) | factor_values == 0)) {
                           warning("All normalization factors are NA or 0. Using 'sum' normalization instead.")
                           row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
                           all_factors <- row_sums
                           return(x_matrix / row_sums)
                         }
                         
                         msg(sprintf("Normalizing using factors from '%s' column...", factor_col))
                         
                         # Normalize biological samples
                         x_matrix[non_qc_indices, ] <- x_matrix[non_qc_indices, ] / factor_values
                         
                         # Normalize QC samples
                         if (qc_normalize == "mean") {
                           qc_factor <- mean(factor_values, na.rm = TRUE)
                           msg("QC samples normalized using mean of biological factors")
                         } else if (qc_normalize == "median") {
                           qc_factor <- median(factor_values, na.rm = TRUE)
                           msg("QC samples normalized using median of biological factors")
                         } else {
                           qc_factor <- 1
                           msg("QC samples not normalized")
                         }
                         
                         x_matrix[qc_indices, ] <- x_matrix[qc_indices, ] / qc_factor
                         
                         all_factors <- rep(NA, nrow(x_matrix))
                         all_factors[non_qc_indices] <- factor_values
                         all_factors[qc_indices] <- qc_factor
                         x_matrix
                       },
                       
                       "pqn_global" = {
                         msg("Applying PQN (global median reference)...")
                         reference_spectrum <- matrixStats::colMedians(x_matrix, na.rm = TRUE)
                         quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
                         median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                         all_factors <- median_quotients
                         x_matrix / median_quotients
                       },
                       
                       "pqn_reference" = {
                         if (is.null(ref_sample)) {
                           stop("'ref_sample' must be specified for pqn_reference method")
                         }
                         if (!sample_id_col %in% colnames(metadata)) {
                           stop(sprintf("'%s' column not found in metadata", sample_id_col))
                         }
                         
                         ref_idx <- which(metadata[[sample_id_col]] == ref_sample)
                         if (length(ref_idx) == 0) {
                           stop(sprintf("Reference sample '%s' not found in '%s' column",
                                        ref_sample, sample_id_col))
                         }
                         
                         msg(sprintf("Applying PQN (reference sample: %s)...", ref_sample))
                         reference_spectrum <- as.numeric(x_matrix[ref_idx[1], ])
                         quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
                         median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                         all_factors <- median_quotients
                         x_matrix / median_quotients
                       },
                       
                       "pqn_group" = {
                         pooled_indices <- metadata$Group_ == "QC"
                         
                         if (sum(pooled_indices) == 0) {
                           stop(sprintf("No QC samples found for group PQN. Check '%s' column and qc_types parameter",
                                        groups))
                         }
                         
                         msg(sprintf("Applying PQN (group: %s)...", group_sample))
                         reference_spectrum <- matrixStats::colMedians(
                           x_matrix[pooled_indices, , drop = FALSE], na.rm = TRUE
                         )
                         quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
                         median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                         all_factors <- median_quotients
                         x_matrix / median_quotients
                       },
                       
                       "quantile" = {
                         if (!requireNamespace("pmp", quietly = TRUE)) {
                           stop("Package 'pmp' required for quantile normalization. ",
                                "Install it with: BiocManager::install('pmp')")
                         }
                         
                         tryCatch({
                           msg("Applying quantile normalization...")
                           normalized <- pmp::pqn_normalisation(
                             t(x_matrix),
                             classes = metadata$Group_,
                             qc_label = "QC",
                             ref_mean = NULL,
                             qc_frac = 0,
                             sample_frac = 0,
                             ref_method = reference_method
                           )
                           all_factors <- rep(NA, nrow(x_matrix))  # Not directly available
                           t(normalized)
                         }, error = function(e) {
                           warning("Quantile normalization failed: ", e$message,
                                   ". Using 'sum' normalization instead.")
                           row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
                           all_factors <- row_sums
                           x_matrix / row_sums
                         })
                       },
                       
                       {
                         stop("Unknown normalization method '", method, "'. ",
                              "Supported: 'sum', 'median', 'specific_factor', 'pqn_global', ",
                              "'pqn_reference', 'pqn_group', 'quantile'")
                       }
    )
    
    method_description <- method
  }
  
  # Handle Inf values
  if (any(is.infinite(x_matrix))) {
    msg("Warning: Normalization produced Inf values. Replacing with NA.")
    x_matrix[is.infinite(x_matrix)] <- NA
  }
  
  msg("Normalization complete.")
  
  # Return in original class
  if (!was_matrix) {
    x_matrix <- as.data.frame(x_matrix)
  }
  
  result <- list(
    data = x_matrix,
    normalization_factors = all_factors,
    method_used = method_description,
    parameters = list(
      method = if (is.numeric(method)) "custom" else method,
      factor_col = factor_col,
      ref_sample = ref_sample,
      group_sample = group_sample,
      qc_normalize = qc_normalize,
      groups = groups,
      qc_types = qc_types,
      reference_method = reference_method
    )
  )
  
  class(result) <- c("run_normalize", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_normalize <- function(x, ...) {
  cat("=== Normalization Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Data dimensions: %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  invisible(x)
}