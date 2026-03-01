#' Filter Features by Missing Value Threshold
#'
#' @description
#' Removes features (columns) that exceed a specified missing value threshold.
#' Can assess missingness globally or by group, with optional inclusion of QC samples.
#' Supports treating zero values as missing, which is common in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   Missing values should be represented as NA. Zero values can optionally be treated
#'   as missing via the `zero_as_missing` parameter.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain a column specified by `group_col` when `filter_by_group = TRUE`.
#' @param threshold Numeric. Maximum proportion of missing values allowed (0-1).
#'   Features with missing proportion >= threshold are removed. Default: 0.2 (20%).
#' @param filter_by_group Logical. If TRUE, assess missingness within each group separately.
#'   A feature is kept if it passes the threshold in at least one group. Default: TRUE.
#' @param include_QC Logical. If TRUE, include QC samples in missingness calculation.
#'   Default: FALSE (QC samples are excluded from assessment).
#' @param group_col Character. Name of the column in `metadata` containing group information.
#'   Required when `filter_by_group = TRUE`. Default: "Group".
#' @param qc_types Character vector. Group labels that identify QC samples. These samples
#'   are excluded from missingness calculation when `include_QC = FALSE`. Default: c("QC", "SQC", "EQC").
#' @param zero_as_missing Logical. If TRUE, treat zero values as missing (NA) before
#'   applying the filter. This is common in metabolomics where zeros often represent
#'   values below detection limit rather than true zeros. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_filtermissing" containing:
#'   \item{data}{Matrix or data frame of filtered data (same class as input `x`)}
#'   \item{features_removed}{Character vector of removed feature names}
#'   \item{features_kept}{Character vector of retained feature names}
#'   \item{n_features_before}{Integer. Number of features before filtering}
#'   \item{n_features_after}{Integer. Number of features after filtering}
#'   \item{n_features_removed}{Integer. Number of features removed}
#'   \item{n_zeros_converted}{Integer. Number of zero values converted to NA (if applicable)}
#'   \item{parameters}{List of parameters used}
#'   \item{missingness_summary}{Data frame with per-feature missingness statistics}
#'
#' @details
#' This function filters features based on the proportion of missing values (NAs).
#'
#' \strong{Zero values treated as missing}
#'
#' In metabolomics data, zero values often represent one of the following:
#' \itemize{
#'   \item Peak intensities below the detection limit
#'   \item Failed peak integration
#'   \item Metabolites not detected in the sample
#' }
#'
#' When \code{zero_as_missing = TRUE} (default), zero values are converted to
#' missing values before filtering. This is recommended for most untargeted
#' metabolomics datasets.
#'
#' Set \code{zero_as_missing = FALSE} only if one of the following is true:
#' \itemize{
#'   \item Your preprocessing pipeline already converted zeros to NA
#'   \item Zero represents a true biological absence with quantitative meaning
#'   \item Your dataset contains negative values (for example after baseline correction)
#' }
#'
#' \strong{Filtering strategies}
#'
#' Two filtering modes are available:
#'
#' \itemize{
#'   \item \strong{Global filtering} (\code{filter_by_group = FALSE}):
#'   Missingness is calculated across all eligible samples. Features with
#'   missing proportion greater than or equal to \code{threshold} are removed.
#'
#'   \item \strong{Group-wise filtering} (\code{filter_by_group = TRUE}):
#'   Missingness is calculated separately within each group. A feature is
#'   retained if it satisfies the threshold in at least one group. This approach
#'   is less stringent and helps preserve group-specific features.
#' }
#'
#' \strong{QC sample handling}
#'
#' When \code{include_QC = FALSE}, samples whose group labels match
#' \code{qc_types} are excluded from missingness calculations. This is
#' recommended because QC samples often exhibit different missingness
#' patterns compared with biological samples.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control, Batch Correction,
#' Exploration & Data Cleaning (Version 2.1) Zendono. \doi{10.5281/zenodo.16824822}.
#' Retrieved from \url{https://github.com/broadhurstdavid/QC-MXP}.
#' 
#' Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
#' Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data.
#' Scientific Reports, 8(1), 663. \doi{10.1038/s41598-017-19120-0}
#'
#' @importFrom matrixStats colMeans2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with zeros and NAs
#' set.seed(123)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
#' x[sample(length(x), 300)] <- 0    # Add zeros
#' x[sample(length(x), 200)] <- NA   # Add NAs
#' colnames(x) <- paste0("Feature", 1:50)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
#' )
#' 
#' # Treat zeros as missing (default)
#' result1 <- run_filtermissing(x, metadata, threshold = 0.3)
#' 
#' # Keep zeros as valid values
#' result2 <- run_filtermissing(x, metadata, threshold = 0.3, 
#'                               zero_as_missing = FALSE)
#' 
#' # Group-wise filtering with zeros as missing
#' result3 <- run_filtermissing(x, metadata, threshold = 0.3, 
#'                               filter_by_group = TRUE,
#'                               zero_as_missing = TRUE)
#' }
run_filtermissing <- function(
    x,
    metadata,
    threshold = 0.2,
    filter_by_group = TRUE,
    include_QC = FALSE,
    group_col = "Group",
    qc_types = c("QC", "SQC", "EQC"),
    zero_as_missing = TRUE,
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
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0 || threshold > 1) {
    stop("'threshold' must be a numeric value between 0 and 1")
  }
  if (filter_by_group && !group_col %in% colnames(metadata)) {
    stop(sprintf("'%s' column not found in metadata", group_col))
  }
  if (!is.logical(zero_as_missing) || length(zero_as_missing) != 1) {
    stop("'zero_as_missing' must be a single logical value (TRUE or FALSE)")
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  n_features_before <- ncol(x_matrix)
  feature_names <- colnames(x_matrix)
  
  msg(sprintf("Starting missing value filtering (threshold: %.1f%%)...", threshold * 100))
  
  # Convert zeros to NA if requested
  n_zeros_converted <- 0
  if (zero_as_missing) {
    n_zeros_before <- sum(x_matrix == 0, na.rm = TRUE)
    if (n_zeros_before > 0) {
      x_matrix[x_matrix == 0] <- NA
      n_zeros_converted <- n_zeros_before
      msg(sprintf("Converted %d zero values to NA (zero_as_missing = TRUE)", n_zeros_converted))
    } else {
      msg("No zero values found in data (zero_as_missing = TRUE)")
    }
  } else {
    msg("Zeros treated as valid values (zero_as_missing = FALSE)")
  }
  
  # Create Group_ column for QC identification
  if (!include_QC && filter_by_group && group_col %in% colnames(metadata)) {
    metadata$Group_ <- ifelse(metadata[[group_col]] %in% qc_types, "QC", metadata[[group_col]])
  }
  
  # Determine which samples to use for missingness calculation
  if (include_QC || !filter_by_group) {
    assessment_indices <- seq_len(nrow(x_matrix))
  } else {
    qc_indices <- metadata[[group_col]] %in% qc_types
    assessment_indices <- which(!qc_indices)
  }
  
  if (length(assessment_indices) == 0) {
    stop("No samples available for missingness assessment after QC exclusion")
  }
  
  # Calculate missingness
  if (filter_by_group) {
    groups <- unique(metadata[[group_col]][assessment_indices])
    
    missing_by_group <- sapply(groups, function(g) {
      group_indices <- which(metadata[[group_col]] == g & seq_len(nrow(metadata)) %in% assessment_indices)
      if (length(group_indices) == 0) return(rep(0, ncol(x_matrix)))
      colMeans(is.na(x_matrix[group_indices, , drop = FALSE]))
    })
    
    if (is.vector(missing_by_group)) {
      features_to_keep <- missing_by_group < threshold
    } else {
      # Keep feature if it passes threshold in at least one group
      features_to_keep <- !apply(missing_by_group, 1, function(x) all(x >= threshold))
    }
    
    missingness_summary <- as.data.frame(missing_by_group)
    missingness_summary$Feature <- feature_names
    missingness_summary$Kept <- features_to_keep
    
  } else {
    missing_overall <- colMeans(is.na(x_matrix[assessment_indices, , drop = FALSE]))
    features_to_keep <- missing_overall < threshold
    
    missingness_summary <- data.frame(
      Feature = feature_names,
      Missingness = missing_overall,
      Kept = features_to_keep
    )
  }
  
  # Apply filter
  x_filtered <- x_matrix[, features_to_keep, drop = FALSE]
  n_features_after <- ncol(x_filtered)
  n_removed <- n_features_before - n_features_after
  
  msg(sprintf("Removed %d features (%.1f%%) with missingness >= %.1f%%",
              n_removed, (n_removed / n_features_before) * 100, threshold * 100))
  
  # Return in original class
  if (!was_matrix) {
    x_filtered <- as.data.frame(x_filtered)
  }
  
  result <- list(
    data = x_filtered,
    features_removed = feature_names[!features_to_keep],
    features_kept = feature_names[features_to_keep],
    n_features_before = n_features_before,
    n_features_after = n_features_after,
    n_features_removed = n_removed,
    n_zeros_converted = n_zeros_converted,
    parameters = list(
      threshold = threshold,
      filter_by_group = filter_by_group,
      include_QC = include_QC,
      group_col = group_col,
      qc_types = qc_types,
      zero_as_missing = zero_as_missing
    ),
    missingness_summary = missingness_summary
  )
  
  class(result) <- c("run_filtermissing", "list")
  return(result)
}


# Updated S3 print method
#' @export
print.run_filtermissing <- function(x, ...) {
  cat("=== Missing Value Filtering Results ===\n")
  cat(sprintf("Features before: %d\n", x$n_features_before))
  cat(sprintf("Features after:  %d\n", x$n_features_after))
  cat(sprintf("Features removed: %d (%.1f%%)\n", 
              x$n_features_removed, 
              (x$n_features_removed / x$n_features_before) * 100))
  cat(sprintf("Threshold: %.1f%%\n", x$parameters$threshold * 100))
  if (x$parameters$zero_as_missing && x$n_zeros_converted > 0) {
    cat(sprintf("Zeros converted to NA: %d\n", x$n_zeros_converted))
  }
  invisible(x)
}