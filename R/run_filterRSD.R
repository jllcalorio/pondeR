#' Filter Features by Relative Standard Deviation in QC Samples
#'
#' @description
#' Removes features with high relative standard deviation (RSD, also known as
#' coefficient of variation) in quality control samples. High RSD indicates
#' poor analytical reproducibility.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain column specified by `groups` to identify QC samples.
#' @param max_rsd Numeric. Maximum allowed RSD (0-1). Features with RSD >= this
#'   threshold in QC samples are removed. Default: 0.3 (30%).
#' @param groups Character. Name of column in `metadata` containing sample group labels.
#'   Default: "Group".
#' @param qc_type Character. QC sample type to use for RSD calculation. Options:
#'   \itemize{
#'     \item \code{"EQC"}: Extract QC samples (diluted sample pool)
#'     \item \code{"SQC"}: Sample QC samples (undiluted sample pool)
#'     \item \code{"QC"}: Any QC sample (use when both EQC and SQC present, or no distinction)
#'   }
#'   Default: "EQC".
#' @param qc_types Character vector. Group labels that should be treated as QC samples.
#'   Default: c("QC", "SQC", "EQC").
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_filterRSD" containing:
#'   \item{data}{Matrix or data frame of filtered data (same class as input `x`)}
#'   \item{features_removed}{Character vector of removed feature names}
#'   \item{features_kept}{Character vector of retained feature names}
#'   \item{rsd_values}{Numeric vector of RSD values for all original features}
#'   \item{n_features_before}{Integer. Number of features before filtering}
#'   \item{n_features_after}{Integer. Number of features after filtering}
#'   \item{n_features_removed}{Integer. Number of features removed}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Relative Standard Deviation (RSD):**
#' 
#' RSD is calculated as: RSD = (SD / Mean) × 100%
#' 
#' Also known as coefficient of variation (CV), RSD expresses variability
#' as a percentage of the mean, making it scale-independent and suitable
#' for comparing features with different magnitudes.
#' 
#' **Why Filter by RSD:**
#' 
#' Quality control samples should have low variability since they represent
#' the same biological matrix. High RSD in QC samples indicates:
#' - Poor analytical reproducibility
#' - Instrument instability
#' - Ion suppression/enhancement effects
#' - Features near detection limit
#' 
#' Features with high QC RSD are unreliable and should be removed before
#' biological interpretation.
#' 
#' **Choosing QC Type:**
#' 
#' - **EQC** (default): Most commonly used. Represents diluted pooled samples
#'   analyzed regularly throughout the run.
#'   
#' - **SQC**: Undiluted pooled samples. Use when you have a specific SQC type
#'   distinct from EQC.
#'   
#' - **QC**: Use when you don't distinguish between QC types, or want to
#'   assess RSD across all QC samples combined.
#' 
#' **Typical RSD Thresholds:**
#' 
#' - Stringent: 20% (0.2) - high-quality data only
#' - Standard: 30% (0.3) - recommended for most metabolomics studies
#' - Lenient: 40% (0.4) - when sample size is limited
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control, Batch Correction,
#' Exploration & Data Cleaning (Version 2.1) Zendono. \doi{10.5281/zenodo.16824822}.
#' Retrieved from \url{https://github.com/broadhurstdavid/QC-MXP}.
#' 
#' Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and signal
#' batch correction for metabolomics datasets. \doi{10.18129/B9.bioc.pmp},
#' R package version 1.20.0, \url{https://bioconductor.org/packages/pmp}.
#'
#' @importFrom pmp filter_peaks_by_rsd
#' @importFrom matrixStats colSds colMeans2
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
#' # Add high variability to some features in QC samples
#' qc_rows <- 91:100
#' x[qc_rows, 1:10] <- x[qc_rows, 1:10] * runif(10 * 10, 0.5, 1.5)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:100),
#'   Group = c(rep("Control", 40), rep("Treatment", 40), rep("QC", 10), rep("EQC", 10))
#' )
#' 
#' # Filter by RSD in EQC samples
#' result <- run_filterRSD(x, metadata, max_rsd = 0.3, qc_type = "EQC")
#' }
run_filterRSD <- function(
    x,
    metadata,
    max_rsd = 0.3,
    groups = "Group",
    qc_type = "EQC",
    qc_types = c("QC", "SQC", "EQC"),
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
  if (!is.numeric(max_rsd) || length(max_rsd) != 1 || max_rsd < 0 || max_rsd > 1) {
    stop("'max_rsd' must be a numeric value between 0 and 1")
  }
  if (!groups %in% colnames(metadata)) {
    stop(sprintf("'%s' column not found in metadata", groups))
  }
  
  valid_qc_types <- c("SQC", "EQC", "QC")
  if (!qc_type %in% valid_qc_types) {
    stop("'qc_type' must be one of: ", paste(valid_qc_types, collapse = ", "))
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  feature_names <- colnames(x_matrix)
  n_features_before <- ncol(x_matrix)
  
  msg(sprintf("Filtering features by RSD (threshold: %.1f%%, QC type: %s)...",
              max_rsd * 100, qc_type))
  
  # Create Group_ column for QC identification
  if (qc_type == "QC") {
    metadata$Group_ <- ifelse(metadata[[groups]] %in% qc_types, "QC", metadata[[groups]])
    qc_label <- "QC"
    classes <- metadata$Group_
  } else {
    qc_label <- qc_type
    classes <- metadata[[groups]]
  }
  
  # Check if QC samples exist
  if (!qc_label %in% classes) {
    warning(sprintf("No %s samples found. Returning original data without filtering.", qc_label))
    result <- list(
      data = if (was_matrix) x_matrix else as.data.frame(x_matrix),
      features_removed = character(0),
      features_kept = feature_names,
      rsd_values = rep(NA, n_features_before),
      n_features_before = n_features_before,
      n_features_after = n_features_before,
      n_features_removed = 0,
      parameters = list(max_rsd = max_rsd, groups = groups, qc_type = qc_type)
    )
    class(result) <- c("run_filterRSD", "list")
    return(result)
  }
  
  # Calculate RSD in QC samples manually for reporting
  qc_indices <- which(classes == qc_label)
  x_qc <- x_matrix[qc_indices, , drop = FALSE]
  
  qc_means <- matrixStats::colMeans2(x_qc, na.rm = TRUE)
  qc_sds <- matrixStats::colSds(x_qc, na.rm = TRUE)
  rsd_values <- (qc_sds / qc_means)
  rsd_values[!is.finite(rsd_values)] <- NA
  
  # Use pmp for filtering if available
  if (requireNamespace("pmp", quietly = TRUE)) {
    
    tryCatch({
      
      x_filtered <- pmp::filter_peaks_by_rsd(
        t(x_matrix),  # pmp expects features in rows
        max_rsd = max_rsd * 100,  # pmp uses percentage (0-100)
        class = as.vector(classes),
        qc_label = qc_label
      )
      
      x_filtered <- t(x_filtered)
      features_kept <- colnames(x_filtered)
      
    }, error = function(e) {
      warning("pmp::filter_peaks_by_rsd failed: ", e$message,
              ". Using manual RSD filtering.")
      
      # Manual fallback
      features_kept <<- feature_names[rsd_values < max_rsd | is.na(rsd_values)]
      x_filtered <<- x_matrix[, features_kept, drop = FALSE]
    })
    
  } else {
    # Manual filtering when pmp not available
    msg("Package 'pmp' not available. Using manual RSD filtering.")
    features_kept <- feature_names[rsd_values < max_rsd | is.na(rsd_values)]
    x_filtered <- x_matrix[, features_kept, drop = FALSE]
  }
  
  n_features_after <- length(features_kept)
  n_removed <- n_features_before - n_features_after
  features_removed <- setdiff(feature_names, features_kept)
  
  msg(sprintf("Removed %d features (%.1f%%) with RSD >= %.1f%%",
              n_removed, (n_removed / n_features_before) * 100, max_rsd * 100))
  
  # Return in original class
  if (!was_matrix) {
    x_filtered <- as.data.frame(x_filtered)
  }
  
  result <- list(
    data = x_filtered,
    features_removed = features_removed,
    features_kept = features_kept,
    rsd_values = rsd_values,
    n_features_before = n_features_before,
    n_features_after = n_features_after,
    n_features_removed = n_removed,
    parameters = list(
      max_rsd = max_rsd,
      groups = groups,
      qc_type = qc_type,
      qc_types = qc_types,
      qc_label = qc_label
    )
  )
  
  class(result) <- c("run_filterRSD", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_filterRSD <- function(x, ...) {
  cat("=== RSD Filtering Results ===\n")
  cat(sprintf("Features before: %d\n", x$n_features_before))
  cat(sprintf("Features after:  %d\n", x$n_features_after))
  cat(sprintf("Features removed: %d (%.1f%%)\n", 
              x$n_features_removed,
              (x$n_features_removed / x$n_features_before) * 100))
  cat(sprintf("RSD threshold: %.1f%%\n", x$parameters$max_rsd * 100))
  cat(sprintf("QC type used: %s\n", x$parameters$qc_label))
  invisible(x)
}