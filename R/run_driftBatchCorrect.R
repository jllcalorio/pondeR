#' Correct Signal Drift and Batch Effects Using QC Samples
#'
#' @description
#' Applies Quality Control-based Robust Spline Correction (QCRSC) or batch correction
#' using ComBat to remove systematic signal drift and batch effects in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain columns specified by `injection_sequence`, `batch_numbers`, and `groups`.
#' @param perform_correction Logical. If TRUE, perform correction; if FALSE, return
#'   original data unchanged. Default: TRUE.
#' @param batch_corr_only Logical. If TRUE, perform only batch correction using ComBat
#'   (no drift correction). If FALSE, perform QCRSC (drift + batch correction). Default: FALSE.
#' @param injection_sequence Character. Name of column in `metadata` containing injection
#'   order. Required for QCRSC. Default: "InjectionSequence".
#' @param batch_numbers Character. Name of column in `metadata` containing batch numbers.
#'   Required for both QCRSC and ComBat. Default: "Batches".
#' @param groups Character. Name of column in `metadata` containing sample group labels
#'   (e.g., "Sample", "QC", "SQC", "EQC"). Required for QCRSC. Default: "Group".
#' @param qc_label Character. Group label identifying QC samples for QCRSC. Default: "QC".
#' @param qc_types Character vector. Group labels that should be converted to `qc_label`
#'   internally. Default: c("QC", "SQC", "EQC").
#' @param spline_smooth_param Numeric. Smoothing parameter for spline fitting in QCRSC
#'   (0 = more flexible, 1 = more rigid). Default: 0.
#' @param min_QC Integer. Minimum number of QC samples required per batch for QCRSC.
#'   Batches with fewer QCs will not be corrected. Default: 5.
#' @param spar_limit Numeric vector of length 2. Limits for spline fitting in QCRSC
#'   as c(min, max). Default: c(-1.5, 1.5).
#' @param log_scale Logical. If TRUE, fit QCRSC spline on log-transformed data. Default: TRUE.
#' @param use_parametric Logical. For ComBat only. Use parametric adjustment for batch
#'   effects. Default: TRUE.
#' @param display_plots Logical. For ComBat only. Display diagnostic plots during correction.
#'   Default: FALSE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_driftBatchCorrect" containing:
#'   \item{data}{Matrix or data frame of corrected data (same class as input `x`)}
#'   \item{data_before_correction}{Original data before correction}
#'   \item{correction_applied}{Logical indicating if correction was performed}
#'   \item{method_used}{Character describing correction method ("QCRSC" or "ComBat")}
#'   \item{uncorrected_features}{Character vector of features that could not be corrected}
#'   \item{n_uncorrected}{Integer. Number of uncorrected features}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **QCRSC (Quality Control-based Robust Spline Correction):**
#' 
#' QCRSC is the default and recommended method for metabolomics data with regular
#' QC injections. It works by:
#' 
#' 1. Fitting a smoothing spline through QC sample intensities across injection order
#' 2. Using this spline to model systematic drift and batch effects
#' 3. Normalizing all samples (biological + QC) to remove the fitted drift
#' 
#' **When to use QCRSC:**
#' - Long analytical runs with visible signal drift
#' - Multiple batches with sufficient QC samples (≥ `min_QC` per batch)
#' - QC samples show systematic trends across injection order
#' 
#' **ComBat Batch Correction:**
#' 
#' When `batch_corr_only = TRUE`, uses the ComBat algorithm for batch effect removal
#' only (no drift correction). This is useful when:
#' 
#' - QC samples are insufficient or absent
#' - Primary concern is batch differences, not within-batch drift
#' - Data comes from multiple platforms or sites
#' 
#' ComBat uses empirical Bayes methods to adjust for batch effects while preserving
#' biological variation.
#' 
#' **Uncorrected Features:**
#' 
#' Some features may remain uncorrected due to:
#' - Insufficient QC samples (< `min_QC`) in one or more batches
#' - Poor spline fit or numerical issues
#' - All QC values identical (no variance to model)
#' 
#' Uncorrected features are returned unchanged with their names listed in the output.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Kirwan, J.A., Broadhurst, D.I., Davidson, R.L. et al. (2013).
#' Characterising and correcting batch variation in an automated direct infusion
#' mass spectrometry (DIMS) metabolomics workflow.
#' Analytical and Bioanalytical Chemistry, 405, 5147-5157. \doi{10.1007/s00216-013-6856-7}
#' 
#' Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray
#' expression data using empirical Bayes methods. Biostatistics, 8(1), 118-127.
#' \doi{10.1093/biostatistics/kxj037}
#' 
#' Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and signal
#' batch correction for metabolomics datasets. \doi{10.18129/B9.bioc.pmp},
#' R package version 1.20.0, \url{https://bioconductor.org/packages/pmp}.
#'
#' @importFrom pmp QCRSC
#' @importFrom sva ComBat
#' @importFrom matrixStats colMaxs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with drift and batch effects
#' set.seed(123)
#' n_samples <- 120
#' n_features <- 50
#' x <- matrix(rnorm(n_samples * n_features, mean = 100, sd = 20),
#'             nrow = n_samples, ncol = n_features)
#' colnames(x) <- paste0("Feature", 1:n_features)
#' 
#' # Add systematic drift
#' drift <- seq(0, 30, length.out = n_samples)
#' x <- x + drift
#' 
#' metadata <- data.frame(
#'   Sample = paste0("S", 1:n_samples),
#'   InjectionSequence = 1:n_samples,
#'   Batches = rep(1:3, each = 40),
#'   Group = rep(c(rep("Sample", 9), "QC"), n_samples / 10)
#' )
#' 
#' # Apply QCRSC correction
#' result1 <- run_driftBatchCorrect(x, metadata)
#' 
#' # Apply ComBat correction only
#' result2 <- run_driftBatchCorrect(x, metadata, batch_corr_only = TRUE)
#' }
run_driftBatchCorrect <- function(
    x,
    metadata,
    perform_correction = TRUE,
    batch_corr_only = FALSE,
    injection_sequence = "InjectionSequence",
    batch_numbers = "Batches",
    groups = "Group",
    qc_label = "QC",
    qc_types = c("QC", "SQC", "EQC"),
    spline_smooth_param = 0,
    min_QC = 5,
    spar_limit = c(-1.5, 1.5),
    log_scale = TRUE,
    use_parametric = TRUE,
    display_plots = FALSE,
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
  
  if (!perform_correction) {
    msg("Correction disabled. Returning original data.")
    result <- list(
      data = x,
      data_before_correction = x,
      correction_applied = FALSE,
      method_used = "None",
      uncorrected_features = character(0),
      n_uncorrected = 0,
      parameters = list(perform_correction = FALSE)
    )
    class(result) <- c("run_driftBatchCorrect", "list")
    return(result)
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  x_before <- x_matrix
  
  # Create Group_ column
  if (groups %in% colnames(metadata)) {
    metadata$Group_ <- ifelse(metadata[[groups]] %in% qc_types, qc_label, metadata[[groups]])
  } else {
    stop(sprintf("'%s' column not found in metadata", groups))
  }
  
  uncorrected_features <- character(0)
  method_used <- NULL
  
  if (!batch_corr_only) {
    
    # QCRSC correction
    msg("Applying QCRSC drift and batch correction...")
    
    # Validate required columns
    required_cols <- c(injection_sequence, batch_numbers, groups)
    missing_cols <- setdiff(required_cols, colnames(metadata))
    if (length(missing_cols) > 0) {
      stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
    }
    
    if (!requireNamespace("pmp", quietly = TRUE)) {
      stop("Package 'pmp' is required for QCRSC. Install it with: BiocManager::install('pmp')")
    }
    
    # Validate parameters
    if (!is.numeric(spline_smooth_param) || length(spline_smooth_param) != 1 ||
        spline_smooth_param < 0 || spline_smooth_param > 1) {
      stop("'spline_smooth_param' must be numeric between 0 and 1")
    }
    if (!is.numeric(spar_limit) || length(spar_limit) != 2) {
      stop("'spar_limit' must be a numeric vector of length 2")
    }
    
    tryCatch({
      
      x_corrected <- pmp::QCRSC(
        df = t(x_matrix),  # QCRSC expects features in rows
        order = as.numeric(metadata[[injection_sequence]]),
        batch = as.numeric(metadata[[batch_numbers]]),
        classes = as.vector(metadata$Group_),
        spar = spline_smooth_param,
        log = log_scale,
        minQC = min_QC,
        qc_label = qc_label,
        spar_lim = spar_limit
      )
      
      x_matrix <- t(x_corrected)
      method_used <- "QCRSC"
      
      # Identify uncorrected features
      tryCatch({
        if (identical(dim(x_before), dim(x_matrix)) &&
            identical(colnames(x_before), colnames(x_matrix))) {
          
          tolerance <- .Machine$double.eps^0.5
          diff_matrix <- abs(x_before - x_matrix)
          col_max_diffs <- matrixStats::colMaxs(diff_matrix, na.rm = TRUE)
          
          uncorrected_indices <- (col_max_diffs <= tolerance)
          uncorrected_features <- colnames(x_before)[uncorrected_indices]
          
          if (length(uncorrected_features) > 0) {
            msg(sprintf("Warning: %d features could not be corrected (insufficient QC samples)",
                        length(uncorrected_features)))
          }
        }
      }, error = function(e) {
        warning("Could not identify uncorrected features: ", e$message)
      })
      
      msg("QCRSC correction completed successfully.")
      
    }, error = function(e) {
      stop("QCRSC correction failed: ", e$message,
           "\nConsider using batch_corr_only = TRUE or checking your QC samples.")
    })
    
  } else {
    
    # ComBat batch correction only
    msg("Applying ComBat batch correction...")
    
    if (!batch_numbers %in% colnames(metadata)) {
      stop(sprintf("'%s' column not found in metadata", batch_numbers))
    }
    
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("Package 'sva' is required for ComBat. Install it with: BiocManager::install('sva')")
    }
    
    batch_vector <- as.numeric(metadata[[batch_numbers]])
    
    if (length(unique(batch_vector)) < 2) {
      warning("Only one batch detected. Batch correction not applicable. Returning original data.")
      x_matrix <- x_before
      method_used <- "None (single batch)"
    } else {
      
      tryCatch({
        
        x_corrected <- sva::ComBat(
          dat = t(x_matrix),  # ComBat expects features in rows
          batch = batch_vector,
          mod = NULL,
          par.prior = use_parametric,
          prior.plots = display_plots,
          mean.only = FALSE,
          ref.batch = NULL
        )
        
        x_matrix <- t(x_corrected)
        method_used <- "ComBat"
        
        msg("ComBat correction completed successfully.")
        
      }, error = function(e) {
        stop("ComBat correction failed: ", e$message)
      })
    }
  }
  
  # Handle any new NAs introduced by correction
  new_na_count <- sum(is.na(x_matrix)) - sum(is.na(x_before))
  if (new_na_count > 0) {
    msg(sprintf("Warning: Correction introduced %d new missing values", new_na_count))
  }
  
  # Return in original class
  if (!was_matrix) {
    x_matrix <- as.data.frame(x_matrix)
    x_before <- as.data.frame(x_before)
  }
  
  result <- list(
    data = x_matrix,
    data_before_correction = x_before,
    correction_applied = TRUE,
    method_used = method_used,
    uncorrected_features = uncorrected_features,
    n_uncorrected = length(uncorrected_features),
    parameters = list(
      perform_correction = perform_correction,
      batch_corr_only = batch_corr_only,
      injection_sequence = injection_sequence,
      batch_numbers = batch_numbers,
      groups = groups,
      qc_label = qc_label,
      qc_types = qc_types,
      spline_smooth_param = spline_smooth_param,
      min_QC = min_QC,
      spar_limit = spar_limit,
      log_scale = log_scale,
      use_parametric = use_parametric,
      display_plots = display_plots
    )
  )
  
  class(result) <- c("run_driftBatchCorrect", "list")
  return(result)
}

# S3 print methods for cleaner output
#' @export
print.run_driftBatchCorrect <- function(x, ...) {
  cat("=== Drift/Batch Correction Results ===\n")
  cat(sprintf("Correction applied: %s\n", x$correction_applied))
  cat(sprintf("Method used: %s\n", x$method_used))
  if (x$n_uncorrected > 0) {
    cat(sprintf("Uncorrected features: %d\n", x$n_uncorrected))
  }
  invisible(x)
}