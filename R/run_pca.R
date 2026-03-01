#' Perform Principal Component Analysis
#'
#' @description
#' Performs Principal Component Analysis (PCA) on preprocessed metabolomics data.
#' Optionally applies scaling and transformation before PCA if the input data
#' has not been pre-scaled. Supports sample exclusion by group.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   Can be output from `run_scale()` or raw data (which will be scaled if `scale_method`
#'   is specified).
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain columns for grouping, batch information, and sample identifiers.
#' @param scale_method Character or NULL. Scaling method to apply before PCA if `x` is not
#'   already scaled. Options: "mean", "auto", "pareto". If NULL (default) and `x` is not
#'   from `run_scale()`, a warning is issued and PCA proceeds without additional scaling.
#' @param transform_method Character or NULL. Transformation method to apply before scaling
#'   and PCA. Options: "log2", "log10", "sqrt", "cbrt", "vsn", "glog". Default: NULL (no transformation).
#' @param group Character. Name of column in `metadata` containing group labels for sample
#'   exclusion. Required if `exclude` is specified. Default: "Group".
#' @param exclude Character vector or NULL. Group labels (from the `group` column) to exclude
#'   from PCA analysis. For example, c("QC", "Blank") excludes all samples with these group
#'   labels. Default: NULL (include all samples).
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class "run_pca" containing:
#'   \item{scores}{Matrix of PC scores (samples × PCs)}
#'   \item{loadings}{Matrix of PC loadings (features × PCs)}
#'   \item{variance_explained}{Numeric vector of variance explained (%) per PC}
#'   \item{eigenvalues}{Numeric vector of eigenvalues per PC}
#'   \item{cumulative_variance}{Numeric vector of cumulative variance explained (%) per PC}
#'   \item{center_values}{Numeric vector of centering values used (if any)}
#'   \item{scale_values}{Numeric vector of scaling values used (if any)}
#'   \item{pca_object}{Original prcomp object from stats::prcomp}
#'   \item{data_used}{Matrix of data used for PCA (after any transformation/scaling)}
#'   \item{metadata}{Data frame of metadata aligned with PCA scores}
#'   \item{n_samples}{Integer. Number of samples used in PCA}
#'   \item{n_samples_excluded}{Integer. Number of samples excluded}
#'   \item{n_features}{Integer. Number of features}
#'   \item{n_pcs}{Integer. Number of PCs computed}
#'   \item{excluded_groups}{Character vector of excluded group labels}
#'   \item{preprocessing}{List documenting preprocessing steps applied}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Sample Exclusion:**
#' 
#' The `exclude` parameter allows removal of specific sample groups before PCA:
#' - Common use: Exclude QC samples (`exclude = c("QC", "SQC", "EQC")`)
#' - Multiple groups can be excluded: `exclude = c("QC", "Blank", "Failed")`
#' - Exclusion is performed before any preprocessing steps
#' - Excluded samples are not included in transformation/scaling calculations
#' 
#' **Input Data Requirements:**
#' 
#' The function accepts data in two forms:
#' 
#' 1. **Pre-scaled data** from `run_scale()`: PCA is performed directly without
#'    additional preprocessing. This is the recommended workflow.
#'    
#' 2. **Raw data**: If `scale_method` is specified, the function will first apply
#'    optional transformation (via `transform_method`), then scaling (via `scale_method`),
#'    before performing PCA.
#' 
#' **PCA Methodology:**
#' 
#' PCA is performed using singular value decomposition (SVD) via `stats::prcomp()`.
#' The function assumes data is appropriately preprocessed (normalized, transformed,
#' and scaled) before analysis. PCA is performed without additional centering or
#' scaling (`center = FALSE, scale. = FALSE`) as these should have been applied
#' during the scaling step.
#' 
#' **Interpretation:**
#' 
#' - **PC Scores**: Coordinates of samples in the principal component space.
#'   Used for visualization and identifying patterns/clusters.
#'   
#' - **PC Loadings**: Contribution of each feature to each principal component.
#'   Features with high absolute loadings drive separation along that PC.
#'   
#' - **Variance Explained**: Percentage of total variance captured by each PC.
#'   First few PCs typically capture most variance in metabolomics data.
#'   
#' - **Eigenvalues**: Variance captured by each PC in original units.
#'   Eigenvalue > 1 (Kaiser criterion) often used to determine significant PCs.
#' 
#' **Scaling Recommendations:**
#' 
#' - **Auto-scaling** ("auto"): Recommended for PCA. Gives equal weight to all
#'   features regardless of magnitude.
#'   
#' - **Pareto-scaling** ("pareto"): Alternative that partially preserves variance
#'   structure. Can be useful when feature magnitudes carry biological meaning.
#'   
#' - **Mean-centering** ("mean"): Only centers data. Use when features are already
#'   on comparable scales.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). The New S Language.
#' Wadsworth & Brooks/Cole. \doi{10.1201/9781351074988}
#' 
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). Multivariate Analysis.
#' London: Academic Press.
#' 
#' Venables, W.N., & Ripley, B.D. (2002). Modern Applied Statistics with S.
#' Springer-Verlag.
#' 
#' van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., & van der Werf, M.J. (2006).
#' Centering, scaling, and transformations: improving the biological information content of
#' metabolomics data. BMC Genomics, 7, 142. \doi{10.1186/1471-2164-7-142}
#'
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 50, mean = 100, sd = 20),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' rownames(x) <- paste0("Sample", 1:100)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("Sample", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
#'   Batch = rep(1:4, each = 25),
#'   Time = rep(c(0, 6, 12, 24), 25)
#' )
#' 
#' # PCA on pre-scaled data (recommended)
#' scaled_result <- run_scale(x, method = "auto")
#' pca_result <- run_pca(scaled_result$data, metadata)
#' 
#' # PCA excluding QC samples
#' pca_result2 <- run_pca(scaled_result$data, metadata, 
#'                        group = "Group",
#'                        exclude = "QC")
#' 
#' # PCA with automatic scaling, excluding multiple groups
#' pca_result3 <- run_pca(x, metadata, 
#'                        scale_method = "auto",
#'                        group = "Group",
#'                        exclude = c("QC", "Blank"))
#' 
#' # PCA with transformation and scaling
#' pca_result4 <- run_pca(x, metadata, 
#'                        transform_method = "log2",
#'                        scale_method = "pareto",
#'                        group = "Group",
#'                        exclude = "QC")
#' 
#' # View results
#' print(pca_result)
#' summary(pca_result)
#' }
run_pca <- function(
    x,
    metadata,
    scale_method = NULL,
    transform_method = NULL,
    group = "Group",
    exclude = NULL,
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
  if (!is.character(group) || length(group) != 1) {
    stop("'group' must be a single character string specifying a column name in metadata")
  }
  if (!group %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", group))
  }
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      stop("'exclude' must be NULL or a character vector of group labels")
    }
    # Check if excluded groups exist in the data
    available_groups <- unique(metadata[[group]])
    invalid_groups <- exclude[!exclude %in% available_groups]
    if (length(invalid_groups) > 0) {
      warning("The following groups in 'exclude' were not found in metadata: ",
              paste(invalid_groups, collapse = ", "))
    }
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  
  # Handle sample exclusion
  n_samples_original <- nrow(x_matrix)
  excluded_groups <- character(0)
  
  if (!is.null(exclude)) {
    exclude_indices <- metadata[[group]] %in% exclude
    n_excluded <- sum(exclude_indices)
    
    if (n_excluded > 0) {
      msg(sprintf("Excluding %d samples from groups: %s", 
                  n_excluded, paste(exclude, collapse = ", ")))
      x_matrix <- x_matrix[!exclude_indices, , drop = FALSE]
      metadata <- metadata[!exclude_indices, , drop = FALSE]
      excluded_groups <- exclude[exclude %in% unique(metadata[[group]][exclude_indices])]
    } else {
      msg("No samples found matching exclusion criteria")
    }
  } else {
    n_excluded <- 0
    msg("No samples excluded")
  }
  
  if (nrow(x_matrix) == 0) {
    stop("All samples were excluded. No data remaining for PCA.")
  }
  
  # Check if data is from run_scale
  came_from_run_scale <- inherits(x, "run_scale") || 
    (!is.null(attr(x, "scaled:center")) && !is.null(attr(x, "scaled:scale")))
  
  preprocessing_log <- list(
    came_from_run_scale = came_from_run_scale,
    transformation_applied = FALSE,
    scaling_applied = FALSE,
    transform_method = NA,
    scale_method_used = NA
  )
  
  # Handle preprocessing
  if (!came_from_run_scale) {
    if (is.null(scale_method) && is.null(transform_method)) {
      warning("Data does not appear to be from run_scale() and no scale_method specified. ",
              "PCA will proceed on raw data. Consider using scale_method='auto' for proper scaling.")
    }
    
    # Apply transformation if requested
    if (!is.null(transform_method)) {
      msg(sprintf("Applying '%s' transformation before PCA...", transform_method))
      
      transform_result <- run_transform(
        x = x_matrix,
        method = transform_method,
        metadata = metadata,
        verbose = verbose
      )
      
      x_matrix <- as.matrix(transform_result$data)
      preprocessing_log$transformation_applied <- TRUE
      preprocessing_log$transform_method <- transform_method
    }
    
    # Apply scaling if requested
    if (!is.null(scale_method)) {
      msg(sprintf("Applying '%s' scaling before PCA...", scale_method))
      
      scale_result <- run_scale(
        x = x_matrix,
        method = scale_method,
        verbose = verbose
      )
      
      x_matrix <- as.matrix(scale_result$data)
      preprocessing_log$scaling_applied <- TRUE
      preprocessing_log$scale_method_used <- scale_method
    }
  } else {
    msg("Data appears to be pre-scaled. Proceeding with PCA...")
    if (!is.null(scale_method) || !is.null(transform_method)) {
      warning("Data is from run_scale() but scale_method/transform_method specified. ",
              "Additional preprocessing will be skipped to avoid double-scaling.")
    }
  }
  
  # Check for missing values
  if (any(is.na(x_matrix))) {
    stop("Data contains missing values (NA). PCA cannot proceed. ",
         "Please impute missing values using run_mvimpute() first.")
  }
  
  # Check for zero variance features
  feature_vars <- apply(x_matrix, 2, var, na.rm = TRUE)
  zero_var_features <- sum(feature_vars == 0 | is.na(feature_vars))
  if (zero_var_features > 0) {
    warning(sprintf("Found %d features with zero or NA variance. These may cause issues in PCA.",
                    zero_var_features))
  }
  
  msg("Performing PCA...")
  
  # Perform PCA (data should already be centered and scaled appropriately)
  pca_result <- tryCatch({
    stats::prcomp(x_matrix, center = FALSE, scale. = FALSE)
  }, error = function(e) {
    stop("PCA failed: ", e$message, 
         "\nThis may be due to:\n",
         "  - Insufficient samples relative to features\n",
         "  - Zero variance features\n",
         "  - Numerical instabilities\n",
         "Consider feature filtering or checking data quality.")
  })
  
  # Extract results
  scores <- pca_result$x
  loadings <- pca_result$rotation
  
  # Calculate variance metrics
  eigenvalues <- pca_result$sdev^2
  total_variance <- sum(eigenvalues)
  variance_explained <- (eigenvalues / total_variance) * 100
  cumulative_variance <- cumsum(variance_explained)
  
  n_pcs <- ncol(scores)
  n_samples <- nrow(scores)
  n_features <- ncol(x_matrix)
  
  msg(sprintf("PCA complete: %d PCs computed from %d samples and %d features",
              n_pcs, n_samples, n_features))
  msg(sprintf("First 5 PCs explain %.1f%% of total variance",
              sum(variance_explained[1:min(5, n_pcs)])))
  
  # Align metadata with scores
  if (!is.null(rownames(x_matrix))) {
    rownames(scores) <- rownames(x_matrix)
  }
  
  # Return results
  result <- list(
    scores = scores,
    loadings = loadings,
    variance_explained = variance_explained,
    eigenvalues = eigenvalues,
    cumulative_variance = cumulative_variance,
    center_values = NULL,
    scale_values = NULL,
    pca_object = pca_result,
    data_used = x_matrix,
    metadata = metadata,
    n_samples = n_samples,
    n_samples_excluded = n_excluded,
    n_features = n_features,
    n_pcs = n_pcs,
    excluded_groups = excluded_groups,
    preprocessing = preprocessing_log,
    parameters = list(
      scale_method = scale_method,
      transform_method = transform_method,
      group = group,
      exclude = exclude
    )
  )
  
  class(result) <- c("run_pca", "list")
  return(result)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' @export
print.run_pca <- function(x, ...) {
  cat("=== Principal Component Analysis ===\n")
  cat("Samples:     ", x$n_samples, "\n")
  if (x$n_samples_excluded > 0) {
    cat("Excluded:    ", x$n_samples_excluded, " (", 
        paste(x$excluded_groups, collapse = ", "), ")\n", sep = "")
  }
  cat("Features:    ", x$n_features, "\n")
  cat("PCs computed:", x$n_pcs, "\n")
  cat("\nVariance explained by first 5 PCs:\n")
  for (i in 1:min(5, x$n_pcs)) {
    cat(sprintf("  PC%d: %.1f%% (cumulative: %.1f%%)\n",
                i, x$variance_explained[i], x$cumulative_variance[i]))
  }
  
  if (x$preprocessing$transformation_applied) {
    cat(sprintf("\nTransformation: %s\n", x$preprocessing$transform_method))
  }
  if (x$preprocessing$scaling_applied) {
    cat(sprintf("Scaling: %s\n", x$preprocessing$scale_method_used))
  }
  
  cat("\nUse plot_scree() to visualize variance explained\n")
  cat("Use plot_score() to visualize sample positions\n")
  invisible(x)
}

#' @export
summary.run_pca <- function(object, ...) {
  var_table <- data.frame(
    PC = paste0("PC", 1:min(10, object$n_pcs)),
    Variance_Explained = object$variance_explained[1:min(10, object$n_pcs)],
    Cumulative = object$cumulative_variance[1:min(10, object$n_pcs)],
    Eigenvalue = object$eigenvalues[1:min(10, object$n_pcs)]
  )
  
  ans <- list(
    n_samples = object$n_samples,
    n_samples_excluded = object$n_samples_excluded,
    excluded_groups = object$excluded_groups,
    n_features = object$n_features,
    n_pcs = object$n_pcs,
    variance_table = var_table,
    preprocessing = object$preprocessing
  )
  
  class(ans) <- "summary.run_pca"
  return(ans)
}

#' @export
print.summary.run_pca <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("PCA Summary\n")
  cat("---------------------------------------\n")
  cat(sprintf("Samples: %d | Features: %d | PCs: %d\n",
              x$n_samples, x$n_features, x$n_pcs))
  if (x$n_samples_excluded > 0) {
    cat(sprintf("Excluded: %d samples (%s)\n", 
                x$n_samples_excluded, 
                paste(x$excluded_groups, collapse = ", ")))
  }
  cat("\nVariance Explained (Top 10 PCs):\n")
  print(x$variance_table, row.names = FALSE)
  
  if (x$preprocessing$transformation_applied || x$preprocessing$scaling_applied) {
    cat("\n-- Preprocessing Applied --\n")
    if (x$preprocessing$transformation_applied) {
      cat("Transformation:", x$preprocessing$transform_method, "\n")
    }
    if (x$preprocessing$scaling_applied) {
      cat("Scaling:", x$preprocessing$scale_method_used, "\n")
    }
  }
  invisible(x)
}