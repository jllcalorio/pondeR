#' Perform Partial Least Squares (PLS) Analysis
#'
#' @description
#' Performs Partial Least Squares Discriminant Analysis (PLS-DA), Orthogonal 
#' PLS-DA (OPLS-DA), or sparse PLS-DA (sPLS-DA) on preprocessed metabolomics data.
#' Optionally applies scaling and transformation before analysis if the input data
#' has not been pre-scaled. Designed to match the structure of `run_pca()`.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   Can be output from `run_scale()` or raw data (which will be scaled if `scale_method`
#'   is specified).
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain columns for grouping and sample identifiers.
#' @param method Character. The PLS method to use. Options are `"oplsda"` (default),
#'   `"plsda"`, or `"splsda"`.
#' @param group Character. Name of column in `metadata` containing group labels for
#'   classification and sample exclusion. Required parameter. Default: "Group".
#' @param exclude Character vector or NULL. Group labels (from the `group` column) to exclude
#'   from the analysis. For example, c("QC", "Blank") excludes all samples with these group
#'   labels. Default: NULL (include all samples).
#' @param scale_method Character or NULL. Scaling method to apply before PLS if `x` is not
#'   already scaled. Options: "mean", "auto", "pareto". If NULL (default), the function assumes 
#'   the data is appropriately scaled. Pareto scaling is generally recommended for PLS methods.
#' @param transform_method Character or NULL. Transformation method to apply before scaling
#'   and PLS. Options: "log2", "log10", "sqrt", "cbrt", "vsn", "glog". Default: NULL (no transformation).
#' @param ncomp Integer. Number of predictive components to compute for PLS-DA and sPLS-DA. 
#'   Default: 2.
#' @param orthoI Integer or NA. Number of orthogonal components for OPLS-DA. 
#'   If NA (default), the optimal number is determined automatically by cross-validation.
#' @param crossvalI Integer. Number of cross-validation folds. Default: 10.
#' @param permI Integer. Number of permutations for model validation. Default: 20.
#' @param keepX Integer vector. For sPLS-DA only: specifies the number of variables to keep
#'   on each component. Default: NULL (automatic selection or keeps all).
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#'   \item{scores}{Matrix of scores (samples × components), aliased as PC1, PC2... for `plot_score` compatibility}
#'   \item{loadings}{Matrix of feature loadings (features × components)}
#'   \item{variance_explained}{Numeric vector of variance explained (%) per component}
#'   \item{vip_scores}{Numeric vector of Variable Importance in Projection (VIP) scores (PLS-DA/OPLS-DA only)}
#'   \item{pls_object}{Original model object from `ropls` or `mixOmics`}
#'   \item{data_used}{Matrix of data used for PLS (after any transformation/scaling)}
#'   \item{metadata}{Data frame of metadata aligned with scores}
#'   \item{n_samples}{Integer. Number of samples used}
#'   \item{n_samples_excluded}{Integer. Number of samples excluded}
#'   \item{n_features}{Integer. Number of features}
#'   \item{n_pcs}{Integer. Number of components computed}
#'   \item{excluded_groups}{Character vector of excluded group labels}
#'   \item{preprocessing}{List documenting preprocessing steps applied}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Compatibility:**
#' 
#' The output of this function is specifically designed to be compatible with `plot_score()`. 
#' Component names in the `scores` matrix are aliased as `PC1`, `PC2`, etc., to allow seamless 
#' plotting. For OPLS-DA, "PC1" represents the Predictive Component, and "PC2" represents the 
#' first Orthogonal Component.
#' 
#' **Methods:**
#' 
#' - **oplsda**: Orthogonal PLS-DA separates predictive variation (correlated to class) from 
#'   orthogonal variation (uncorrelated). Excellent for biomarker discovery. Utilizes `ropls::opls()`.
#' - **plsda**: Standard PLS-DA. Maximizes covariance between features and class labels.
#'   Utilizes `ropls::opls()`.
#' - **splsda**: Sparse PLS-DA incorporates L1 penalization to perform feature selection 
#'   during model building. Utilizes `mixOmics::splsda()`.
#'
#' **Sample Exclusion:**
#' 
#' QC samples should typically be excluded from discriminant analysis model building. 
#' Use the `exclude` parameter (e.g., `exclude = "QC"`) to drop them prior to scaling and fitting.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Thevenot, E. A., Roux, A., Xu, Y., Ezan, E., & Junot, C. (2015). Analysis of the Human Adult Urinary 
#' Metabolome Variations with Age, Body Mass Index, and Gender by Implementing a Comprehensive Workflow 
#' for Univariate and OPLS Statistical Analyses. Journal of Proteome Research, 14(8), 3322-3335.
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. (2017) mixOmics: An R package for 'omics feature selection 
#' and multiple data integration. PLoS Comput Biol 13(11): e1005752.
#'
#' @importFrom stats var
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 50, mean = 100, sd = 20), nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' rownames(x) <- paste0("Sample", 1:100)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("Sample", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
#' )
#' 
#' # Run PLS-DA excluding QC samples, with automatic pareto scaling
#' pls_result <- run_pls(x, metadata, 
#'                       method = "plsda",
#'                       scale_method = "pareto",
#'                       group = "Group",
#'                       exclude = "QC")
#' 
#' # Visualize with plot_score
#' plot_score(pls_result, color_by = "Group", points_from = "Sample")
#' 
#' # Run OPLS-DA
#' opls_result <- run_pls(x, metadata, 
#'                        method = "oplsda",
#'                        scale_method = "pareto",
#'                        group = "Group",
#'                        exclude = "QC")
#' }
run_pls <- function(
    x,
    metadata,
    method = "oplsda",
    group = "Group",
    exclude = NULL,
    scale_method = NULL,
    transform_method = NULL,
    ncomp = 2,
    orthoI = NA,
    crossvalI = 10,
    permI = 20,
    keepX = NULL,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  if (nrow(x) != nrow(metadata)) {
    stop("Number of rows in 'x' must equal number of rows in 'metadata'")
  }
  
  method <- tolower(method)
  valid_methods <- c("oplsda", "plsda", "splsda")
  if (!method %in% valid_methods) {
    stop(sprintf("'method' must be one of: %s", paste(valid_methods, collapse = ", ")))
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
    available_groups <- unique(metadata[[group]])
    invalid_groups <- exclude[!exclude %in% available_groups]
    if (length(invalid_groups) > 0) {
      warning("The following groups in 'exclude' were not found in metadata: ",
              paste(invalid_groups, collapse = ", "))
    }
  }
  
  if (method %in% c("oplsda", "plsda") && !requireNamespace("ropls", quietly = TRUE)) {
    stop("Package 'ropls' is required for PLS-DA and OPLS-DA. Install it from Bioconductor.")
  }
  if (method == "splsda" && !requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required for sPLS-DA. Install it from Bioconductor.")
  }
  
  was_matrix <- is.matrix(x)
  x_matrix <- as.matrix(x)
  
  # ============================================================================
  # SAMPLE EXCLUSION & PREPROCESSING
  # ============================================================================
  
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
    stop("All samples were excluded. No data remaining for analysis.")
  }
  
  # Check if data is from run_scale (which leaves scale attributes)
  came_from_run_scale <- inherits(x, "run_scale") || 
    (!is.null(attr(x, "scaled:center")) && !is.null(attr(x, "scaled:scale")))
  
  preprocessing_log <- list(
    came_from_run_scale = came_from_run_scale,
    transformation_applied = FALSE,
    scaling_applied = FALSE,
    transform_method = NA,
    scale_method_used = NA
  )
  
  if (!came_from_run_scale) {
    if (is.null(scale_method) && is.null(transform_method)) {
      warning("Data does not appear to be from run_scale() and no scale_method specified. ",
              "PLS will proceed on raw data. Consider using scale_method='pareto'.")
    }
    
    # Transformation
    if (!is.null(transform_method)) {
      msg(sprintf("Applying '%s' transformation...", transform_method))
      if (!exists("run_transform", mode = "function")) {
         stop("Helper function 'run_transform' not found. Please load it into your environment.")
      }
      transform_result <- run_transform(x = x_matrix, method = transform_method, 
                                        metadata = metadata, verbose = verbose)
      x_matrix <- as.matrix(transform_result$data)
      preprocessing_log$transformation_applied <- TRUE
      preprocessing_log$transform_method <- transform_method
    }
    
    # Scaling
    if (!is.null(scale_method)) {
      msg(sprintf("Applying '%s' scaling...", scale_method))
      if (!exists("run_scale", mode = "function")) {
        stop("Helper function 'run_scale' not found. Please load it into your environment.")
      }
      scale_result <- run_scale(x = x_matrix, method = scale_method, verbose = verbose)
      x_matrix <- as.matrix(scale_result$data)
      preprocessing_log$scaling_applied <- TRUE
      preprocessing_log$scale_method_used <- scale_method
    }
  } else {
    msg("Data appears to be pre-scaled. Proceeding with PLS...")
    if (!is.null(scale_method) || !is.null(transform_method)) {
      warning("Data is pre-scaled but scale_method/transform_method specified. ",
              "Additional preprocessing will be skipped to avoid double-scaling.")
    }
  }
  
  if (any(is.na(x_matrix))) stop("Data contains missing values (NA). Impute values first.")
  
  feature_vars <- apply(x_matrix, 2, stats::var, na.rm = TRUE)
  if (sum(feature_vars == 0 | is.na(feature_vars)) > 0) {
    warning("Found features with zero variance. These may cause numerical instabilities.")
  }
  
  # ============================================================================
  # MODEL FITTING
  # ============================================================================
  
  Y <- as.factor(metadata[[group]])
  
  if (nlevels(Y) < 2) {
    stop("At least 2 unique groups are required for Discriminant Analysis.")
  }
  
  msg(sprintf("Performing %s...", toupper(method)))
  
  pls_model <- NULL
  scores <- NULL
  loadings <- NULL
  variance_explained <- NULL
  vip_scores <- NULL
  
  if (method == "oplsda") {
    
    if (nlevels(Y) > 2) {
      warning("OPLS-DA typically requires exactly 2 classes. 'ropls' may fail or fall back to PLS-DA.")
    }
    
    pls_model <- tryCatch({
      ropls::opls(x_matrix, Y, predI = 1, orthoI = orthoI, crossvalI = crossvalI, 
                  permI = permI, scaleC = "none", 
                  #fig.mac = FALSE, 
                  info.txtC = "none")
    }, error = function(e) stop("ropls::opls failed: ", e$message))
    
    # Combine Predictive and Orthogonal scores
    if (is.null(pls_model@orthoScoreMN)) {
      warning("No orthogonal components found. The model is functionally equivalent to PLS-DA.")
      scores <- pls_model@scoreMN
      loadings <- pls_model@loadingMN
      variance_explained <- pls_model@modelDF$R2X * 100
    } else {
      scores <- cbind(pls_model@scoreMN, pls_model@orthoScoreMN)
      loadings <- cbind(pls_model@loadingMN, pls_model@orthoLoadingMN)
      
      # Extract R2X (Predictive then Orthogonal)
      r2x_pred <- pls_model@modelDF["p1", "R2X"] * 100
      r2x_ortho <- pls_model@modelDF[grep("^o", rownames(pls_model@modelDF)), "R2X"] * 100
      variance_explained <- c(r2x_pred, r2x_ortho)
    }
    vip_scores <- ropls::getVipVn(pls_model)
    
  } else if (method == "plsda") {
    
    pls_model <- tryCatch({
      ropls::opls(x_matrix, Y, predI = ncomp, orthoI = 0, crossvalI = crossvalI, 
                  permI = permI, scaleC = "none", 
                  #fig.mac = FALSE, 
                  info.txtC = "none")
    }, error = function(e) stop("ropls::opls failed: ", e$message))
    
    scores <- pls_model@scoreMN
    loadings <- pls_model@loadingMN
    variance_explained <- pls_model@modelDF$R2X * 100
    vip_scores <- ropls::getVipVn(pls_model)
    
  } else if (method == "splsda") {
    
    pls_model <- tryCatch({
      mixOmics::splsda(X = x_matrix, Y = Y, ncomp = ncomp, keepX = keepX)
    }, error = function(e) stop("mixOmics::splsda failed: ", e$message))
    
    scores <- pls_model$variates$X
    loadings <- pls_model$loadings$X
    variance_explained <- pls_model$prop_expl_var$X * 100
    vip_scores <- NULL # sPLS-DA uses loadings for variable selection instead of VIP
  }
  
  # ============================================================================
  # FORMAT OUTPUT FOR plot_score COMPATIBILITY
  # ============================================================================
  
  # n_pcs <- ncol(scores)
  # n_samples <- nrow(scores)
  # n_features <- ncol(x_matrix)
  
  # # Crucial Step: plot_score expects columns addressable via 'paste0("PC", pc)' 
  # colnames(scores) <- paste0("PC", seq_len(n_pcs))
  
  # msg(sprintf("%s complete: %d components computed from %d samples and %d features",
  #             toupper(method), n_pcs, n_samples, n_features))

  # Ensure scores is a matrix/df before checking ncol
  if (is.null(scores)) {
    stop("The PLS model failed to produce scores. Check if groups are valid or if data has variance.")
  }
  
  # If scores is a vector (single component), force it to a matrix
  if (is.null(dim(scores))) {
    scores <- matrix(scores, ncol = 1)
  }

  n_pcs <- ncol(scores)
  
  # Safety: Only rename if we actually have columns
  if (n_pcs > 0) {
    colnames(scores) <- paste0("PC", seq_len(n_pcs))
  } else {
    stop("No components were computed. The model may be invalid.")
  }
  
  n_samples <- nrow(scores)
  n_features <- ncol(x_matrix)
  
  msg(sprintf("%s complete: %d components computed from %d samples and %d features",
              toupper(method), n_pcs, n_samples, n_features))
  
  if (!is.null(rownames(x_matrix))) rownames(scores) <- rownames(x_matrix)
  
  result <- list(
    scores = as.data.frame(scores),
    loadings = loadings,
    variance_explained = unname(variance_explained),
    vip_scores = vip_scores,
    pls_object = pls_model,
    data_used = x_matrix,
    metadata = metadata,
    n_samples = n_samples,
    n_samples_excluded = n_excluded,
    n_features = n_features,
    n_pcs = n_pcs,
    excluded_groups = excluded_groups,
    preprocessing = preprocessing_log,
    method_used = method,
    parameters = list(
      scale_method = scale_method,
      transform_method = transform_method,
      group = group,
      exclude = exclude,
      ncomp = ncomp,
      orthoI = orthoI
    )
  )
  
  class(result) <- c("run_pls", "list")
  return(result)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' @export
print.run_pls <- function(x, ...) {
  cat(sprintf("=== Partial Least Squares Analysis (%s) ===\n", toupper(x$method_used)))
  cat("Samples:     ", x$n_samples, "\n")
  if (x$n_samples_excluded > 0) {
    cat("Excluded:    ", x$n_samples_excluded, " (", 
        paste(x$excluded_groups, collapse = ", "), ")\n", sep = "")
  }
  cat("Features:    ", x$n_features, "\n")
  cat("Components:  ", x$n_pcs, "\n")
  cat("\nVariance explained by components:\n")
  
  comp_names <- if (x$method_used == "oplsda") {
    c("Predictive", rep("Orthogonal", x$n_pcs - 1))
  } else {
    paste("Component", seq_len(x$n_pcs))
  }
  
  for (i in seq_len(min(5, x$n_pcs))) {
    cat(sprintf("  %s: %.1f%%\n", comp_names[i], x$variance_explained[i]))
  }
  
  if (x$preprocessing$transformation_applied) {
    cat(sprintf("\nTransformation: %s\n", x$preprocessing$transform_method))
  }
  if (x$preprocessing$scaling_applied) {
    cat(sprintf("Scaling: %s\n", x$preprocessing$scale_method_used))
  }
  
  cat("\nUse plot_score() to visualize sample positions.\n")
  invisible(x)
}

#' @export
summary.run_pls <- function(object, ...) {
  comp_names <- if (object$method_used == "oplsda") {
    c("Predictive", paste0("Orthogonal_", seq_len(object$n_pcs - 1)))
  } else {
    paste0("Comp_", seq_len(object$n_pcs))
  }
  
  var_table <- data.frame(
    Component = comp_names[1:min(10, object$n_pcs)],
    Variance_Explained = object$variance_explained[1:min(10, object$n_pcs)]
  )
  
  ans <- list(
    method = object$method_used,
    n_samples = object$n_samples,
    n_samples_excluded = object$n_samples_excluded,
    excluded_groups = object$excluded_groups,
    n_features = object$n_features,
    n_pcs = object$n_pcs,
    variance_table = var_table,
    preprocessing = object$preprocessing
  )
  
  class(ans) <- "summary.run_pls"
  return(ans)
}

#' @export
print.summary.run_pls <- function(x, ...) {
  cat("---------------------------------------\n")
  cat(sprintf("%s Summary\n", toupper(x$method)))
  cat("---------------------------------------\n")
  cat(sprintf("Samples: %d | Features: %d | Components: %d\n",
              x$n_samples, x$n_features, x$n_pcs))
  if (x$n_samples_excluded > 0) {
    cat(sprintf("Excluded: %d samples (%s)\n", 
                x$n_samples_excluded, 
                paste(x$excluded_groups, collapse = ", ")))
  }
  cat("\nVariance Explained:\n")
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