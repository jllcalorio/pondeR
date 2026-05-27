#' Perform LEfSe (Linear Discriminant Analysis Effect Size) analysis
#'
#' This function performs LEfSe analysis to identify differentially abundant features
#' between classes (e.g., biomarkers) with effect size estimation.
#'
#' @param x A data frame or matrix of features (samples x features). Should be
#'   relative abundance data (proportions summing to 1 per sample).
#' @param metadata A data frame containing sample metadata. Must include the
#'   column specified by \code{class} and optionally \code{subclass}.
#' @param class A single character string specifying the column name in \code{metadata}
#'   that defines the main classes for comparison (e.g., disease status).
#' @param subclass Optional. A single character string specifying the column name
#'   in \code{metadata} that defines subclasses (e.g., age groups, gender) for
#'   stratified analysis. If provided, the analysis is performed within each subclass.
#' @param alpha Numeric value between 0 and 1 specifying the significance level
#'   for the Kruskal-Wallis test (default: 0.05).
#' @param wilcox_alpha Numeric value between 0 and 1 specifying the significance
#'   level for the pairwise Wilcoxon test (default: 0.05).
#' @param lda_threshold Numeric value specifying the threshold on the absolute
#'   logarithmic LDA score for discriminative features (default: 2.0).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#' @return A list of class \code{run_lefse} containing:
#'   \item{lda}{A data frame with LDA results (features, means per class, LDA scores).}
#'   \item{summary_table}{A data frame summarizing significant features with
#'     effect sizes, p-values, and classifications.}
#'   \item{parameters}{A list of the parameters used in the analysis.}
#'   \item{call}{The matched call.}
#' @details
#' The LEfSe algorithm consists of three steps:
#' \enumerate{
#'   \item Kruskal-Wallis test to detect significantly abundant features between classes.
#'   \item Pairwise Wilcoxon test to check consistency of biological trends
#'     across subclasses within each class.
#'   \item Linear Discriminant Analysis to estimate the effect size of each
#'     significantly abundant feature.
#' }
#' The input data should ideally be relative abundances (e.g., from
#' \code{run_relabund} or similar). The function assumes that higher values
#' indicate higher abundance.
#' @examples
#' # Load example data (if available) or use your own
#' # data(metagenomicexample, package = "pondeR")
#' \dontrun{
#'   # Assuming `otu_table` is a feature table and `meta` is metadata
#'   lefse_res <- run_lefse(
#'     x = otu_table,
#'     metadata = meta,
#'     class = "Condition",
#'     subclass = "AgeGroup",
#'     alpha = 0.05,
#'     wilcox_alpha = 0.05,
#'     lda_threshold = 2.0
#'   )
#' }
#' @seealso \code{\link{plot_lefse}} for visualization (if implemented),
#'   \code{\link{run_assoc}} for association testing.
#' @encoding UTF-8
#' @export
run_lefse <- function(x,
                      metadata,
                      class = "Group",
                      subclass = NULL,
                      alpha = 0.05,
                      wilcox_alpha = 0.05,
                      lda_threshold = 2.0,
                      verbose = TRUE) {
  
  # Define internal message function
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("'x' must be a data frame or matrix", call. = FALSE)
  }
  
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("'x' is a matrix but has no column names", call. = FALSE)
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  
  if (nrow(x) == 0) {
    stop("'x' has zero rows", call. = FALSE)
  }
  
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame", call. = FALSE)
  }
  
  if (nrow(x) != nrow(metadata)) {
    stop("Number of rows in 'x' must equal number of rows in 'metadata'", call. = FALSE)
  }
  
  if (!is.character(class) || length(class) != 1) {
    stop("'class' must be a single character string specifying a column name in metadata", call. = FALSE)
  }
  
  if (!class %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", class), call. = FALSE)
  }
  
  if (!is.null(subclass)) {
    if (!is.character(subclass) || length(subclass) != 1) {
      stop("'subclass' must be NULL or a single character string specifying a column name in metadata", call. = FALSE)
    }
    if (!subclass %in% colnames(metadata)) {
      stop(sprintf("Column '%s' not found in metadata", subclass), call. = FALSE)
    }
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value in (0, 1]", call. = FALSE)
  }
  
  if (!is.numeric(wilcox_alpha) || length(wilcox_alpha) != 1 || is.na(wilcox_alpha) || wilcox_alpha <= 0 || wilcox_alpha > 1) {
    stop("'wilcox_alpha' must be a single numeric value in (0, 1]", call. = FALSE)
  }
  
  if (!is.numeric(lda_threshold) || length(lda_threshold) != 1 || is.na(lda_threshold) || lda_threshold <= 0) {
    stop("'lda_threshold' must be a single positive numeric value", call. = FALSE)
  }
  
  # Check for missing values
  if (any(is.na(x))) {
    warning("Input data contains missing values. Consider imputation before running LEfSe.", call. = FALSE)
  }
  
  # Check for negative values (relative abundances should be non-negative)
  if (any(x < 0, na.rm = TRUE)) {
    warning("Input data contains negative values. LEfSe assumes non-negative abundances.", call. = FALSE)
  }
  
  msg("Starting LEfSe analysis...")
  msg("Features: ", ncol(x), ", Samples: ", nrow(x))
  
  # Prepare data
  class_vec <- metadata[[class]]
  if (!is.null(subclass)) {
    subclass_vec <- metadata[[subclass]]
  } else {
    subclass_vec <- rep("NA", nrow(metadata))
  }
  
  # Step 1: Kruskal-Wallis test for main class effect
  msg("Step 1: Kruskal-Wallis test for class effects...")
  kw_pvals <- apply(x, 2, function(feature) {
    # Remove NA values for the test
    complete_cases <- !is.na(feature)
    if (sum(complete_cases) < 3) return(NA)  # Need at least 3 observations
    
    # Group by class
    groups <- split(feature[complete_cases], class_vec[complete_cases])
    # Remove empty groups
    groups <- groups[sapply(groups, length) > 0]
    
    if (length(groups) < 2) return(NA)  # Need at least 2 groups
    
    # Perform Kruskal-Wallis test
    kruskal.test(feature ~ class_vec[complete_cases])$p.value
  })
  
  # Features passing KW test
  kw_sig <- !is.na(kw_pvals) & (kw_pvals < alpha)
  msg("  Features passing Kruskal-Wallis test (p < ", alpha, "): ", sum(kw_sig, na.rm = TRUE))
  
  if (sum(kw_sig, na.rm = TRUE) == 0) {
    warning("No features passed the Kruskal-Wallis test. Consider adjusting alpha.", call = FALSE)
    return(structure(list(
      lda = data.frame(),
      summary_table = data.frame(),
      parameters = list(class = class, subclass = subclass, alpha = alpha,
                        wilcox_alpha = wilcox_alpha, lda_threshold = lda_threshold),
      call = match.call()
    ), class = "run_lefse"))
  }
  
  # Step 2: Pairwise Wilcoxon test for subclass consistency (if subclass provided)
  msg("Step 2: Pairwise Wilcoxon test for subclass consistency...")
  wilcox_pass <- rep(TRUE, ncol(x))  # Initialize all as TRUE
  
  if (!is.null(subclass) && any(!is.na(subclass_vec))) {
    # Get unique subclasses
    unique_subclasses <- unique(subclass_vec[!is.na(subclass_vec)])
    msg("  Checking consistency across subclasses: ", paste(unique_subclasses, collapse = ", "))
    
    # For each feature passing KW, check pairwise Wilcoxon within each main class
    for (feat_idx in which(kw_sig)) {
      feature_vals <- x[, feat_idx]
      passed <- TRUE
      
      # For each main class
      unique_classes <- unique(class_vec[!is.na(class_vec)])
      for (cl in unique_classes) {
        cl_mask <- (class_vec == cl) & !is.na(class_vec)
        if (sum(cl_mask) < 2) next  # Need at least 2 samples
        
        # Get subclass values for this class
        cl_feature <- feature_vals[cl_mask]
        cl_subclass <- subclass_vec[cl_mask]
        
        # Check each pair of subclasses
        sub_classes <- unique(cl_subclass[!is.na(cl_subclass)])
        if (length(sub_classes) >= 2) {
          # Perform pairwise Wilcoxon tests
          for (i in 1:(length(sub_classes)-1)) {
            for (j in (i+1):length(sub_classes)) {
              sub1_mask <- (cl_subclass == sub_classes[i]) & !is.na(cl_subclass)
              sub2_mask <- (cl_subclass == sub_classes[j]) & !is.na(cl_subclass)
              
              if (sum(sub1_mask) >= 2 && sum(sub2_mask) >= 2) {
                # Wilcoxon rank sum test
                wilcox_p <- tryCatch({
                  wilcox.test(cl_feature[sub1_mask], cl_feature[sub2_mask])$p.value
                }, error = function(e) NA)
                
                if (!is.na(wilcox_p) && (wilcox_p >= wilcox_alpha)) {
                  passed <- FALSE
                  break
                }
              }
            }
            if (!passed) break
          }
        }
        if (!passed) break
      }
      wilcox_pass[feat_idx] <- passed
    }
    
    msg("  Features passing Wilcoxon test: ", sum(wilcox_pass & kw_sig, na.rm = TRUE))
  }
  
  # Features passing both tests
  final_sig <- kw_sig & wilcox_pass
  msg("  Final significant features: ", sum(final_sig, na.rm = TRUE))
  
  if (sum(final_sig, na.rm = TRUE) == 0) {
    warning("No features passed both statistical tests. Consider adjusting parameters.", call = FALSE)
    return(structure(list(
      lda = data.frame(),
      summary_table = data.frame(),
      parameters = list(class = class, subclass = subclass, alpha = alpha,
                        wilcox_alpha = wilcox_alpha, lda_threshold = lda_threshold),
      call = match.call()
    ), class = "run_lefse"))
  }
  
  # Step 3: Linear Discriminant Analysis
  msg("Step 3: Linear Discriminant Analysis for effect size...")
  
  # Prepare data for LDA (use features passing tests)
  lda_features <- x[, final_sig, drop = FALSE]
  lda_class <- class_vec
  
  # Remove samples with missing class
  complete_mask <- !is.na(lda_class)
  if (sum(complete_mask) < 2) {
    warning("Insufficient samples with class information for LDA.", call = FALSE)
    return(structure(list(
      lda = data.frame(),
      summary_table = data.frame(),
      parameters = list(class = class, subclass = subclass, alpha = alpha,
                        wilcox_alpha = wilcox_alpha, lda_threshold = lda_threshold),
      call = match.call()
    ), class = "run_lefse"))
  }
  
  lda_features <- lda_features[complete_mask, , drop = FALSE]
  lda_class <- lda_class[complete_mask]
  
  # Check if we have enough classes
  unique_classes <- unique(lda_class)
  if (length(unique_classes) < 2) {
    warning("Need at least 2 classes for LDA.", call = FALSE)
    return(structure(list(
      lda = data.frame(),
      summary_table = data.frame(),
      parameters = list(class = class, subclass = subclass, alpha = alpha,
                        wilcox_alpha = wilcox_alpha, lda_threshold = lda_threshold),
      call = match.call()
    ), class = "run_lefse"))
  }
  
  # Simple LDA implementation using MASS::lda if available, otherwise manual
  lda_results <- tryCatch({
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for LDA. Install with: install.packages('MASS')", call = FALSE)
    }
    
    # MASS::lda expects a formula or matrix and grouping
    lda_obj <- MASS::lda(lda_class ~ ., data = as.data.frame(lda_features))
    
    # Extract results
    lda_df <- data.frame(
      feature = colnames(lda_features),
      stringsAsFactors = FALSE
    )
    
    # Add class means
    for (cl in unique_classes) {
      cl_mean <- sapply(lda_features[lda_class == cl, , drop = FALSE], mean, na.rm = TRUE)
      lda_df[[paste0("mean_", cl)]] <- cl_mean
    }
    
    # Add LDA scores (scaling)
    lda_scores <- predict(lda_obj)$x  # This gives the discriminant scores
    # For multi-class, we take the absolute value of the first discriminant
    # or compute something similar to LEfSe's approach
    if (ncol(lda_scores) >= 1) {
      lda_df$lda_score <- abs(lda_scores[, 1])  # First discriminant component
    } else {
      lda_df$lda_score <- 0
    }
    
    lda_df
  }, error = function(e) {
    warning("LDA failed: ", e$message, ". Using effect size based on fold changes.", call = FALSE)
    # Fallback: compute effect size as log2(fold change) between group means
    lda_df <- data.frame(
      feature = colnames(lda_features),
      stringsAsFactors = FALSE
    )
    
    # Calculate mean for each class
    class_means <- sapply(unique_classes, function(cl) {
      colMeans(lda_features[lda_class == cl, , drop = FALSE], na.rm = TRUE)
    })
    colnames(class_means) <- unique_classes
    
    # Add to dataframe
    for (cl in unique_classes) {
      lda_df[[paste0("mean_", cl)]] <- class_means[, cl]
    }
    
    # Compute effect size as max(abs(log2(fold change))) between any two classes
    effect_sizes <- apply(lda_features, 2, function(feature_vals) {
      # Get means for each class
      means <- tapply(feature_vals, lda_class, mean, na.rm = TRUE)
      # Compute all pairwise log2 fold changes
      fc_vals <- c()
      if (length(means) >= 2) {
        for (i in 1:(length(means)-1)) {
          for (j in (i+1):length(means)) {
            if (means[i] > 0 && means[j] > 0) {
              fc_vals <- c(fc_vals, abs(log2(means[i] / means[j])))
            }
          }
        }
      }
      if (length(fc_vals) > 0) {
        max(fc_vals, na.rm = TRUE)
      } else {
        0
      }
    })
    
    lda_df$lda_score <- effect_sizes
    lda_df
  })
  
  # Filter by LDA threshold
  lda_significant <- lda_results$lda_score >= lda_threshold
  msg("  Features with LDA score >= ", lda_threshold, ": ", sum(lda_significant))
  
  # Prepare summary table
  summary_table <- lda_results[lda_significant, , drop = FALSE]
  
  # Add p-value from KW test
  summary_table$pvalue <- kw_pvals[final_sig][lda_significant]
  
# Order results by LDA score (largest first)
  ord <- order(-summary_table$lda_score)
  summary_table <- summary_table[ord, ]
  
  # Create final structure
  result <- list(
    lda = lda_results,
    summary_table = summary_table,
    parameters = list(
      class = class,
      subclass = subclass,
      alpha = alpha,
      wilcox_alpha = wilcox_alpha,
      lda_threshold = lda_threshold
    ),
    call = match.call()
  )
  class(result) <- "run_lefse"
  return(result)
}