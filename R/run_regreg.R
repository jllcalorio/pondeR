#' Perform Regularized Regression Analysis
#'
#' @title Regularized Regression with Cross-Validation (Ridge, Elastic Net, LASSO)
#'
#' @description
#' Fits regularized regression models — Ridge (alpha = 0), LASSO (alpha = 1), or
#' Elastic Net (0 < alpha < 1) — with k-fold cross-validation via
#' \code{\link[glmnet]{cv.glmnet}}. Supports both categorical (binary and multinomial
#' classification) and continuous (Gaussian regression) dependent variables.
#' When multiple alpha values are supplied, all are evaluated on identical
#' cross-validation folds and the best-performing model is returned alongside
#' a comparison table.
#'
#' Unpenalized covariates (e.g., batch indicators, demographic confounders) can be
#' included through the \code{not_penalized} argument; their coefficients are
#' estimated freely while all other predictors remain subject to regularization.
#'
#' @details
#' \strong{Data assembly:}
#' The feature matrix is taken from \code{x}. If \code{not_penalized} column names
#' belonging to \code{metadata} are provided, those columns are appended to the
#' feature matrix. A \code{penalty.factor} vector of length \code{ncol(x)} is set to
#' \code{1} for all penalized predictors and \code{0} for unpenalized ones, and is
#' passed to \code{glmnet::cv.glmnet}.
#'
#' \strong{Cross-validation consistency:}
#' Fold IDs are generated once (using \code{seed}) and reused across all tested alpha
#' values, ensuring a fair, apple-to-apple comparison of regularization paths.
#'
#' \strong{Multinomial coefficients:}
#' For multinomial models the raw coefficient matrix is normalized relative to
#' the reference level so that each comparison reflects log-odds vs. the reference.
#' Custom pairwise contrasts are not supported in this version; use the
#' \code{Coefficients} data frame and compute differences manually.
#'
#' \strong{Why no p-values?}
#' P-values are not provided because they are statistically invalid after
#' regularized variable selection. Report cross-validated performance metrics and
#' coefficient magnitudes instead.
#'
#' @param x A \code{data.frame}, \code{tibble}, or \code{matrix} of predictor
#'   variables. Must have column names (special characters are supported). All
#'   columns are treated as features to be penalized unless overridden by
#'   \code{not_penalized}. Rows must correspond to the same samples as
#'   \code{metadata}.
#' @param metadata A \code{data.frame} or \code{tibble} containing sample-level
#'   metadata. Must have the same number of rows as \code{x} and in the same
#'   order. Must contain the column specified by \code{pred}.
#' @param pred A single character string naming a column in \code{metadata} to
#'   use as the dependent variable. If the column is numeric, Gaussian regression
#'   is performed. If it is character or factor, classification is performed.
#' @param ref A single character string specifying the reference level for
#'   categorical \code{pred}. Required when \code{pred} is categorical; ignored
#'   otherwise. Must be a value present in the \code{pred} column.
#' @param not_penalized Optional character vector of column names from
#'   \code{metadata} and/or \code{x} to include in the model without
#'   penalization (penalty factor = 0). Useful for known confounders or
#'   covariates that must be controlled. The column named in \code{pred} is
#'   automatically excluded with a warning if accidentally included (e.g.
#'   when passing \code{names(metadata)}). Default: \code{NULL}.
#' @param is_numeric Optional character vector of column names that are listed
#'   in \code{not_penalized} and should be treated as already numeric, bypassing
#'   the automatic non-numeric detection and integer-code conversion. Use this
#'   when a column is genuinely numeric but is stored as \code{character} or
#'   \code{factor} (e.g. due to upstream data import). All names must exist in
#'   \code{metadata} or \code{x} and must also appear in \code{not_penalized}.
#'   Default: \code{NULL}.
#' @param train_percent A numeric scalar in \code{(0, 1)} specifying the
#'   proportion of samples used for model training. The remaining
#'   \code{1 - train_percent} is held out for evaluation. Default: \code{0.8}.
#' @param alpha A numeric scalar or vector with values in \code{[0, 1]}:
#'   \describe{
#'     \item{0}{Ridge regression (L2 penalty; retains all predictors)}
#'     \item{1}{LASSO regression (L1 penalty; performs variable selection)}
#'     \item{(0, 1)}{Elastic Net (blend of L1 and L2)}
#'   }
#'   When a vector is supplied all values are evaluated and the best model is
#'   selected. Default: \code{0.5}.
#' @param lambda A single string controlling lambda selection after
#'   cross-validation:
#'   \describe{
#'     \item{\code{"1se"}}{Lambda within one standard error of the minimum CV
#'       error (more conservative; fewer selected predictors). Default.}
#'     \item{\code{"min"}}{Lambda minimizing CV error (less conservative; more
#'       selected predictors).}
#'   }
#' @param ref A single character string specifying the reference category for
#'   categorical \code{pred}. Required for classification; ignored for numeric
#'   \code{pred}.
#' @param cv_folds A positive integer between 3 and 20 specifying the number of
#'   cross-validation folds. Default: \code{10}.
#' @param type_measure A character string specifying the CV loss function.
#'   \code{NULL} (default) selects automatically based on \code{pred} type:
#'   \describe{
#'     \item{Categorical}{\code{"class"} (misclassification error). Also accepts
#'       \code{"auc"} (binary only; falls back to \code{"class"} for multinomial)
#'       or \code{"deviance"}.}
#'     \item{Numeric}{\code{"mse"} (mean squared error). Also accepts
#'       \code{"mae"} or \code{"deviance"}.}
#'   }
#' @param parallel Logical. If \code{TRUE}, cross-validation is parallelized
#'   using a registered \pkg{foreach} backend (e.g., via
#'   \code{doParallel::registerDoParallel()}). A warning is issued if
#'   \code{TRUE} but no backend is registered. Default: \code{FALSE}.
#' @param standardize Logical. Whether to standardize predictors inside
#'   \code{glmnet} prior to fitting. Set to \code{FALSE} (default) when data
#'   have already been preprocessed/scaled. Default: \code{FALSE}.
#' @param maxit A positive integer specifying the maximum number of iterations
#'   for coordinate descent convergence. Default: \code{100000}.
#' @param seed A single numeric value passed to \code{set.seed()} for
#'   reproducibility of the train/test split and cross-validation fold
#'   assignment. Set to \code{NULL} to disable. Default: \code{123}.
#'
#' @return A named list of class \code{"run_regreg"} with the following elements:
#' \describe{
#'   \item{\code{modelsummary}}{A \code{data.frame} summarising performance
#'     metrics (and a rank column) for every alpha tested.}
#'   \item{\code{sampledistribution}}{A \code{data.frame} showing the number of
#'     training and testing samples per group (classification) or overall
#'     (regression), with a \code{Total} row appended.}
#'   \item{\code{BestModel}}{A named list for the best-performing alpha:
#'     \describe{
#'       \item{\code{Model}}{The fitted \code{cv.glmnet} object.}
#'       \item{\code{Performance}}{A one-row \code{data.frame} of metrics
#'         (Accuracy + Kappa for classification; MSE, RMSE, MAE, R2 for
#'         regression).}
#'       \item{\code{Coefficients}}{A \code{data.frame} of non-zero coefficients
#'         (intercept included if non-zero). Columns: \code{Feature},
#'         \code{Coefficient}, and \code{OddsRatio} (classification) or just
#'         \code{Feature} and \code{Coefficient} (regression). Multinomial
#'         models include a leading \code{Comparison} column.}
#'       \item{\code{N_Coefficients}}{Integer; total non-zero coefficients
#'         including intercept(s).}
#'       \item{\code{ConfusionMatrix}}{A \code{confusionMatrix} object
#'         (\pkg{caret}); \code{NULL} for regression.}
#'       \item{\code{Predictions}}{A \code{data.frame} with columns
#'         \code{Actual} and \code{Predicted}.}
#'       \item{\code{Lambda}}{Numeric; the selected lambda value.}
#'       \item{\code{Alpha}}{Numeric; the alpha value of this model.}
#'     }
#'   }
#'   \item{\code{AllModels}}{A named list (\code{"alpha_<value>"}) with the
#'     same structure as \code{BestModel} for every alpha tested.}
#'   \item{\code{AlphaComparison}}{A \code{data.frame} comparing all tested
#'     alpha values by performance metric and number of non-zero coefficients.
#'     \code{NULL} when only one alpha is tested.}
#' }
#'
#' @importFrom stats predict
#' 
#' @author John Lennon L. Calorio
#'
#' @references
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for
#' Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical
#' Software}, 33(1), 1--22. \doi{10.18637/jss.v033.i01}.
#'
#' Zou, H. and Hastie, T. (2005) Regularization and variable selection via the
#' elastic net. \emph{Journal of the Royal Statistical Society: Series B},
#' 67(2), 301--320.
#'
#' Tibshirani, R. (1996) Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society: Series B}, 58(1), 267--288.
#'
#' Kuhn, M. (2008) Building predictive models in R using the caret package.
#' \emph{Journal of Statistical Software}. \doi{10.18637/jss.v028.i05}.
#'
#' Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization
#' Paths for Cox's Proportional Hazards Model via Coordinate Descent.
#' \emph{Journal of Statistical Software}, 39(5), 1--13.
#' \doi{10.18637/jss.v039.i05}.
#'
#' @examples
#' \dontrun{
#' ## ---- Simulate data -------------------------------------------------------
#' set.seed(42)
#' n  <- 80
#' p  <- 30
#'
#' feat <- as.data.frame(
#'   matrix(rnorm(n * p), nrow = n,
#'          dimnames = list(NULL, paste0("feat_", seq_len(p))))
#' )
#'
#' meta <- data.frame(
#'   SampleID = paste0("S", seq_len(n)),
#'   Group    = sample(c("Control", "Case"), n, replace = TRUE),
#'   Age      = round(runif(n, 20, 65))
#' )
#'
#' ## ---- Binary classification (Elastic Net, single alpha) -------------------
#' res_bin <- run_regreg(
#'   x             = feat,
#'   metadata      = meta,
#'   pred          = "Group",
#'   ref           = "Control",
#'   alpha         = 0.5,
#'   train_percent = 0.8,
#'   seed          = 123
#' )
#' print(res_bin)
#' res_bin$modelsummary
#' res_bin$BestModel$Coefficients
#'
#' ## ---- Binary classification (multiple alpha sweep) ------------------------
#' res_sweep <- run_regreg(
#'   x             = feat,
#'   metadata      = meta,
#'   pred          = "Group",
#'   ref           = "Control",
#'   alpha         = c(0, 0.25, 0.5, 0.75, 1),
#'   lambda        = "1se",
#'   seed          = 123
#' )
#' res_sweep$AlphaComparison
#'
#' ## ---- LASSO with an unpenalized covariate ---------------------------------
#' res_lasso <- run_regreg(
#'   x             = feat,
#'   metadata      = meta,
#'   pred          = "Group",
#'   ref           = "Control",
#'   not_penalized = "Age",
#'   alpha         = 1,
#'   seed          = 123
#' )
#'
#' ## ---- Gaussian regression (numeric pred) ----------------------------------
#' res_num <- run_regreg(
#'   x             = feat,
#'   metadata      = meta,
#'   pred          = "Age",
#'   alpha         = c(0, 0.5, 1),
#'   train_percent = 0.75,
#'   seed          = 7
#' )
#' res_num$BestModel$Performance
#'
#' ## ---- Multinomial classification ------------------------------------------
#' meta3 <- meta
#' meta3$Group <- sample(c("A", "B", "C"), n, replace = TRUE)
#'
#' res_multi <- run_regreg(
#'   x             = feat,
#'   metadata      = meta3,
#'   pred          = "Group",
#'   ref           = "A",
#'   alpha         = 0.5,
#'   seed          = 1
#' )
#' res_multi$BestModel$Coefficients
#' }
#'
#' @export
run_regreg <- function(
    x,
    metadata,
    pred,
    ref           = NULL,
    not_penalized = NULL,
    is_numeric    = NULL,
    train_percent = 0.8,
    alpha         = 0.5,
    lambda        = "1se",
    cv_folds      = 10L,
    type_measure  = NULL,
    parallel      = FALSE,
    standardize   = FALSE,
    maxit         = 100000L,
    seed          = 123
) {

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  .rr_check_x(x)
  .rr_check_metadata(metadata, x)
  .rr_check_pred(pred, metadata)
  .rr_check_scalar_numeric(train_percent, "train_percent", lo = 0, hi = 1,
                           exclusive = TRUE)
  .rr_check_alpha(alpha)
  .rr_check_lambda(lambda)
  .rr_check_cv_folds(cv_folds)
  .rr_check_type_measure(type_measure)
  .rr_check_flag(parallel,    "parallel")
  .rr_check_flag(standardize, "standardize")
  .rr_check_maxit(maxit)
  .rr_check_seed(seed)
  .rr_check_is_numeric(is_numeric, not_penalized, metadata, x)

  # ---------------------------------------------------------------------------
  # 2. Build feature matrix and response vector
  # ---------------------------------------------------------------------------
  prep <- .rr_prepare(x, metadata, pred, ref, not_penalized, is_numeric)

  feat_mat  <- prep$feat_mat    # numeric matrix; rows = samples (non-NA)
  resp      <- prep$resp        # factor or numeric
  pen_fac   <- prep$pen_fac     # penalty.factor vector
  resp_type <- prep$resp_type   # "categorical" | "numeric"
  resp_info <- prep$resp_info   # list with levels, ref_level, etc.

  # ---------------------------------------------------------------------------
  # 3. Train / test split
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  n_total    <- nrow(feat_mat)
  train_idx  <- .rr_split(resp, resp_type, train_percent)
  x_train    <- feat_mat[train_idx,  , drop = FALSE]
  x_test     <- feat_mat[-train_idx, , drop = FALSE]
  y_train    <- resp[train_idx]
  y_test     <- resp[-train_idx]

  # Guard: non-empty splits
  if (nrow(x_train) == 0L || nrow(x_test) == 0L)
    stop("Train/test split produced an empty set. Adjust 'train_percent'.",
         call. = FALSE)

  # ---------------------------------------------------------------------------
  # 4. Adjust cv_folds if necessary (stratified)
  # ---------------------------------------------------------------------------
  cv_folds <- .rr_adjust_folds(y_train, resp_type, cv_folds)

  # ---------------------------------------------------------------------------
  # 5. Sample distribution summary
  # ---------------------------------------------------------------------------
  sampledist <- .rr_sample_distribution(y_train, y_test, resp_type,
                                        train_percent)

  # ---------------------------------------------------------------------------
  # 6. Fixed fold IDs (generated once; reused across all alpha values)
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  fold_ids <- if (resp_type == "categorical") {
    caret::createFolds(y_train, k = cv_folds, list = FALSE)
  } else {
    n_tr <- nrow(x_train)
    sample(rep(seq_len(cv_folds), length.out = n_tr))
  }

  # ---------------------------------------------------------------------------
  # 7. Parallel backend check
  # ---------------------------------------------------------------------------
  if (isTRUE(parallel) && !foreach::getDoParRegistered()) {
    warning(
      "'parallel = TRUE' but no parallel backend is registered. ",
      "Cross-validation will run sequentially. ",
      "Register a backend first, e.g. doParallel::registerDoParallel().",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 8. Determine family and type_measure
  # ---------------------------------------------------------------------------
  family_type  <- .rr_family(resp_info)
  type_measure <- .rr_type_measure(type_measure, resp_type, family_type)

  # ---------------------------------------------------------------------------
  # 9. Fit models for each alpha
  # ---------------------------------------------------------------------------
  alpha_vals <- unique(alpha)
  all_models <- vector("list", length(alpha_vals))
  names(all_models) <- paste0("alpha_", alpha_vals)

  alpha_perf <- vector("list", length(alpha_vals))

  for (i in seq_along(alpha_vals)) {
    av <- alpha_vals[i]
    fit <- .rr_fit_one(
      x_train      = x_train,
      y_train      = y_train,
      x_test       = x_test,
      y_test       = y_test,
      av           = av,
      lambda_rule  = lambda,
      family_type  = family_type,
      type_measure = type_measure,
      fold_ids     = fold_ids,
      pen_fac      = pen_fac,
      parallel     = parallel,
      standardize  = standardize,
      maxit        = maxit,
      resp_type    = resp_type,
      resp_info    = resp_info
    )
    all_models[[i]] <- fit

    alpha_perf[[i]] <- data.frame(
      Alpha          = av,
      fit$Performance,
      N_Coefficients = fit$N_Coefficients,
      stringsAsFactors = FALSE
    )
  }

  # Combine performance table
  perf_df <- do.call(rbind, alpha_perf)
  rownames(perf_df) <- NULL

  # ---------------------------------------------------------------------------
  # 10. Select best model
  # ---------------------------------------------------------------------------
  best_idx <- if (resp_type == "categorical") {
    which.max(perf_df$Accuracy)
  } else {
    which.min(perf_df$RMSE)
  }
  best_model <- all_models[[best_idx]]

  # ---------------------------------------------------------------------------
  # 11. Model summary with rank
  # ---------------------------------------------------------------------------
  modelsummary <- perf_df
  if (resp_type == "categorical") {
    modelsummary$AccuracyRank <- rank(-modelsummary$Accuracy,
                                     ties.method = "min")
  } else {
    modelsummary$RMSERank <- rank(modelsummary$RMSE, ties.method = "min")
  }
  modelsummary <- modelsummary[order(modelsummary$Alpha), ]
  rownames(modelsummary) <- NULL

  # ---------------------------------------------------------------------------
  # 12. Alpha comparison (only meaningful for >1 alpha)
  # ---------------------------------------------------------------------------
  alpha_comparison <- if (length(alpha_vals) > 1L) {
    perf_df[order(perf_df$Alpha), ]
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # 13. Assemble and return
  # ---------------------------------------------------------------------------
  out <- structure(
    list(
      modelsummary       = modelsummary,
      sampledistribution = sampledist,
      BestModel          = best_model,
      AllModels          = all_models,
      AlphaComparison    = alpha_comparison,
      x_names            = colnames(x)   # stored for get_sigfeatures()
    ),
    class = "run_regreg"
  )
  out
}


# =============================================================================
# Internal helpers
# =============================================================================

# ---- Input checkers ---------------------------------------------------------

.rr_check_x <- function(x) {
  if (!is.data.frame(x) && !is.matrix(x))
    stop("'x' must be a data.frame, tibble, or matrix.", call. = FALSE)
  if (is.null(colnames(x)))
    stop("'x' must have column names.", call. = FALSE)
  if (nrow(x) < 4L)
    stop("'x' must have at least 4 rows.", call. = FALSE)
  if (ncol(x) < 1L)
    stop("'x' must have at least 1 column.", call. = FALSE)
}

.rr_check_metadata <- function(metadata, x) {
  if (!is.data.frame(metadata) && !is.matrix(metadata))
    stop("'metadata' must be a data.frame or tibble.", call. = FALSE)
  if (nrow(metadata) != nrow(x))
    stop(
      "'metadata' must have the same number of rows as 'x'. ",
      "Got nrow(x) = ", nrow(x), " but nrow(metadata) = ", nrow(metadata), ".",
      call. = FALSE
    )
}

.rr_check_pred <- function(pred, metadata) {
  if (!is.character(pred) || length(pred) != 1L || nchar(trimws(pred)) == 0L)
    stop("'pred' must be a single non-empty character string.", call. = FALSE)
  if (!pred %in% colnames(metadata))
    stop(
      "'pred' column \"", pred, "\" not found in 'metadata'. ",
      "Available columns: ", paste(colnames(metadata), collapse = ", "), ".",
      call. = FALSE
    )
}

.rr_check_scalar_numeric <- function(val, nm, lo = -Inf, hi = Inf,
                                     exclusive = FALSE) {
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val))
    stop("'", nm, "' must be a single finite numeric value.", call. = FALSE)
  ok <- if (exclusive) val > lo && val < hi else val >= lo && val <= hi
  if (!ok)
    stop(
      "'", nm, "' must be ", if (exclusive) "strictly " else "",
      "between ", lo, " and ", hi, ". Got: ", val, ".",
      call. = FALSE
    )
}

.rr_check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) < 1L)
    stop("'alpha' must be a numeric scalar or vector.", call. = FALSE)
  if (any(!is.finite(alpha)) || any(alpha < 0) || any(alpha > 1))
    stop("All 'alpha' values must be finite and in [0, 1].", call. = FALSE)
}

.rr_check_lambda <- function(lambda) {
  if (!is.character(lambda) || length(lambda) != 1L ||
      !lambda %in% c("1se", "min"))
    stop("'lambda' must be either \"1se\" or \"min\".", call. = FALSE)
}

.rr_check_cv_folds <- function(cv_folds) {
  if (!is.numeric(cv_folds) || length(cv_folds) != 1L ||
      cv_folds != round(cv_folds) || cv_folds < 3L || cv_folds > 20L)
    stop("'cv_folds' must be an integer between 3 and 20.", call. = FALSE)
}

.rr_check_type_measure <- function(type_measure) {
  if (is.null(type_measure)) return(invisible(NULL))
  valid <- c("class", "auc", "deviance", "mse", "mae")
  if (!is.character(type_measure) || length(type_measure) != 1L ||
      !type_measure %in% valid)
    stop(
      "'type_measure' must be one of: ",
      paste(valid, collapse = ", "), ", or NULL.",
      call. = FALSE
    )
}

.rr_check_flag <- function(val, nm) {
  if (!is.logical(val) || length(val) != 1L || is.na(val))
    stop("'", nm, "' must be TRUE or FALSE.", call. = FALSE)
}

.rr_check_maxit <- function(maxit) {
  if (!is.numeric(maxit) || length(maxit) != 1L ||
      maxit != round(maxit) || maxit < 1L)
    stop("'maxit' must be a single positive integer.", call. = FALSE)
}

.rr_check_seed <- function(seed) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L))
    stop("'seed' must be a single numeric value or NULL.", call. = FALSE)
}

# ---- Data preparation -------------------------------------------------------

.rr_check_is_numeric <- function(is_numeric, not_penalized, metadata, x) {
  if (is.null(is_numeric)) return(invisible(NULL))
  if (!is.character(is_numeric))
    stop("'is_numeric' must be a character vector of column names, or NULL.",
         call. = FALSE)
  if (is.null(not_penalized))
    stop(
      "'is_numeric' was supplied but 'not_penalized' is NULL. ",
      "'is_numeric' only applies to columns listed in 'not_penalized'.",
      call. = FALSE
    )
  all_cols   <- c(colnames(metadata), colnames(x))
  bad        <- setdiff(is_numeric, all_cols)
  if (length(bad) > 0L)
    stop(
      "The following column(s) in 'is_numeric' were not found in 'metadata' ",
      "or 'x': ", paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  not_in_np  <- setdiff(is_numeric, not_penalized)
  if (length(not_in_np) > 0L)
    warning(
      "The following column(s) in 'is_numeric' are not listed in ",
      "'not_penalized' and will be ignored: ",
      paste(not_in_np, collapse = ", "), ".",
      call. = FALSE
    )
}

.rr_prepare <- function(x, metadata, pred, ref, not_penalized,
                        is_numeric = NULL) {

  # Convert x to matrix (handles special-character column names safely)
  feat_mat <- if (is.matrix(x)) x else as.matrix(x)
  storage.mode(feat_mat) <- "double"

  # Response variable
  resp_raw <- metadata[[pred]]

  # Remove NA in response
  na_mask <- is.na(resp_raw)
  if (any(na_mask)) {
    n_na <- sum(na_mask)
    warning(n_na, " sample(s) with NA in '", pred, "' removed.",
            call. = FALSE)
    feat_mat <- feat_mat[!na_mask, , drop = FALSE]
    metadata <- metadata[!na_mask, , drop = FALSE]
    resp_raw <- resp_raw[!na_mask]
  }

  # Determine response type
  resp_type <- if (is.numeric(resp_raw)) "numeric" else "categorical"

  resp_info <- list(type = resp_type)

  if (resp_type == "categorical") {

    # Validate ref
    uniq_vals <- unique(as.character(resp_raw))
    if (is.null(ref))
      stop(
        "'ref' must be specified when 'pred' is categorical. ",
        "Unique values in '", pred, "': ",
        paste(sort(uniq_vals), collapse = ", "), ".",
        call. = FALSE
      )
    if (!is.character(ref) || length(ref) != 1L)
      stop("'ref' must be a single character string.", call. = FALSE)
    if (!ref %in% uniq_vals)
      stop(
        "'ref' value \"", ref, "\" not found in '", pred, "'. ",
        "Available values: ", paste(sort(uniq_vals), collapse = ", "), ".",
        call. = FALSE
      )

    other_lvls <- sort(setdiff(uniq_vals, ref))
    lvls       <- c(ref, other_lvls)
    resp       <- factor(resp_raw, levels = lvls)

    group_counts <- table(resp)
    if (any(group_counts < 3L))
      warning(
        "One or more groups have fewer than 3 samples: ",
        paste(names(group_counts)[group_counts < 3], collapse = ", "),
        ". Model performance may be unreliable.",
        call. = FALSE
      )

    resp_info$levels      <- lvls
    resp_info$ref_level   <- ref
    resp_info$group_counts <- as.integer(group_counts)
    names(resp_info$group_counts) <- names(group_counts)

  } else {

    resp <- resp_raw

    resp_info$min  <- min(resp, na.rm = TRUE)
    resp_info$max  <- max(resp, na.rm = TRUE)
    resp_info$mean <- mean(resp, na.rm = TRUE)
    resp_info$sd   <- sd(resp, na.rm = TRUE)

  }

  # Remove constant features (zero variance)
  feat_var    <- apply(feat_mat, 2L, var, na.rm = TRUE)
  drop_const  <- is.na(feat_var) | feat_var < .Machine$double.eps
  if (any(drop_const)) {
    n_drop <- sum(drop_const)
    warning(n_drop, " constant or near-constant feature(s) removed from 'x'.",
            call. = FALSE)
    feat_mat <- feat_mat[, !drop_const, drop = FALSE]
  }

  # Handle missing values in features
  if (anyNA(feat_mat)) {
    complete <- complete.cases(feat_mat)
    n_inc    <- sum(!complete)
    warning(n_inc, " sample(s) with missing feature values removed.",
            call. = FALSE)
    feat_mat <- feat_mat[complete, , drop = FALSE]
    resp     <- if (resp_type == "categorical") resp[complete] else resp[complete]
    metadata <- metadata[complete, , drop = FALSE]
  }

  if (nrow(feat_mat) < 4L)
    stop("Fewer than 4 samples remain after removing NAs. ",
         "Check 'x' and 'metadata' for missing values.", call. = FALSE)

  # penalty.factor: 1 (penalized) for all x columns by default
  pen_fac <- rep(1.0, ncol(feat_mat))

  # Append and zero-penalize not_penalized columns
  if (!is.null(not_penalized)) {
    if (!is.character(not_penalized))
      stop("'not_penalized' must be a character vector of column names.",
           call. = FALSE)

    # --- Guard 1: remove 'pred' if it snuck in (e.g. names(metadata)) -------
    if (pred %in% not_penalized) {
      warning(
        "The dependent variable '", pred, "' was found in 'not_penalized' ",
        "and has been automatically removed. The response variable must not ",
        "be used as a predictor.",
        call. = FALSE
      )
      not_penalized <- setdiff(not_penalized, pred)
    }

    if (length(not_penalized) == 0L) {
      warning(
        "After removing 'pred' from 'not_penalized', no columns remain. ",
        "'not_penalized' will be ignored.",
        call. = FALSE
      )
    } else {

      # Resolve from metadata and/or x
      from_meta  <- intersect(not_penalized, colnames(metadata))
      from_x     <- intersect(not_penalized, colnames(x))
      missing_np <- setdiff(not_penalized, c(colnames(metadata), colnames(x)))

      if (length(missing_np) > 0L)
        stop(
          "The following 'not_penalized' column(s) were not found in ",
          "'metadata' or 'x': ", paste(missing_np, collapse = ", "), ".",
          call. = FALSE
        )

      # Columns already in feat_mat (from x) — just flip their penalty to 0
      if (length(from_x) > 0L) {
        idx_in_feat          <- match(from_x, colnames(feat_mat))
        pen_fac[idx_in_feat] <- 0.0
      }

      # Columns from metadata — append to feat_mat with penalty 0
      if (length(from_meta) > 0L) {

        # Apply row subsetting that was applied to feat_mat
        np_sub <- metadata[, from_meta, drop = FALSE]

        np_mat <- matrix(
          NA_real_,
          nrow     = nrow(np_sub),
          ncol     = ncol(np_sub),
          dimnames = list(NULL, colnames(np_sub))
        )

        for (j in seq_len(ncol(np_sub))) {
          col_nm <- colnames(np_sub)[j]
          col_j  <- np_sub[[j]]

          # Determine whether to treat as numeric
          force_numeric <- !is.null(is_numeric) && col_nm %in% is_numeric

          if (force_numeric) {
            # User asserts this column is numeric — coerce directly
            coerced <- suppressWarnings(as.numeric(col_j))
            if (anyNA(coerced) && !anyNA(col_j))
              warning(
                "Column '", col_nm, "' was listed in 'is_numeric' but ",
                "could not be fully coerced to numeric; NAs introduced. ",
                "Check the column values.",
                call. = FALSE
              )
            np_mat[, j] <- coerced

          } else if (is.numeric(col_j)) {
            np_mat[, j] <- col_j

          } else {
            # Non-numeric and not overridden: encode as integer codes
            np_mat[, j] <- as.numeric(as.factor(col_j))
            warning(
              "Column '", col_nm, "' in 'not_penalized' is not numeric; ",
              "converted to integer codes for glmnet. ",
              "If this column is already numeric but stored as character or ",
              "factor, add it to the 'is_numeric' argument to bypass this ",
              "conversion.",
              call. = FALSE
            )
          }
        }

        # --- Guard 2: check for NAs in the appended metadata columns ---------
        na_rows <- which(!complete.cases(np_mat))
        if (length(na_rows) > 0L) {
          warning(
            length(na_rows), " sample(s) have missing values in the ",
            "'not_penalized' metadata column(s) and will be removed: ",
            paste(colnames(np_mat)[apply(is.na(np_mat), 2L, any)],
                  collapse = ", "), ".",
            call. = FALSE
          )
          feat_mat <- feat_mat[-na_rows, , drop = FALSE]
          np_mat   <- np_mat[-na_rows,   , drop = FALSE]
          resp     <- if (is.factor(resp)) resp[-na_rows] else resp[-na_rows]
        }

        feat_mat <- cbind(feat_mat, np_mat)
        pen_fac  <- c(pen_fac, rep(0.0, ncol(np_mat)))
      }
    }
  }

  list(
    feat_mat  = feat_mat,
    resp      = resp,
    pen_fac   = pen_fac,
    resp_type = resp_type,
    resp_info = resp_info
  )
}

# ---- Train / test split -----------------------------------------------------

.rr_split <- function(resp, resp_type, train_percent) {
  if (resp_type == "categorical") {
    idx <- caret::createDataPartition(resp, p = train_percent,
                                      list = FALSE, times = 1L)
    as.integer(idx)
  } else {
    n <- length(resp)
    sample(n, size = floor(train_percent * n))
  }
}

# ---- Adjust cv_folds --------------------------------------------------------

.rr_adjust_folds <- function(y_train, resp_type, cv_folds) {
  if (resp_type != "categorical") return(as.integer(cv_folds))
  min_grp <- min(table(y_train))
  if (cv_folds > min_grp) {
    new_folds <- max(3L, as.integer(min_grp))
    warning(
      "'cv_folds' reduced from ", cv_folds, " to ", new_folds,
      " because the smallest training group has only ", min_grp, " sample(s).",
      call. = FALSE
    )
    return(new_folds)
  }
  as.integer(cv_folds)
}

# ---- Resolve glmnet family --------------------------------------------------

.rr_family <- function(resp_info) {
  if (resp_info$type == "numeric") return("gaussian")
  if (length(resp_info$levels) == 2L) "binomial" else "multinomial"
}

# ---- Resolve type_measure ---------------------------------------------------

.rr_type_measure <- function(type_measure, resp_type, family_type) {
  if (is.null(type_measure)) {
    return(if (resp_type == "categorical") "class" else "mse")
  }
  if (resp_type == "categorical") {
    valid <- c("class", "auc", "deviance")
    if (!type_measure %in% valid) {
      warning(
        "'type_measure = \"", type_measure, "\"' is not valid for categorical ",
        "responses. Falling back to \"class\". ",
        "Valid options: ", paste(valid, collapse = ", "), ".",
        call. = FALSE
      )
      return("class")
    }
    if (type_measure == "auc" && family_type == "multinomial") {
      warning(
        "'type_measure = \"auc\"' requires binary classification. ",
        "Falling back to \"class\" for multinomial response.",
        call. = FALSE
      )
      return("class")
    }
  } else {
    valid <- c("mse", "mae", "deviance")
    if (!type_measure %in% valid) {
      warning(
        "'type_measure = \"", type_measure, "\"' is not valid for numeric ",
        "responses. Falling back to \"mse\". ",
        "Valid options: ", paste(valid, collapse = ", "), ".",
        call. = FALSE
      )
      return("mse")
    }
  }
  type_measure
}

# ---- Sample distribution summary --------------------------------------------

.rr_sample_distribution <- function(y_train, y_test, resp_type,
                                     train_percent) {
  pct_train <- round(train_percent * 100)
  pct_test  <- 100L - pct_train
  col_tr    <- paste0("Training (", pct_train, "%)")
  col_te    <- paste0("Testing (",  pct_test,  "%)")

  if (resp_type == "categorical") {
    lvls       <- levels(y_train)
    tr_counts  <- table(factor(y_train, levels = lvls))
    te_counts  <- table(factor(y_test,  levels = lvls))

    body <- data.frame(
      Group         = lvls,
      Training      = as.integer(tr_counts),
      Testing       = as.integer(te_counts),
      Total         = as.integer(tr_counts) + as.integer(te_counts),
      stringsAsFactors = FALSE,
      check.names   = FALSE
    )
    names(body)[2L:3L] <- c(col_tr, col_te)

    total_row <- data.frame(
      Group    = "Total",
      Training = sum(body[[col_tr]]),
      Testing  = sum(body[[col_te]]),
      Total    = sum(body$Total),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(total_row)[2L:3L] <- c(col_tr, col_te)

    rbind(body, total_row)

  } else {
    data.frame(
      Group         = "Numeric response",
      Training      = length(y_train),
      Testing       = length(y_test),
      Total         = length(y_train) + length(y_test),
      stringsAsFactors = FALSE,
      check.names   = FALSE
    ) |> setNames(c("Group", col_tr, col_te, "Total"))
  }
}

# ---- Fit one alpha ----------------------------------------------------------

.rr_fit_one <- function(x_train, y_train, x_test, y_test,
                         av, lambda_rule, family_type, type_measure,
                         fold_ids, pen_fac, parallel, standardize, maxit,
                         resp_type, resp_info) {

  cv_fit <- glmnet::cv.glmnet(
    x              = x_train,
    y              = y_train,
    alpha          = av,
    family         = family_type,
    type.measure   = type_measure,
    foldid         = fold_ids,
    penalty.factor = pen_fac,
    parallel       = parallel,
    standardize    = standardize,
    maxit          = as.integer(maxit)
  )

  # Select lambda
  sel_lambda <- if (lambda_rule == "1se") {
    lam <- cv_fit$lambda.1se
    if (is.na(lam)) {
      warning(
        "lambda.1se is NA for alpha = ", av, ". Falling back to lambda.min.",
        call. = FALSE
      )
      cv_fit$lambda.min
    } else {
      lam
    }
  } else {
    cv_fit$lambda.min
  }

  # Predictions and performance
  if (resp_type == "categorical") {

    pred_class <- as.vector(
      predict(cv_fit, newx = x_test, s = sel_lambda, type = "class")
    )
    pred_class <- factor(pred_class, levels = resp_info$levels)
    y_test_f   <- factor(y_test,     levels = resp_info$levels)

    cm        <- caret::confusionMatrix(pred_class, y_test_f)
    accuracy  <- unname(cm$overall["Accuracy"])
    kappa     <- unname(cm$overall["Kappa"])
    perf      <- data.frame(Accuracy = accuracy, Kappa = kappa,
                             stringsAsFactors = FALSE)

    pred_df <- data.frame(Actual    = as.character(y_test),
                          Predicted = as.character(pred_class),
                          stringsAsFactors = FALSE)
    conf_mat <- cm

  } else {

    pred_num <- as.vector(
      predict(cv_fit, newx = x_test, s = sel_lambda, type = "response")
    )
    residuals <- y_test - pred_num
    mse_val  <- mean(residuals^2)
    rmse_val <- sqrt(mse_val)
    mae_val  <- mean(abs(residuals))
    sst      <- sum((y_test - mean(y_test))^2)
    sse      <- sum(residuals^2)
    r2       <- 1 - sse / sst

    perf <- data.frame(MSE = mse_val, RMSE = rmse_val,
                        MAE = mae_val, R2   = r2,
                        stringsAsFactors = FALSE)

    pred_df  <- data.frame(Actual = y_test, Predicted = pred_num,
                            stringsAsFactors = FALSE)
    conf_mat <- NULL
  }

  # Coefficients
  coef_out <- .rr_extract_coefs(cv_fit, sel_lambda, family_type, resp_info)

  list(
    Model           = cv_fit,
    Performance     = perf,
    ConfusionMatrix = conf_mat,
    Predictions     = pred_df,
    Coefficients    = coef_out$coef_df,
    N_Coefficients  = coef_out$n_coefs,
    Lambda          = sel_lambda,
    Alpha           = av
  )
}

# ---- Extract coefficients ---------------------------------------------------

.rr_extract_coefs <- function(cv_fit, sel_lambda, family_type, resp_info) {

  raw <- glmnet::coef.glmnet(cv_fit, s = sel_lambda)

  if (family_type == "multinomial") {

    lvls      <- resp_info$levels
    ref_level <- resp_info$ref_level

    # Coerce each sparse matrix to a dense column, bind
    coef_mat           <- do.call(cbind, lapply(raw, as.matrix))
    colnames(coef_mat) <- lvls

    # Normalize relative to reference
    ref_vec   <- coef_mat[, ref_level, drop = TRUE]
    norm_mat  <- coef_mat - ref_vec  # broadcasts over columns

    # Build comparison data frames for each non-reference group
    target_lvls <- setdiff(lvls, ref_level)
    rows <- lapply(target_lvls, function(grp) {
      betas     <- norm_mat[, grp]
      active    <- which(betas != 0)
      if (length(active) == 0L) return(NULL)
      data.frame(
        Comparison       = paste0(grp, " vs ", ref_level),
        Feature          = rownames(norm_mat)[active],
        Coefficient      = betas[active],
        OddsRatio        = exp(betas[active]),
        stringsAsFactors = FALSE
      )
    })
    coef_df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
    rownames(coef_df) <- NULL

    return(list(coef_df = coef_df, n_coefs = nrow(coef_df)))
  }

  # Binomial or Gaussian
  coef_df           <- as.data.frame(as.matrix(raw))
  colnames(coef_df) <- "Coefficient"
  coef_df$Feature   <- rownames(coef_df)
  coef_df           <- coef_df[coef_df$Coefficient != 0, , drop = FALSE]

  if (family_type == "binomial") {
    coef_df$OddsRatio <- exp(coef_df$Coefficient)
    coef_df <- coef_df[, c("Feature", "Coefficient", "OddsRatio"),
                        drop = FALSE]
  } else {
    coef_df <- coef_df[, c("Feature", "Coefficient"), drop = FALSE]
  }
  rownames(coef_df) <- NULL

  list(coef_df = coef_df, n_coefs = nrow(coef_df))
}


# =============================================================================
# S3 methods
# =============================================================================

#' @export
print.run_regreg <- function(x, ...) {
  bm   <- x$BestModel
  rtype <- if ("Accuracy" %in% names(bm$Performance)) "categorical" else "numeric"

  cat("=== run_regreg: Regularized Regression ===\n")
  cat("Response type  :", rtype, "\n")
  cat("Best alpha     :", bm$Alpha, "\n")
  cat("Selected lambda:", signif(bm$Lambda, 4L), "\n")
  cat("Non-zero coefs :", bm$N_Coefficients, "\n")
  if (rtype == "categorical") {
    cat("Test Accuracy  :", round(bm$Performance$Accuracy, 4L), "\n")
    cat("Test Kappa     :", round(bm$Performance$Kappa,    4L), "\n")
  } else {
    cat("Test RMSE      :", round(bm$Performance$RMSE, 4L), "\n")
    cat("Test R2        :", round(bm$Performance$R2,   4L), "\n")
  }
  if (!is.null(x$AlphaComparison)) {
    cat("\nAlpha comparison:\n")
    print(x$AlphaComparison, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.run_regreg <- function(object, ...) {
  cat("=== run_regreg Model Summary ===\n")
  print(object$modelsummary, row.names = FALSE)
  cat("\nSample distribution:\n")
  print(object$sampledistribution, row.names = FALSE)
  invisible(object)
}