#' Extract Significant Features from a pondeR Model Result
#'
#' @title Get Features with Non-Zero Coefficients from a Regularized Regression
#'   Result
#'
#' @description
#' Extracts the names of predictor features that received non-zero coefficients
#' from a fitted regularized regression model produced by \code{\link{run_regreg}}.
#' By default the best-performing model (as determined internally by
#' \code{run_regreg}) is used. A specific alpha value may be requested instead.
#'
#' Intercept terms are excluded from the returned vector unless
#' \code{include_intercept = TRUE}. For multinomial models, where the same
#' feature may appear in multiple group comparisons, the union of all non-zero
#' features across comparisons is returned by default; set
#' \code{by_comparison = TRUE} to obtain a named list broken down by comparison.
#'
#' @param res A list of class \code{"run_regreg"} produced by
#'   \code{\link{run_regreg}}.
#' @param alpha A single numeric value in \code{[0, 1]} specifying which
#'   fitted model to query. Must match one of the alpha values that was passed
#'   to \code{run_regreg}. When \code{NULL} (default), the best-performing
#'   model selected by \code{run_regreg} is used.
#' @param include_intercept Logical. If \code{TRUE}, intercept term(s)
#'   (\code{"(Intercept)"}) are retained in the output. Default: \code{FALSE}.
#' @param by_comparison Logical. Only relevant for multinomial classification
#'   results, which contain a \code{Comparison} column in their coefficient
#'   table. When \code{TRUE}, a named list is returned where each element
#'   contains the non-zero feature names for one group comparison (e.g.,
#'   \code{"B vs A"}). When \code{FALSE} (default), the union of non-zero
#'   features across all comparisons is returned as a single character vector.
#'   Ignored for binary classification and regression results.
#' @param sort_by A single string controlling how the returned feature names are
#'   ordered. One of:
#'   \describe{
#'     \item{\code{"name"}}{Alphabetical order (default).}
#'     \item{\code{"coef"}}{Descending absolute coefficient magnitude. Not
#'       available when \code{by_comparison = TRUE} and the result is
#'       multinomial (the union step discards magnitudes); a warning is issued
#'       and the function falls back to \code{"name"} in that case.}
#'     \item{\code{"none"}}{Preserve the original row order from the
#'       coefficient table.}
#'   }
#'
#' @return
#' \describe{
#'   \item{Default (\code{by_comparison = FALSE})}{A character vector of
#'     feature names with non-zero coefficients, in the order specified by
#'     \code{sort_by}. Returns \code{character(0)} with a message when no
#'     non-zero features are found (e.g., the null model).}
#'   \item{\code{by_comparison = TRUE} (multinomial only)}{A named list of
#'     character vectors, one per group comparison. Each element is named by
#'     its comparison label (e.g., \code{"B vs A"}).}
#' }
#'
#' @examples
#' \dontrun{
#' ## ---- Setup ---------------------------------------------------------------
#' set.seed(42)
#' n  <- 80
#' p  <- 30
#'
#' feat <- as.data.frame(
#'   matrix(rnorm(n * p), nrow = n,
#'          dimnames = list(NULL, paste0("feat_", seq_len(p))))
#' )
#' meta <- data.frame(
#'   Group = sample(c("Control", "Case"), n, replace = TRUE)
#' )
#'
#' ## ---- Binary classification -----------------------------------------------
#' res_bin <- run_regreg(
#'   x    = feat, metadata = meta,
#'   pred = "Group", ref = "Control",
#'   alpha = c(0, 0.5, 1), seed = 123
#' )
#'
#' # Best model features (default)
#' get_sigfeatures(res_bin)
#'
#' # Features from LASSO specifically
#' get_sigfeatures(res_bin, alpha = 1)
#'
#' # Sorted by absolute coefficient magnitude
#' get_sigfeatures(res_bin, sort_by = "coef")
#'
#' # Include the intercept term
#' get_sigfeatures(res_bin, include_intercept = TRUE)
#'
#' ## ---- Multinomial classification ------------------------------------------
#' meta3 <- data.frame(
#'   Group = sample(c("A", "B", "C"), n, replace = TRUE)
#' )
#' res_multi <- run_regreg(
#'   x    = feat, metadata = meta3,
#'   pred = "Group", ref = "A",
#'   alpha = 0.5, seed = 1
#' )
#'
#' # Union of non-zero features across all comparisons
#' get_sigfeatures(res_multi)
#'
#' # Per-comparison breakdown
#' get_sigfeatures(res_multi, by_comparison = TRUE)
#'
#' ## ---- Gaussian regression -------------------------------------------------
#' meta_num <- data.frame(Score = rnorm(n))
#' res_num  <- run_regreg(
#'   x    = feat, metadata = meta_num,
#'   pred = "Score", alpha = 0.5, seed = 7
#' )
#' get_sigfeatures(res_num)
#' }
#'
#' @param x_names An optional character vector of column names from the
#'   original \code{x} that was passed to \code{run_regreg()}. When supplied,
#'   only features whose names appear in \code{x_names} are returned, silently
#'   dropping any unpenalized metadata columns that were appended via
#'   \code{not_penalized}. If \code{NULL} (default), the function first checks
#'   whether \code{res} stores the original \code{x} column names internally
#'   (via \code{res$x_names}); if found, those are used automatically. If
#'   neither is available, no name-based filtering is applied and a note is
#'   issued. The simplest usage is \code{x_names = colnames(x)}.
#'
#' @seealso \code{\link{run_regreg}}
#'
#' @export
get_sigfeatures <- function(
    res,
    alpha             = NULL,
    include_intercept = FALSE,
    by_comparison     = FALSE,
    sort_by           = "name",
    x_names           = NULL
) {

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  .gsf_check_res(res)

  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha < 0 || alpha > 1)
      stop(
        "'alpha' must be a single numeric value in [0, 1], or NULL to use ",
        "the best model.",
        call. = FALSE
      )
  }

  if (!is.logical(include_intercept) || length(include_intercept) != 1L ||
      is.na(include_intercept))
    stop("'include_intercept' must be TRUE or FALSE.", call. = FALSE)

  if (!is.logical(by_comparison) || length(by_comparison) != 1L ||
      is.na(by_comparison))
    stop("'by_comparison' must be TRUE or FALSE.", call. = FALSE)

  valid_sort <- c("name", "coef", "none")
  if (!is.character(sort_by) || length(sort_by) != 1L ||
      !sort_by %in% valid_sort)
    stop(
      "'sort_by' must be one of: ", paste(valid_sort, collapse = ", "), ".",
      call. = FALSE
    )

  if (!is.null(x_names) && (!is.character(x_names) || length(x_names) == 0L))
    stop(
      "'x_names' must be a non-empty character vector of column names from ",
      "the original 'x', or NULL.",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 2. Resolve the allowlist of x column names
  # ---------------------------------------------------------------------------
  # Priority: explicit x_names > stored x_names inside res > no filtering
  allowed <- .gsf_resolve_x_names(res, x_names)

  # ---------------------------------------------------------------------------
  # 3. Pull the requested coefficient table
  # ---------------------------------------------------------------------------
  coef_df <- .gsf_get_coef_df(res, alpha)

  # ---------------------------------------------------------------------------
  # 4. Optionally remove intercept rows
  # ---------------------------------------------------------------------------
  if (!include_intercept) {
    coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
  }

  # ---------------------------------------------------------------------------
  # 5. Filter to x-only features (drop metadata / not_penalized columns)
  # ---------------------------------------------------------------------------
  if (!is.null(allowed)) {
    n_before <- nrow(coef_df)
    coef_df  <- coef_df[coef_df$Feature %in% allowed, , drop = FALSE]
    n_dropped <- n_before - nrow(coef_df)
    if (n_dropped > 0L)
      message(
        n_dropped, " non-zero coefficient(s) belonging to 'not_penalized' ",
        "metadata column(s) were excluded from the output. ",
        "They are still part of the fitted model."
      )
  }

  # Guard: nothing left
  if (nrow(coef_df) == 0L) {
    message(
      "No non-zero features from 'x' found for the selected model",
      if (!is.null(alpha)) paste0(" (alpha = ", alpha, ")") else " (best model)",
      ". The model may have selected no predictors at the chosen lambda. ",
      "Consider using lambda = \"min\" or a different alpha."
    )
    return(character(0L))
  }

  # ---------------------------------------------------------------------------
  # 6. Is this a multinomial result? (has a Comparison column)
  # ---------------------------------------------------------------------------
  is_multinomial <- "Comparison" %in% colnames(coef_df)

  # ---------------------------------------------------------------------------
  # 7. by_comparison = TRUE  (multinomial only)
  # ---------------------------------------------------------------------------
  if (by_comparison) {

    if (!is_multinomial) {
      warning(
        "'by_comparison = TRUE' is only meaningful for multinomial results. ",
        "This model has no 'Comparison' column — returning a plain vector.",
        call. = FALSE
      )
      # fall through to the plain-vector path below
    } else {
      comps <- unique(coef_df$Comparison)
      out   <- lapply(comps, function(cmp) {
        sub_df <- coef_df[coef_df$Comparison == cmp, , drop = FALSE]
        .gsf_sort(sub_df, sort_by, context = paste0("comparison '", cmp, "'"))
      })
      names(out) <- comps
      return(out)
    }
  }

  # ---------------------------------------------------------------------------
  # 8. Plain vector path
  # ---------------------------------------------------------------------------

  # For multinomial: union across comparisons; magnitude sort is ambiguous here
  if (is_multinomial && sort_by == "coef") {
    warning(
      "'sort_by = \"coef\"' is ambiguous for multinomial results when ",
      "'by_comparison = FALSE', because the same feature can have different ",
      "magnitudes in different comparisons. Falling back to sort_by = \"name\".",
      call. = FALSE
    )
    sort_by <- "name"
  }

  .gsf_sort(coef_df, sort_by, context = "selected model")
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.gsf_resolve_x_names <- function(res, x_names) {
  # Explicit argument always wins
  if (!is.null(x_names)) return(x_names)
  # Fall back to names stored inside the result object by run_regreg
  if (!is.null(res[["x_names"]])) return(res[["x_names"]])
  # No allowlist available — cannot filter
  NULL
}

#' @keywords internal
.gsf_check_res <- function(res) {

  if (!is.list(res))
    stop(
      "'res' must be a list returned by a pondeR modelling function ",
      "(e.g., run_regreg()).",
      call. = FALSE
    )

  # Class-based dispatch: currently only run_regreg is supported
  if (!inherits(res, "run_regreg"))
    stop(
      "'res' does not appear to be the output of a supported pondeR function. ",
      "Currently supported: run_regreg(). ",
      "Make sure you pass the result object directly, e.g. get_sigfeatures(res).",
      call. = FALSE
    )

  # Structural sanity checks
  if (!"BestModel" %in% names(res))
    stop(
      "'res' is missing the 'BestModel' element. ",
      "The object may be incomplete or corrupted.",
      call. = FALSE
    )

  if (!"Coefficients" %in% names(res$BestModel))
    stop(
      "'res$BestModel' is missing the 'Coefficients' element. ",
      "The model may have failed internally.",
      call. = FALSE
    )
}

#' @keywords internal
.gsf_get_coef_df <- function(res, alpha) {

  if (is.null(alpha)) {
    return(res$BestModel$Coefficients)
  }

  # Locate requested alpha in AllModels
  if (!"AllModels" %in% names(res) || is.null(res$AllModels)) {
    # Single-alpha run; check if the only alpha matches
    fitted_alpha <- res$BestModel$Alpha
    if (!isTRUE(all.equal(alpha, fitted_alpha)))
      stop(
        "Alpha ", alpha, " was not found. ",
        "This result contains only one fitted alpha: ", fitted_alpha, ". ",
        "Use alpha = NULL to select it, or re-run run_regreg() with the ",
        "desired alpha values.",
        call. = FALSE
      )
    return(res$BestModel$Coefficients)
  }

  key <- paste0("alpha_", alpha)

  if (!key %in% names(res$AllModels)) {
    available <- sort(as.numeric(gsub("alpha_", "", names(res$AllModels))))
    stop(
      "Alpha ", alpha, " was not found among the fitted models. ",
      "Available alpha values: ", paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }

  coef_df <- res$AllModels[[key]]$Coefficients

  if (is.null(coef_df))
    stop(
      "The coefficient table for alpha = ", alpha, " is NULL. ",
      "The model at this alpha may have failed internally.",
      call. = FALSE
    )

  coef_df
}

#' @keywords internal
.gsf_sort <- function(coef_df, sort_by, context = "") {

  # Deduplicate feature names (union across comparisons for multinomial)
  if ("Comparison" %in% colnames(coef_df)) {
    # For multinomial union path: just unique feature names
    feats <- unique(coef_df$Feature)

    if (sort_by == "name") {
      return(sort(feats))
    }
    # "none" or "coef" (coef already warned above at the call site)
    return(feats)
  }

  # Binomial / Gaussian: Feature + Coefficient columns present
  if (sort_by == "name") {
    return(sort(coef_df$Feature))
  }

  if (sort_by == "coef") {
    if (!"Coefficient" %in% colnames(coef_df)) {
      warning(
        "Cannot sort by 'coef': no 'Coefficient' column found for ", context,
        ". Falling back to sort_by = \"name\".",
        call. = FALSE
      )
      return(sort(coef_df$Feature))
    }
    ord <- order(abs(coef_df$Coefficient), decreasing = TRUE)
    return(coef_df$Feature[ord])
  }

  # sort_by == "none"
  coef_df$Feature
}