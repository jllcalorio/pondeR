#' @title Perform Binary Logistic Regression (Standard or Firth-Corrected) with Iterative Feature Combinations
#'
#' @description
#' A robust wrapper for \code{\link[stats]{glm}} (standard logistic regression) and
#' \code{\link[logistf]{logistf}} (Firth's bias-reduced penalized-likelihood logistic regression)
#' that performs binary logistic regression. It automatically handles column names with spaces,
#' computes advanced model fit metrics, maps specific reference categories, removes specified
#' categories from predictor variables, and optionally iterates through all possible combinations
#' of independent variables. When \code{apply_firth = TRUE}, Firth's correction is applied
#' selectively — only to models where at least one predictor-outcome cross-tabulation contains a
#' zero cell, which is the primary condition warranting penalized likelihood. A column in the
#' output indicates whether Firth's correction was applied to each model. Output is formatted as
#' a clean data frame.
#'
#' When \code{crude = TRUE}, a separate simple logistic regression is fitted for each variable
#' in \code{indep} against \code{y}, producing unadjusted (crude) odds ratios. Results are
#' stacked into a single tidy data frame. The \code{confounders} and \code{iterate} arguments
#' are ignored in this mode.
#'
#' @param x A data frame, matrix, or tibble containing the dataset.
#' @param y A string specifying the valid, categorical column name in \code{x} to be used as the dependent variable.
#' @param confounders A character vector of valid column names in \code{x} to be used as confounding factors.
#'   Ignored when \code{crude = TRUE} (a warning is emitted). Default is \code{NULL}.
#' @param indep A character vector of valid column names in \code{x} to be used as independent variables.
#'   When \code{crude = TRUE}, each variable is modelled individually against \code{y}.
#'   Default is \code{NULL}.
#' @param ref A string specifying the valid category of \code{y} to be used as the reference category.
#' @param ref_levels A list of two-sided formulas mapping categorical predictor column names to their
#'   desired reference categories. The left-hand side must be a valid column name of \code{x} (quoted
#'   as a string), and the right-hand side must be a valid category of that column
#'   (e.g., \code{list("Sex" ~ "Female", "District" ~ "North")}). Applied in both standard and
#'   crude modes. Default is \code{NULL}.
#' @param remove A list of two-sided formulas specifying categories to drop from predictor variables
#'   before analysis. The left-hand side must be a valid column name of \code{x} (quoted as a string),
#'   and the right-hand side must be a valid category of that column to be removed. Multiple mappings
#'   may target the same column to remove multiple categories
#'   (e.g., \code{list("Sex" ~ "Unknown", "District" ~ "Other")}). Applied in both standard and
#'   crude modes. Default is \code{NULL}.
#' @param exclude A character vector of categories in \code{y} to drop before analysis. Applied in
#'   both standard and crude modes. Default is \code{NULL}.
#' @param iterate Logical; if \code{FALSE} (default), combines \code{confounders} and \code{indep}
#'   into a single model. If \code{TRUE}, fits separate models for all possible combinations of
#'   \code{indep} added to \code{confounders}. Ignored when \code{crude = TRUE}.
#' @param crude Logical; if \code{TRUE}, fits a separate simple (unadjusted) logistic regression
#'   for each variable in \code{indep} against \code{y}, producing crude odds ratios. Results are
#'   stacked into a single tidy data frame with the same column structure as the standard output.
#'   \code{confounders} and \code{iterate} are ignored in this mode (a warning is emitted if
#'   \code{confounders} is non-\code{NULL}). All other parameters (\code{ref}, \code{ref_levels},
#'   \code{remove}, \code{exclude}, \code{add_vif}, \code{apply_firth},
#'   \code{force_apply_firth}) remain applicable. Default is \code{FALSE}.
#' @param add_vif Logical; if \code{TRUE} (default), adds Variance Inflation Factor (VIF) values
#'   to the output. VIFs are computed via \code{car::vif()}, which returns GVIF for categorical
#'   predictors and standard VIF for continuous predictors, consistent with jamovi's output.
#'   Note that VIF is computed from a standard \code{glm()} fit even when Firth's correction is
#'   applied, as \code{car::vif()} does not support \code{logistf} objects. In crude mode,
#'   VIF is not meaningful for simple (one-predictor) models and will be \code{NA} throughout.
#' @param apply_firth Logical; if \code{TRUE} (default), Firth's bias-reduced penalized-likelihood
#'   logistic regression (\code{logistf::logistf()}) is applied \emph{selectively} — only to models
#'   where at least one predictor variable has a zero cell count in its cross-tabulation with the
#'   outcome. Models with no zero cells are fitted using standard \code{glm()} regardless of this
#'   setting. A logical column \code{Firth Corrected} is added to the output indicating whether
#'   Firth's correction was applied to each model. Applicable in both standard and crude modes.
#' @param force_apply_firth Logical; if \code{TRUE} (default \code{FALSE}), overrides the dynamic
#'   zero-cell detection of \code{apply_firth} and applies Firth's correction to \emph{all} models
#'   regardless of whether zero cells are present. Ignored if \code{apply_firth = FALSE}.
#'   Applicable in both standard and crude modes.
#'
#' @details
#' \strong{Crude odds ratios (\code{crude = TRUE}):} Each variable in \code{indep} is fitted in
#' a separate simple logistic regression model with \code{y} as the only outcome and no
#' covariates, yielding unadjusted (crude) odds ratios. All preprocessing steps — category
#' exclusion (\code{exclude}), predictor category removal (\code{remove}), dependent variable
#' releveling (\code{ref}), and predictor releveling (\code{ref_levels}) — are applied before
#' each model is fitted. The results are stacked row-wise into a single tidy data frame. Each
#' variable's rows are separated by a blank spacer row to visually group the predictors.
#' Inline model-fit metrics (Deviance, AIC, BIC, and R2 measures) are placed in the first row
#' of each variable's block (aligned with its intercept row), and \code{NA} elsewhere. Because
#' each model contains only one predictor, VIF is not meaningful and the \code{VIF} column
#' will be \code{NA} for all rows. \code{confounders} and \code{iterate} are silently ignored
#' after emitting a warning.
#'
#' When \code{iterate = TRUE}, the function generates all possible combinations (the power set)
#' of the \code{indep} vector. For example, if \code{indep = c("A", "B", "C")}, the combinations
#' evaluated alongside \code{confounders} will be: \code{A}, \code{B}, \code{C}, \code{A+B},
#' \code{A+C}, \code{B+C}, and \code{A+B+C}.
#'
#' \strong{Firth's correction (\code{apply_firth}):} Firth's penalized-likelihood method
#' (\code{logistf::logistf()}) addresses the monotone likelihood problem caused by complete or
#' quasi-complete separation, which arises when a predictor level perfectly predicts one outcome
#' group — typically evidenced by a zero cell in the cross-tabulation. When \code{apply_firth = TRUE},
#' the function checks each model's predictors against the outcome for zero cells. If any are found,
#' \code{logistf()} is used in place of \code{glm()} for that model. If \code{force_apply_firth = TRUE},
#' this check is skipped and all models are fitted with \code{logistf()}. The \code{Firth Corrected}
#' column in the output records which models received the correction. Note that Firth-corrected
#' models report profile likelihood confidence intervals rather than Wald intervals, consistent with
#' \code{logistf}'s default behavior.
#'
#' \strong{Category removal (\code{remove}):} Rows where the specified predictor variable equals
#' the given category are dropped from \code{x} before any model is fitted. This is applied after
#' \code{exclude} but before factor releveling. Multiple formulas targeting the same column are
#' each applied in sequence, allowing several categories to be removed from one variable.
#'
#' \strong{Aliased coefficients:} When standard \code{glm()} cannot estimate a coefficient due to
#' perfect or quasi-complete separation, it sets that coefficient to \code{NA}. \code{summary()}
#' silently omits these aliased rows while \code{confint.default()} retains them, causing a
#' row-count mismatch. \code{run_logreg()} detects this automatically, aligns both objects to
#' their shared row names, and emits a warning naming every dropped aliased term. This situation
#' is largely avoided when \code{apply_firth = TRUE} since such models will be routed to
#' \code{logistf()} instead.
#'
#' \strong{VIF computation:} This function uses \code{car::vif()} to match jamovi's output. For
#' categorical predictors with more than one degree of freedom, \code{car::vif()} returns the
#' Generalized VIF (GVIF) and its degree-freedom-adjusted form \eqn{GVIF^{1/(2\cdot df)}}. The
#' value stored in the \code{VIF} column is the raw GVIF (or standard VIF for continuous/binary
#' predictors). Because \code{car::vif()} does not support \code{logistf} objects, VIF for
#' Firth-corrected models is computed from a parallel standard \code{glm()} fit on the same
#' formula and data.
#'
#' \strong{Significance codes:} The \code{Significance} column uses the conventional scheme:
#' \code{***} if \eqn{p < .001}, \code{**} if \eqn{p < .01}, \code{*} if \eqn{p < .05},
#' and blank otherwise.
#'
#' \strong{Inline model metrics:} Model-fit statistics (Deviance, AIC, BIC, McFadden R2,
#' CoxSnell R2, Nagelkerke R2, Tjur R2, Firth Corrected) are appended as additional columns to
#' the right of the \code{VIF} column. Values are placed only in the first row of each model
#' block (aligned with the intercept); all subsequent rows contain \code{NA} for these columns.
#'
#' @return
#' If \code{crude = TRUE}, returns a single tidy data frame with rows for all \code{indep}
#' variables stacked sequentially (separated by blank spacer rows). Columns: Predictor, Log Odds,
#' Std Error, p-value, Significance, OR, 95\% CI Low, 95\% CI High, VIF (all \code{NA}),
#' Deviance, AIC, BIC, McFadden R2, CoxSnell R2, Nagelkerke R2, Tjur R2, Firth Corrected.
#' Inline model metrics are placed in the first (intercept) row of each variable's block.
#'
#' If \code{crude = FALSE} and \code{iterate = FALSE}, returns a data frame containing Predictor,
#' Log Odds, Std Error, p-value, Significance, OR, 95\% CIs, VIF, and inline model metrics
#' (Deviance, AIC, BIC, McFadden R2, CoxSnell R2, Nagelkerke R2, Tjur R2, Firth Corrected) in
#' the first row. Model metrics are also attached as \code{attr(result, "model_metrics")}.
#'
#' If \code{crude = FALSE} and \code{iterate = TRUE}, returns a named list containing:
#' \describe{
#'   \item{Tables}{A named list of result data frames, each with inline model metrics in row 1.}
#'   \item{Metrics}{A data frame comparing model metrics across all combinations.}
#'   \item{filtered_Tables}{A subset of \code{Tables} excluding only models where \emph{all}
#'     non-intercept \code{95\% CI Low} values or \emph{all} non-intercept \code{95\% CI High}
#'     values are \code{Inf}.}
#'   \item{p_filtered_Tables}{A subset of \code{filtered_Tables} retaining only models where at
#'     least one non-intercept \code{p-value} is less than .05.}
#'   \item{vif_filtered_Tables}{A subset of \code{p_filtered_Tables} retaining only models where
#'     all non-\code{NA} VIF values are less than 5, ensuring no violation of multicollinearity.}
#' }
#'
#' @examples
#' \dontrun{
#' # --- Example 1: Single model with Firth's correction applied dynamically ---
#' result <- run_logreg(
#'   x           = my_data,
#'   y           = "Outcome",
#'   confounders = c("Age", "Sex"),
#'   indep       = c("Biomarker_A", "Biomarker_B"),
#'   ref         = "Control",
#'   ref_levels  = list("Sex" ~ "Female"),
#'   apply_firth = TRUE
#' )
#' result
#' attr(result, "model_metrics")
#'
#' # --- Example 2: Iterative models with forced Firth and VIF filtering -------
#' results <- run_logreg(
#'   x                 = my_data,
#'   y                 = "Outcome",
#'   confounders       = c("Age", "Sex", "Site"),
#'   indep             = c("Biomarker_A", "Biomarker_B", "Biomarker_C"),
#'   ref               = "Control",
#'   ref_levels        = list("Sex" ~ "Female", "Site" ~ "Main"),
#'   exclude           = "Indeterminate",
#'   iterate           = TRUE,
#'   apply_firth       = TRUE,
#'   force_apply_firth = TRUE
#' )
#' results$Metrics
#' results$vif_filtered_Tables
#'
#' # --- Example 3: Crude (unadjusted) odds ratios for each biomarker ----------
#' crude_results <- run_logreg(
#'   x           = my_data,
#'   y           = "Outcome",
#'   indep       = c("Biomarker_A", "Biomarker_B", "Biomarker_C"),
#'   ref         = "Control",
#'   ref_levels  = list("Sex" ~ "Female"),
#'   crude       = TRUE,
#'   apply_firth = TRUE
#' )
#' crude_results
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Fox J, Weisberg S (2019). \emph{An R Companion to Applied Regression}, 3rd ed. Sage.
#'
#' Firth D (1993). "Bias reduction of maximum likelihood estimates." \emph{Biometrika}, 80(1),
#' 27–38. \doi{10.1093/biomet/80.1.27}
#'
#' Heinze G, Schemper M (2002). "A solution to the problem of separation in logistic regression."
#' \emph{Statistics in Medicine}, 21(16), 2409–2419. \doi{10.1002/sim.1047}
#'
#' Heinze G, Ploner M, Jiricka L, Steuer A (2024). \emph{logistf: Firth's Bias-Reduced Logistic
#' Regression}. R package. \url{https://CRAN.R-project.org/package=logistf}
#'
#' Lüdecke D, Ben-Shachar MS, Patil I, Waggoner P, Makowski D (2021). "performance: An R Package
#' for Assessment, Comparison and Testing of Statistical Models." \emph{Journal of Open Source Software},
#' 6(60), 3139. \doi{10.21105/joss.03139}
#'
#' @importFrom stats glm binomial as.formula confint.default nobs AIC BIC terms model.matrix relevel
#' @importFrom performance r2_mcfadden r2_coxsnell r2_nagelkerke r2_tjur
#' @importFrom car vif
#' @importFrom logistf logistf
#' @importFrom utils combn
#' @export
run_logreg <- function(x, y, confounders = NULL, indep = NULL, ref, ref_levels = NULL,
                       remove = NULL, exclude = NULL, iterate = FALSE, crude = FALSE,
                       add_vif = TRUE, apply_firth = TRUE, force_apply_firth = FALSE) {

  # --- Integration with run_DIpreprocess ---
  if (inherits(x, "run_DIpreprocess")) {
    target_meta <- if (!is.null(x$metadata_merged)) x$metadata_merged else x$metadata
    target_data <- if (!is.null(x$data_nonpls_merged)) x$data_nonpls_merged else x$data_nonpls
    x <- cbind(target_meta, target_data)
  }

  # ---------------------------------------------------------------------------
  # 1. Error Handling & Input Validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x)) stop("Error: 'x' must be a data frame, matrix, or tibble.")
  x <- as.data.frame(x)

  if (!(y %in% names(x))) stop(sprintf("Error: Dependent variable '%s' not found in data.", y))

  if (crude) {
    if (is.null(indep)) {
      stop("Error: 'indep' must be specified when 'crude = TRUE'.")
    }
    if (!is.null(confounders)) {
      warning(
        "Constructive Warning: 'confounders' is ignored when 'crude = TRUE'. ",
        "Crude models are unadjusted (simple) regressions for each variable in 'indep'.",
        call. = FALSE
      )
    }
    if (iterate) {
      warning(
        "Constructive Warning: 'iterate' is ignored when 'crude = TRUE'.",
        call. = FALSE
      )
    }
  } else {
    if (is.null(confounders) && is.null(indep)) {
      stop("Error: Both 'confounders' and 'indep' cannot be NULL at the same time.")
    }
  }

  all_vars <- c(confounders, indep)
  missing_vars <- setdiff(all_vars, names(x))
  if (length(missing_vars) > 0) {
    stop(sprintf("Error: Predictors not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # 2. Data Preparation
  # ---------------------------------------------------------------------------

  # Drop specified categories from the dependent variable
  if (!is.null(exclude)) {
    x <- x[!x[[y]] %in% exclude, , drop = FALSE]
  }

  # Drop specified categories from predictor variables via 'remove'
  if (!is.null(remove)) {
    if (!is.list(remove)) {
      stop("Error: 'remove' must be a list of formulas (e.g., list(\"Sex\" ~ \"Unknown\")).")
    }

    for (mapping in remove) {
      if (!inherits(mapping, "formula") || length(mapping) != 3) {
        stop("Error: Each element in 'remove' must be a two-sided formula.")
      }

      lhs      <- mapping[[2]]
      col_name <- if (is.character(lhs)) lhs else as.character(lhs)

      rhs     <- mapping[[3]]
      cat_val <- if (is.character(rhs)) rhs else as.character(rhs)

      if (!(col_name %in% names(x))) {
        stop(sprintf("Error: Column '%s' from 'remove' not found in data.", col_name))
      }

      col_vals <- as.character(x[[col_name]])
      if (!(cat_val %in% col_vals)) {
        stop(sprintf("Error: Category '%s' not found in column '%s'.", cat_val, col_name))
      }

      x <- x[col_vals != cat_val, , drop = FALSE]
    }
  }

  # Factorize and set reference level for dependent variable
  x[[y]] <- as.factor(x[[y]])
  valid_levels <- levels(x[[y]])

  if (!(ref %in% valid_levels)) {
    stop(sprintf("Error: Reference category '%s' not found in '%s'.", ref, y))
  }
  if (length(valid_levels) != 2) {
    stop(sprintf(
      "Error: Dependent variable '%s' must have exactly 2 categories after exclusions.", y
    ))
  }
  x[[y]] <- stats::relevel(x[[y]], ref = ref)

  # Apply reference levels to predictor variables
  if (!is.null(ref_levels)) {
    if (!is.list(ref_levels)) {
      stop("Error: 'ref_levels' must be a list of formulas (e.g., list(\"Sex\" ~ \"Female\")).")
    }

    for (mapping in ref_levels) {
      if (!inherits(mapping, "formula") || length(mapping) != 3) {
        stop("Error: Each element in 'ref_levels' must be a two-sided formula.")
      }

      lhs      <- mapping[[2]]
      col_name <- if (is.character(lhs)) lhs else as.character(lhs)

      rhs     <- mapping[[3]]
      ref_val <- if (is.character(rhs)) rhs else as.character(rhs)

      if (!(col_name %in% names(x))) {
        stop(sprintf("Error: Column '%s' from ref_levels not found in data.", col_name))
      }

      x[[col_name]] <- as.factor(x[[col_name]])
      if (!(ref_val %in% levels(x[[col_name]]))) {
        stop(sprintf("Error: Category '%s' not found in column '%s'.", ref_val, col_name))
      }
      x[[col_name]] <- stats::relevel(x[[col_name]], ref = ref_val)
    }
  }

  # ---------------------------------------------------------------------------
  # Helper: wrap column names in backticks for safe formula building
  # ---------------------------------------------------------------------------
  wrap_backticks <- function(vars) {
    if (is.null(vars) || length(vars) == 0) return(NULL)
    paste0("`", vars, "`")
  }

  # ---------------------------------------------------------------------------
  # Helper: significance stars from p-value
  # ---------------------------------------------------------------------------
  sig_stars <- function(p) {
    ifelse(p < .001, "***",
      ifelse(p < .01, "**",
        ifelse(p < .05, "*", "")))
  }

  # ---------------------------------------------------------------------------
  # Helper: detect zero cells between a set of predictors and the outcome.
  # ---------------------------------------------------------------------------
  has_zero_cells <- function(pred_vars, data, y_var) {
    for (v in pred_vars) {
      col <- data[[v]]
      if (is.factor(col) || is.character(col)) {
        tbl <- table(col, data[[y_var]], useNA = "no")
        if (any(tbl == 0L)) return(TRUE)
      }
    }
    FALSE
  }

  # ---------------------------------------------------------------------------
  # Helper: extract model-fit metrics into a one-row data frame.
  # ---------------------------------------------------------------------------
  get_metrics <- function(model, model_name, firth_used) {
    suppressWarnings({
      r2_mcf <- tryCatch(performance::r2_mcfadden(model),  error = function(e) list(NA_real_))
      r2_cs  <- tryCatch(performance::r2_coxsnell(model),  error = function(e) list(NA_real_))
      r2_nag <- tryCatch(performance::r2_nagelkerke(model), error = function(e) list(NA_real_))
      r2_tjr <- tryCatch(performance::r2_tjur(model),       error = function(e) list(NA_real_))
    })

    extract_num <- function(val) {
      if (is.null(val) || length(val) == 0) return(NA_real_)
      as.numeric(val[[1]])
    }

    n_obs        <- tryCatch(as.numeric(stats::nobs(model)), error = function(e) NA_real_)
    deviance_val <- tryCatch(as.numeric(model$deviance),     error = function(e) NA_real_)
    if (is.null(deviance_val) || length(deviance_val) == 0) deviance_val <- NA_real_
    aic_val      <- tryCatch(as.numeric(stats::AIC(model)),  error = function(e) NA_real_)
    bic_val      <- tryCatch(as.numeric(stats::BIC(model)),  error = function(e) NA_real_)

    data.frame(
      Model              = model_name,
      N                  = n_obs,
      Deviance           = deviance_val,
      AIC                = aic_val,
      BIC                = bic_val,
      `McFadden R2`      = extract_num(r2_mcf),
      `CoxSnell R2`      = extract_num(r2_cs),
      `Nagelkerke R2`    = extract_num(r2_nag),
      `Tjur R2`          = extract_num(r2_tjr),
      `Firth Corrected`  = firth_used,
      stringsAsFactors   = FALSE,
      check.names        = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Helper: append inline metric columns to a result table.
  # Metric values are placed only in row 1 (intercept row); all other rows
  # receive NA. The Model and N columns from the metrics row are intentionally
  # excluded — only the fit statistics are appended.
  # ---------------------------------------------------------------------------
  append_metrics <- function(tbl, metrics_row) {
    metric_cols <- c("Deviance", "AIC", "BIC",
                     "McFadden R2", "CoxSnell R2", "Nagelkerke R2",
                     "Tjur R2", "Firth Corrected")

    n_rows <- nrow(tbl)

    for (col in metric_cols) {
      val <- metrics_row[[col]]

      if (is.logical(val)) {
        tbl[[col]] <- NA
        tbl[[col]] <- as.logical(tbl[[col]])
      } else {
        tbl[[col]] <- NA_real_
      }

      tbl[[col]][1] <- val
    }

    tbl
  }

  # ---------------------------------------------------------------------------
  # Helper: build the result table from a fitted model (glm or logistf).
  # ---------------------------------------------------------------------------
  build_table <- function(mod, combo_name = NULL, pred_vars, data, firth_used) {

    is_firth <- inherits(mod, "logistf")

    if (is_firth) {
      coef_est <- mod$coefficients
      ci_low   <- exp(mod$ci.lower)
      ci_high  <- exp(mod$ci.upper)
      p_vals   <- mod$prob
      std_err  <- sqrt(diag(mod$var))
      terms_used <- names(coef_est)

      res_df <- data.frame(
        Predictor    = terms_used,
        `Log Odds`   = as.numeric(coef_est),
        `Std Error`  = as.numeric(std_err),
        `p-value`    = as.numeric(p_vals),
        Significance = sig_stars(as.numeric(p_vals)),
        OR           = exp(as.numeric(coef_est)),
        `95% CI Low` = as.numeric(ci_low),
        `95% CI High`= as.numeric(ci_high),
        stringsAsFactors = FALSE,
        check.names  = FALSE
      )

    } else {
      coefs <- summary(mod)$coefficients
      ci    <- suppressMessages(stats::confint.default(mod))

      # --- Aliasing fix -------------------------------------------------------
      shared_rows   <- intersect(rownames(coefs), rownames(ci))
      aliased_terms <- setdiff(rownames(ci), rownames(coefs))

      if (length(aliased_terms) > 0) {
        model_label <- if (!is.null(combo_name)) sprintf(" [Model: %s]", combo_name) else ""
        warning(
          sprintf(
            "Constructive Warning%s: %d aliased term(s) detected and excluded from output due to ",
            model_label, length(aliased_terms)
          ),
          "perfect or quasi-complete separation. These levels appear in only one outcome group. ",
          "Excluded term(s): ", paste(aliased_terms, collapse = ", "), ".",
          call. = FALSE
        )
        coefs <- coefs[shared_rows, , drop = FALSE]
        ci    <- ci[shared_rows,    , drop = FALSE]
      }
      # ------------------------------------------------------------------------

      res_df <- data.frame(
        Predictor    = rownames(coefs),
        `Log Odds`   = coefs[, "Estimate"],
        `Std Error`  = coefs[, "Std. Error"],
        `p-value`    = coefs[, "Pr(>|z|)"],
        Significance = sig_stars(coefs[, "Pr(>|z|)"]),
        OR           = exp(coefs[, "Estimate"]),
        `95% CI Low` = exp(ci[, 1]),
        `95% CI High`= exp(ci[, 2]),
        stringsAsFactors = FALSE,
        check.names  = FALSE
      )
    }

    # --- VIF ------------------------------------------------------------------
    if (add_vif) {
      res_df$VIF <- NA_real_

      # VIF requires >= 2 predictors; single-predictor models (crude mode) will
      # always produce NA without attempting a car::vif() call.
      if (length(pred_vars) >= 2) {
        vif_mod <- if (is_firth) {
          tryCatch(
            stats::glm(mod$formula, data = data, family = stats::binomial(link = "logit")),
            error = function(e) NULL
          )
        } else {
          mod
        }

        if (!is.null(vif_mod)) {
          tryCatch({
            vif_raw <- car::vif(vif_mod)

            if (is.matrix(vif_raw)) {
              vif_vec <- setNames(vif_raw[, "GVIF"], rownames(vif_raw))
            } else {
              vif_vec <- vif_raw
            }

            mm          <- stats::model.matrix(vif_mod)
            assign_idx  <- attr(mm, "assign")
            term_labels <- attr(stats::terms(vif_mod), "term.labels")

            clean_labels    <- gsub("`", "", term_labels)
            clean_vif_names <- gsub("`", "", names(vif_vec))

            for (i in seq_len(nrow(res_df))) {
              cname <- res_df$Predictor[i]
              if (cname == "(Intercept)") next

              mm_col_idx <- match(cname, colnames(mm))
              if (is.na(mm_col_idx)) next

              t_idx <- assign_idx[mm_col_idx]
              if (t_idx == 0) next

              actual_label <- clean_labels[t_idx]
              v_idx        <- match(actual_label, clean_vif_names)

              if (!is.na(v_idx)) res_df$VIF[i] <- vif_vec[v_idx]
            }
          }, error = function(e) {
            warning(
              "Constructive Warning: Could not compute VIF. ",
              "Perfect multicollinearity or a single-predictor model is likely.",
              call. = FALSE
            )
          })
        }
      }
    }
    # --------------------------------------------------------------------------

    rownames(res_df) <- NULL
    return(res_df)
  }

  # ---------------------------------------------------------------------------
  # Helper: fit a single model (glm or logistf).
  # ---------------------------------------------------------------------------
  fit_model <- function(form_str, pred_vars, data, combo_name) {
    use_firth <- FALSE

    degenerate <- vapply(pred_vars, function(v) {
      col <- data[[v]]
      if (is.factor(col) || is.character(col) || is.numeric(col)) {
        length(unique(col[!is.na(col)])) < 2L
      } else {
        FALSE
      }
    }, logical(1))

    if (any(degenerate)) {
      bad <- pred_vars[degenerate]
      warning(
        sprintf(
          "Constructive Warning: Model '%s' skipped. The following predictor(s) have fewer than ",
          combo_name
        ),
        "2 unique levels in the working data and cannot be used in regression: ",
        paste(bad, collapse = ", "), ".",
        call. = FALSE
      )
      return(list(mod = NULL, firth_used = FALSE))
    }

    if (apply_firth) {
      use_firth <- force_apply_firth || has_zero_cells(pred_vars, data, y)
    }

    if (use_firth) {
      mod <- tryCatch(
        logistf::logistf(stats::as.formula(form_str), data = data, plconf = NULL),
        error = function(e) {
          warning(sprintf(
            "Constructive Warning: Firth model failed for '%s'. Falling back to glm(). Error: %s",
            combo_name, e$message
          ), call. = FALSE)
          NULL
        }
      )
      if (is.null(mod)) {
        mod <- tryCatch(
          stats::glm(stats::as.formula(form_str), data = data,
                     family = stats::binomial(link = "logit")),
          error = function(e) {
            warning(sprintf(
              "Constructive Warning: glm fallback also failed for '%s'. Error: %s",
              combo_name, e$message
            ), call. = FALSE)
            NULL
          }
        )
        use_firth <- FALSE
      }
    } else {
      mod <- tryCatch(
        stats::glm(stats::as.formula(form_str), data = data,
                   family = stats::binomial(link = "logit")),
        error = function(e) {
          warning(sprintf(
            "Constructive Warning: Model failed for combination '%s'. Error: %s",
            combo_name, e$message
          ), call. = FALSE)
          NULL
        }
      )
    }

    list(mod = mod, firth_used = use_firth)
  }

  # ---------------------------------------------------------------------------
  # Helper: build a blank spacer row with the same columns as a result table.
  # Used in crude mode to visually separate blocks across indep variables.
  # ---------------------------------------------------------------------------
  make_spacer_row <- function(tbl) {
    spacer <- tbl[NA_integer_, , drop = FALSE]
    spacer[1, ] <- NA
    spacer$Predictor <- ""
    rownames(spacer) <- NULL
    spacer
  }

  # ---------------------------------------------------------------------------
  # 3. Model Execution
  # ---------------------------------------------------------------------------

  # ============================================================
  # CRUDE MODE: one simple model per indep variable
  # ============================================================
  if (crude) {

    crude_blocks <- vector("list", length(indep))

    for (i in seq_along(indep)) {
      v          <- indep[i]
      pred_vars  <- v
      form_str   <- paste0("`", y, "` ~ `", v, "`")

      fit        <- fit_model(form_str, pred_vars, x, combo_name = v)
      mod        <- fit$mod
      firth_used <- fit$firth_used

      if (is.null(mod)) {
        # Skipped model: emit a minimal NA block so the variable is still
        # represented in the output
        na_row <- data.frame(
          Predictor    = v,
          `Log Odds`   = NA_real_,
          `Std Error`  = NA_real_,
          `p-value`    = NA_real_,
          Significance = NA_character_,
          OR           = NA_real_,
          `95% CI Low` = NA_real_,
          `95% CI High`= NA_real_,
          VIF          = NA_real_,
          Deviance     = NA_real_,
          AIC          = NA_real_,
          BIC          = NA_real_,
          `McFadden R2`   = NA_real_,
          `CoxSnell R2`   = NA_real_,
          `Nagelkerke R2` = NA_real_,
          `Tjur R2`       = NA_real_,
          `Firth Corrected` = NA,
          stringsAsFactors = FALSE,
          check.names  = FALSE
        )
        crude_blocks[[i]] <- na_row
        next
      }

      tbl         <- build_table(mod, combo_name = v, pred_vars = pred_vars,
                                 data = x, firth_used = firth_used)
      metrics_row <- get_metrics(mod, model_name = v, firth_used = firth_used)
      tbl         <- append_metrics(tbl, metrics_row)

      crude_blocks[[i]] <- tbl
    }

    # Stack blocks, inserting a spacer row between each variable's block
    stacked <- vector("list", length(crude_blocks) * 2 - 1)
    for (i in seq_along(crude_blocks)) {
      stacked[[2 * i - 1]] <- crude_blocks[[i]]
      if (i < length(crude_blocks)) {
        stacked[[2 * i]] <- make_spacer_row(crude_blocks[[i]])
      }
    }

    out <- do.call(rbind, stacked)
    rownames(out) <- NULL
    return(out)
  }

  # ============================================================
  # STANDARD MODE (unchanged from original)
  # ============================================================
  if (!iterate) {

    pred_vars  <- c(confounders, indep)
    predictors <- wrap_backticks(pred_vars)
    form_str   <- paste0("`", y, "` ~ ", paste(predictors, collapse = " + "))

    fit        <- fit_model(form_str, pred_vars, x, combo_name = "Base Model")
    mod        <- fit$mod
    firth_used <- fit$firth_used

    tbl_out     <- build_table(mod, combo_name = NULL, pred_vars = pred_vars,
                               data = x, firth_used = firth_used)
    metrics_out <- get_metrics(mod, model_name = "Base Model", firth_used = firth_used)

    tbl_out <- append_metrics(tbl_out, metrics_out)

    attr(tbl_out, "model_metrics") <- metrics_out
    return(tbl_out)

  } else {

    if (is.null(indep)) {
      stop("Error: 'iterate = TRUE' requires 'indep' variables to be specified.")
    }

    combo_list <- unlist(
      lapply(seq_along(indep), function(i) utils::combn(indep, i, simplify = FALSE)),
      recursive = FALSE
    )

    results_tables <- list()
    metrics_list   <- list()

    for (combo in combo_list) {
      combo_name <- paste(combo, collapse = " + ")
      pred_vars  <- c(confounders, combo)
      predictors <- wrap_backticks(pred_vars)
      form_str   <- paste0("`", y, "` ~ ", paste(predictors, collapse = " + "))

      fit        <- fit_model(form_str, pred_vars, x, combo_name)
      mod        <- fit$mod
      firth_used <- fit$firth_used

      if (!is.null(mod)) {
        metrics_row <- get_metrics(mod, model_name = combo_name, firth_used = firth_used)
        metrics_list[[combo_name]] <- metrics_row

        tbl <- build_table(
          mod, combo_name = combo_name, pred_vars = pred_vars,
          data = x, firth_used = firth_used
        )

        tbl <- append_metrics(tbl, metrics_row)

        results_tables[[combo_name]] <- tbl
      }
    }

    metrics_df <- do.call(rbind, metrics_list)
    rownames(metrics_df) <- NULL

    metrics_df$`Has Significant Predictor` <- vapply(
      names(results_tables),
      function(nm) {
        tbl     <- results_tables[[nm]]
        non_int <- tbl[tbl$Predictor != "(Intercept)", ]
        any(non_int$`p-value` < 0.05, na.rm = TRUE)
      },
      logical(1)
    )

    filtered_tables <- Filter(function(tbl) {
      non_int <- tbl[tbl$Predictor != "(Intercept)", ]
      low     <- non_int$`95% CI Low`
      high    <- non_int$`95% CI High`
      !(all(is.infinite(low)  | low  == 0) |
        all(is.infinite(high) | high == 0))
    }, results_tables)

    p_filtered_tables <- Filter(function(tbl) {
      non_int <- tbl[tbl$Predictor != "(Intercept)", ]
      any(non_int$`p-value` < 0.05, na.rm = TRUE)
    }, filtered_tables)

    vif_filtered_tables <- Filter(function(tbl) {
      non_int  <- tbl[tbl$Predictor != "(Intercept)", ]
      vif_vals <- non_int$VIF[!is.na(non_int$VIF)]
      length(vif_vals) > 0 && all(vif_vals < 5)
    }, p_filtered_tables)

    return(list(
      Tables              = results_tables,
      Metrics             = metrics_df,
      filtered_Tables     = filtered_tables,
      p_filtered_Tables   = p_filtered_tables,
      vif_filtered_Tables = vif_filtered_tables
    ))
  }
}