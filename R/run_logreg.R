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
#' @param x A data frame, matrix, or tibble containing the dataset.
#' @param y A string specifying the valid, categorical column name in \code{x} to be used as the dependent variable.
#' @param confounders A character vector of valid column names in \code{x} to be used as confounding factors. Default is \code{NULL}.
#' @param indep A character vector of valid column names in \code{x} to be used as independent variables. Default is \code{NULL}.
#' @param ref A string specifying the valid category of \code{y} to be used as the reference category.
#' @param ref_levels A list of two-sided formulas mapping categorical predictor column names to their
#'   desired reference categories. The left-hand side must be a valid column name of \code{x} (quoted
#'   as a string), and the right-hand side must be a valid category of that column
#'   (e.g., \code{list("Sex" ~ "Female", "District" ~ "North")}). Default is \code{NULL}.
#' @param remove A list of two-sided formulas specifying categories to drop from predictor variables
#'   before analysis. The left-hand side must be a valid column name of \code{x} (quoted as a string),
#'   and the right-hand side must be a valid category of that column to be removed. Multiple mappings
#'   may target the same column to remove multiple categories
#'   (e.g., \code{list("Sex" ~ "Unknown", "District" ~ "Other")}). Default is \code{NULL}.
#' @param exclude A character vector of categories in \code{y} to drop before analysis. Default is \code{NULL}.
#' @param iterate Logical; if \code{FALSE} (default), combines \code{confounders} and \code{indep}
#'   into a single model. If \code{TRUE}, fits separate models for all possible combinations of
#'   \code{indep} added to \code{confounders}.
#' @param add_vif Logical; if \code{TRUE} (default), adds Variance Inflation Factor (VIF) values
#'   to the output. VIFs are computed via \code{car::vif()}, which returns GVIF for categorical
#'   predictors and standard VIF for continuous predictors, consistent with jamovi's output.
#'   Note that VIF is computed from a standard \code{glm()} fit even when Firth's correction is
#'   applied, as \code{car::vif()} does not support \code{logistf} objects.
#' @param apply_firth Logical; if \code{TRUE} (default), Firth's bias-reduced penalized-likelihood
#'   logistic regression (\code{logistf::logistf()}) is applied \emph{selectively} — only to models
#'   where at least one predictor variable has a zero cell count in its cross-tabulation with the
#'   outcome. Models with no zero cells are fitted using standard \code{glm()} regardless of this
#'   setting. A logical column \code{Firth Corrected} is added to the \code{Metrics} table
#'   (and attached metrics for \code{iterate = FALSE}) indicating whether Firth's correction was
#'   applied to each model.
#' @param force_apply_firth Logical; if \code{TRUE} (default \code{FALSE}), overrides the dynamic
#'   zero-cell detection of \code{apply_firth} and applies Firth's correction to \emph{all} models
#'   regardless of whether zero cells are present. Ignored if \code{apply_firth = FALSE}.
#'
#' @details
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
#' @return
#' If \code{iterate = FALSE}, returns a data frame containing Predictor, Log Odds, Std Error,
#' p-value, Significance, OR, 95\% CIs, and (optionally) VIF. Model metrics (including
#' \code{Firth Corrected}) are attached as \code{attr(result, "model_metrics")}.
#'
#' If \code{iterate = TRUE}, returns a named list containing:
#' \describe{
#'   \item{Tables}{A named list of result data frames for each model combination.}
#'   \item{Metrics}{A data frame comparing model metrics (N, Deviance, AIC, BIC, McFadden \eqn{R^2},
#'     Cox-Snell \eqn{R^2}, Nagelkerke \eqn{R^2}, Tjur \eqn{R^2}, Has Significant Predictor, and
#'     Firth Corrected) across all combinations.}
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
#' # Firth's correction is applied automatically only if zero cells are detected.
#' result <- run_logreg(
#'   x          = my_data,
#'   y          = "Outcome",
#'   confounders = c("Age", "Sex"),
#'   indep      = c("Biomarker_A", "Biomarker_B"),
#'   ref        = "Control",
#'   ref_levels = list("Sex" ~ "Female"),
#'   apply_firth = TRUE
#' )
#' result
#' attr(result, "model_metrics")  # includes Firth Corrected column
#'
#' # --- Example 2: Iterative models with forced Firth and VIF filtering -------
#' # All models receive Firth's correction regardless of zero-cell status.
#' results <- run_logreg(
#'   x                = my_data,
#'   y                = "Outcome",
#'   confounders      = c("Age", "Sex", "Site"),
#'   indep            = c("Biomarker_A", "Biomarker_B", "Biomarker_C"),
#'   ref              = "Control",
#'   ref_levels       = list("Sex" ~ "Female", "Site" ~ "Main"),
#'   exclude          = "Indeterminate",
#'   iterate          = TRUE,
#'   apply_firth      = TRUE,
#'   force_apply_firth = TRUE
#' )
#' results$Metrics
#' results$vif_filtered_Tables  # models with p < .05 and all VIF < 5
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
                       remove = NULL, exclude = NULL, iterate = FALSE, add_vif = TRUE,
                       apply_firth = TRUE, force_apply_firth = FALSE) {

  # ---------------------------------------------------------------------------
  # 1. Error Handling & Input Validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x)) stop("Error: 'x' must be a data frame, matrix, or tibble.")
  x <- as.data.frame(x)

  if (!(y %in% names(x))) stop(sprintf("Error: Dependent variable '%s' not found in data.", y))
  if (is.null(confounders) && is.null(indep)) {
    stop("Error: Both 'confounders' and 'indep' cannot be NULL at the same time.")
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
  # Checks all categorical predictors (factors/characters) in the current combo
  # against y. Returns TRUE if any cross-tabulation cell is zero.
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
  # Accepts both glm and logistf objects; extracts what is available from each.
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

    # nobs / deviance / AIC / BIC differ between glm and logistf
    n_obs <- tryCatch(as.numeric(stats::nobs(model)), error = function(e) NA_real_)

    deviance_val <- tryCatch(as.numeric(model$deviance), error = function(e) NA_real_)
    if (is.null(deviance_val) || length(deviance_val) == 0) deviance_val <- NA_real_

    aic_val <- tryCatch(as.numeric(stats::AIC(model)), error = function(e) NA_real_)
    bic_val <- tryCatch(as.numeric(stats::BIC(model)), error = function(e) NA_real_)

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
  # Helper: build the result table from a fitted model (glm or logistf).
  #
  # Aliasing fix (glm only): summary() drops NA coefficients; confint.default()
  # keeps them, causing a row-count mismatch. We align to shared row names and
  # warn about dropped aliased terms. logistf() does not produce aliased NAs,
  # so this check is bypassed for Firth-corrected models.
  #
  # VIF: car::vif() does not support logistf objects. For Firth-corrected models,
  # VIF is computed from a parallel glm() fit on the same formula and data.
  # ---------------------------------------------------------------------------
  build_table <- function(mod, combo_name = NULL, pred_vars, data, firth_used) {

    is_firth <- inherits(mod, "logistf")

    if (is_firth) {
      # logistf stores coefficients and CIs directly
      coef_est <- mod$coefficients
      ci_mat   <- mod$ci.lower  # will be combined below
      ci_low   <- exp(mod$ci.lower)
      ci_high  <- exp(mod$ci.upper)
      p_vals   <- mod$prob
      std_err  <- sqrt(diag(mod$var))

      # Align all vectors to the same term names (logistf includes intercept)
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
      # Standard glm path
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

      # For logistf models, fit a parallel glm for VIF computation only
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
    # --------------------------------------------------------------------------

    rownames(res_df) <- NULL
    return(res_df)
  }

  # ---------------------------------------------------------------------------
  # Helper: fit a single model (glm or logistf) given a formula string,
  # the prepared data, and the current predictor variable names.
  # Returns a list: mod (the fitted model) and firth_used (logical).
  # ---------------------------------------------------------------------------
  fit_model <- function(form_str, pred_vars, data, combo_name) {
    use_firth <- FALSE

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
      # Fallback to glm if logistf fails
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
  # 3. Model Execution
  # ---------------------------------------------------------------------------
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

    attr(tbl_out, "model_metrics") <- metrics_out
    return(tbl_out)

  } else {

    if (is.null(indep)) {
      stop("Error: 'iterate = TRUE' requires 'indep' variables to be specified.")
    }

    # Generate the power set (all non-empty subsets) of indep
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
        results_tables[[combo_name]] <- build_table(
          mod, combo_name = combo_name, pred_vars = pred_vars,
          data = x, firth_used = firth_used
        )
        metrics_list[[combo_name]] <- get_metrics(
          mod, model_name = combo_name, firth_used = firth_used
        )
      }
    }

    # Assemble metrics data frame and append significance flag
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

    # filtered_Tables: exclude only models where ALL non-intercept values in
    # "95% CI Low" are Inf, OR ALL non-intercept values in "95% CI High" are Inf.
    filtered_tables <- Filter(function(tbl) {
      non_int <- tbl[tbl$Predictor != "(Intercept)", ]
      !(all(is.infinite(non_int$`95% CI Low`)) | all(is.infinite(non_int$`95% CI High`)))
    }, results_tables)

    # p_filtered_Tables: from filtered_tables, keep only models with at least
    # one non-intercept p-value < .05
    p_filtered_tables <- Filter(function(tbl) {
      non_int <- tbl[tbl$Predictor != "(Intercept)", ]
      any(non_int$`p-value` < 0.05, na.rm = TRUE)
    }, filtered_tables)

    # vif_filtered_Tables: from p_filtered_tables, keep only models where all
    # non-NA VIF values are < 5 (no multicollinearity violation)
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