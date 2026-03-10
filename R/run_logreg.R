#' @title Perform Binary Logistic Regression with Iterative Feature Combinations
#'
#' @description
#' A robust wrapper for \code{\link[stats]{glm}} that performs binary logistic regression.
#' It automatically handles column names with spaces, computes advanced model fit metrics,
#' maps specific reference categories, removes specified categories from predictor variables,
#' and optionally iterates through all possible combinations of independent variables.
#' Output is formatted as a clean data frame.
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
#' @param iterate Logical; if \code{FALSE} (default), combines \code{confounders} and \code{indep} into a single model.
#'   If \code{TRUE}, fits separate models for all possible combinations of \code{indep} added to \code{confounders}.
#' @param add_vif Logical; if \code{TRUE} (default), adds Variance Inflation Factor (VIF) values to the output.
#'   VIFs are computed via \code{car::vif()}, which returns GVIF for categorical predictors and standard
#'   VIF for continuous predictors, consistent with jamovi's output.
#'
#' @details
#' When \code{iterate = TRUE}, the function generates all possible combinations (the power set)
#' of the \code{indep} vector. For example, if \code{indep = c("A", "B", "C")}, the combinations
#' evaluated alongside \code{confounders} will be: \code{A}, \code{B}, \code{C}, \code{A+B},
#' \code{A+C}, \code{B+C}, and \code{A+B+C}.
#'
#' \strong{Category removal (\code{remove}):} Rows where the specified predictor variable equals
#' the given category are dropped from \code{x} before any model is fitted. This is applied after
#' \code{exclude} but before factor releveling. Multiple formulas targeting the same column are
#' each applied in sequence, allowing several categories to be removed from one variable.
#'
#' \strong{VIF computation:} This function uses \code{car::vif()} to match jamovi's output. For
#' categorical predictors with more than one degree of freedom, \code{car::vif()} returns the
#' Generalized VIF (GVIF) and its degree-freedom-adjusted form \eqn{GVIF^{1/(2\cdot df)}}. The
#' value stored in the \code{VIF} column is the raw GVIF (or standard VIF for continuous/binary
#' predictors).
#'
#' \strong{Significance codes:} The \code{Significance} column uses the conventional scheme:
#' \code{***} if \eqn{p < .001}, \code{**} if \eqn{p < .01}, \code{*} if \eqn{p < .05},
#' and blank otherwise.
#'
#' @return
#' If \code{iterate = FALSE}, returns a data frame containing Predictor, Log Odds, Std Error,
#' p-value, Significance, OR, 95\% CIs, and (optionally) VIF. Model metrics are attached as
#' \code{attr(result, "model_metrics")}.
#'
#' If \code{iterate = TRUE}, returns a named list containing:
#' \describe{
#'   \item{Tables}{A named list of result data frames for each model combination.}
#'   \item{Metrics}{A data frame comparing model metrics (N, Deviance, AIC, BIC, McFadden \eqn{R^2},
#'     Cox-Snell \eqn{R^2}, Nagelkerke \eqn{R^2}, Tjur \eqn{R^2}, and Has Significant Predictor)
#'     across all combinations.}
#'   \item{filtered_Tables}{A subset of \code{Tables} excluding only models where \emph{all}
#'     non-intercept \code{95\% CI Low} values or \emph{all} non-intercept \code{95\% CI High}
#'     values are \code{Inf}.}
#'   \item{p_filtered_Tables}{A subset of \code{filtered_Tables} retaining only models where at
#'     least one non-intercept \code{p-value} is less than .05.}
#'   \item{vif_filtered_Tables}{A subset of \code{p_filtered_Tables} retaining only models where
#'     all non-\code{NA} VIF values are less than 5, ensuring no violation of multicollinearity.}
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Fox J, Weisberg S (2019). \emph{An R Companion to Applied Regression}, 3rd ed. Sage.
#'
#' Lüdecke D, Ben-Shachar MS, Patil I, Waggoner P, Makowski D (2021). "performance: An R Package
#' for Assessment, Comparison and Testing of Statistical Models." \emph{Journal of Open Source Software},
#' 6(60), 3139. \doi{10.21105/joss.03139}
#'
#' @importFrom stats glm binomial as.formula confint.default nobs AIC BIC terms model.matrix relevel
#' @importFrom performance r2_mcfadden r2_coxsnell r2_nagelkerke r2_tjur
#' @importFrom car vif
#' @importFrom utils combn
#' @export
run_logreg <- function(x, y, confounders = NULL, indep = NULL, ref, ref_levels = NULL,
                       remove = NULL, exclude = NULL, iterate = FALSE, add_vif = TRUE) {

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

  # Apply reference levels to predictor variables.
  # Accepts: list("col_name" ~ "ref_value", ...)
  # LHS is a quoted column name; RHS is the desired reference category.
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
  # Helper: extract model-fit metrics into a one-row data frame
  # ---------------------------------------------------------------------------
  get_metrics <- function(model, model_name) {
    suppressWarnings({
      r2_mcf <- performance::r2_mcfadden(model)
      r2_cs  <- performance::r2_coxsnell(model)
      r2_nag <- performance::r2_nagelkerke(model)
      r2_tjr <- performance::r2_tjur(model)
    })

    extract_num <- function(val) {
      if (is.null(val) || length(val) == 0) return(NA_real_)
      as.numeric(val[[1]])
    }

    data.frame(
      Model           = model_name,
      N               = as.numeric(stats::nobs(model)),
      Deviance        = as.numeric(model$deviance),
      AIC             = as.numeric(stats::AIC(model)),
      BIC             = as.numeric(stats::BIC(model)),
      `McFadden R2`   = extract_num(r2_mcf),
      `CoxSnell R2`   = extract_num(r2_cs),
      `Nagelkerke R2` = extract_num(r2_nag),
      `Tjur R2`       = extract_num(r2_tjr),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Helper: build the result table from a fitted glm object
  #
  # VIF note: car::vif() is used because it:
  #   (a) does NOT emit spurious interaction-term warnings,
  #   (b) returns GVIF for multi-df categorical terms (matching jamovi), and
  #   (c) gives numerically identical estimates to jamovi's logistic regression VIF output.
  #
  #   For df == 1 (continuous or binary predictor), car::vif() returns a plain scalar VIF.
  #   For df > 1 (polytomous factor), it returns a 3-column matrix:
  #     [GVIF, Df, GVIF^(1/(2*Df))]. We store the raw GVIF in both cases.
  # ---------------------------------------------------------------------------
  build_table <- function(mod) {
    coefs <- summary(mod)$coefficients
    ci    <- suppressMessages(stats::confint.default(mod))

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
      check.names = FALSE
    )

    if (add_vif) {
      res_df$VIF <- NA_real_

      tryCatch({
        vif_raw <- car::vif(mod)

        # Normalise to a named numeric vector: term_name -> GVIF value.
        # car::vif() returns a named vector for df==1 terms, or a named matrix for df>1 terms.
        if (is.matrix(vif_raw)) {
          vif_vec <- setNames(vif_raw[, "GVIF"], rownames(vif_raw))
        } else {
          vif_vec <- vif_raw
        }

        # Map each coefficient row back to its parent term label, then look up VIF.
        mm          <- stats::model.matrix(mod)
        assign_idx  <- attr(mm, "assign")           # position -> term index (0 = intercept)
        term_labels <- attr(stats::terms(mod), "term.labels")

        # Strip backticks for safe string matching
        clean_labels    <- gsub("`", "", term_labels)
        clean_vif_names <- gsub("`", "", names(vif_vec))

        for (i in seq_len(nrow(res_df))) {
          cname <- res_df$Predictor[i]
          if (cname == "(Intercept)") next

          mm_col_idx <- match(cname, colnames(mm))
          if (is.na(mm_col_idx)) next

          t_idx <- assign_idx[mm_col_idx]
          if (t_idx == 0) next                      # safety: intercept column

          actual_label <- clean_labels[t_idx]
          v_idx        <- match(actual_label, clean_vif_names)

          if (!is.na(v_idx)) {
            res_df$VIF[i] <- vif_vec[v_idx]
          }
        }
      }, error = function(e) {
        warning(
          "Constructive Warning: Could not compute VIF. ",
          "Perfect multicollinearity or a single-predictor model is likely.",
          call. = FALSE
        )
      })
    }

    rownames(res_df) <- NULL
    return(res_df)
  }

  # ---------------------------------------------------------------------------
  # 3. Model Execution
  # ---------------------------------------------------------------------------
  if (!iterate) {

    predictors <- wrap_backticks(c(confounders, indep))
    form_str   <- paste0("`", y, "` ~ ", paste(predictors, collapse = " + "))

    mod <- stats::glm(
      stats::as.formula(form_str), data = x, family = stats::binomial(link = "logit")
    )

    tbl_out     <- build_table(mod)
    metrics_out <- get_metrics(mod, model_name = "Base Model")

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
      predictors <- wrap_backticks(c(confounders, combo))
      form_str   <- paste0("`", y, "` ~ ", paste(predictors, collapse = " + "))

      mod <- tryCatch(
        stats::glm(
          stats::as.formula(form_str), data = x, family = stats::binomial(link = "logit")
        ),
        error = function(e) {
          warning(sprintf(
            "Constructive Warning: Model failed for combination '%s'. Error: %s",
            combo_name, e$message
          ), call. = FALSE)
          NULL
        }
      )

      if (!is.null(mod)) {
        results_tables[[combo_name]] <- build_table(mod)
        metrics_list[[combo_name]]   <- get_metrics(mod, model_name = combo_name)
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
    # A single Inf in one row does not disqualify the model.
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