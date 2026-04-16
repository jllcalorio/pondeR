#' @title Perform Linear Mixed-Effects Model Analysis
#'
#' @description
#' Fits Linear Mixed-Effects Models (LMEMs) across all features (columns) in a
#' data matrix or data frame, using metadata for model covariates and grouping
#' structure. Results are structured to support downstream diagnostic and effect
#' visualizations including residuals vs. fitted, Q-Q plots, observed vs.
#' predicted, effect plots, and random effects caterpillar plots.
#'
#' @param x A data frame, tibble, or matrix of numeric features to be modeled
#'   (e.g., metabolites, genes). Rows must correspond to observations; columns
#'   to features. Special characters in column names are supported.
#' @param metadata A data frame or tibble with the same number of rows as
#'   \code{x}, containing covariates, grouping variables, and subject
#'   identifiers used in the model formula.
#' @param formula A character string specifying the LMEM formula. The response
#'   variable should be written as \code{.feature} as a placeholder, which will
#'   be substituted for each column in \code{x} during iteration (e.g.,
#'   \code{".feature ~ time + (1 | subject_id)"}). All variables referenced in
#'   the formula (other than \code{.feature}) must be present in
#'   \code{metadata}.
#' @param subject_col A single character string specifying the column name in
#'   \code{metadata} for subject identifiers used as random effects. Defaults
#'   to \code{NULL}, which triggers automatic detection by searching for columns
#'   named \code{"ID"}, \code{"id"}, \code{"SubjectID"}, or \code{"subject_id"}
#'   in that order. A warning is issued if no match is found and
#'   \code{subject_col} remains \code{NULL}.
#' @param group A single character string specifying a column name in
#'   \code{metadata} that encodes group membership (e.g., treatment arm,
#'   timepoint). Used for structuring summaries and downstream visualizations.
#'   Defaults to \code{NULL}.
#' @param p_adjust_method A character string specifying the p-value adjustment
#'   method for the combined fixed effects summary. Passed to
#'   \code{\link[stats]{p.adjust}}. One of \code{"holm"}, \code{"hochberg"},
#'   \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"},
#'   \code{"fdr"}, or \code{"none"}. Default is \code{"BH"}.
#' @param exclude A character vector of column names present in \code{x},
#'   \code{metadata}, or both, to be removed before any analysis. Defaults to
#'   \code{NULL} (no exclusions).
#' @param verbose Logical. If \code{TRUE}, emits a message for each feature as
#'   it is processed. Default is \code{FALSE}.
#'
#' @details
#' The function iterates over all columns in \code{x} (after exclusions),
#' binding each as the response variable into a combined analysis data frame
#' alongside \code{metadata}. The model formula is constructed by substituting
#' the literal token \code{.feature} with the backtick-quoted column name of
#' the current feature. This allows formulas with arbitrary fixed and random
#' effects structures to be reused across all features without modification.
#'
#' Models are fit using \code{\link[lmerTest:lmer]{lmerTest::lmer}} with REML
#' estimation, which provides Satterthwaite-approximated degrees of freedom and
#' p-values for fixed effects.
#'
#' For each successfully fitted model, the following diagnostic data are
#' extracted and stored to support downstream visualization:
#' \itemize{
#'   \item \strong{Residuals vs. Fitted}: raw residuals and fitted values.
#'   \item \strong{Q-Q plot}: standardized residuals and theoretical quantiles.
#'   \item \strong{Observed vs. Predicted}: response vector alongside fitted values.
#'   \item \strong{Effect plots}: fixed effect estimates with 95\% Wald CIs.
#'   \item \strong{Random effects (caterpillar)}: conditional modes and
#'     conditional standard deviations per random effect grouping factor.
#' }
#'
#' Adjusted p-values in \code{Combined_Fixed_Effects_Summary} are computed
#' separately within each unique fixed effect term across all features, using
#' the method specified in \code{p_adjust_method}.
#'
#' @return A named list of class \code{c("run_mixedmodel", "list")} with the 
#' following top-level components:
#' \describe{
#'   \item{\code{<feature_name>}}{
#'     One entry per feature in \code{x}. Each is a list with:
#'     \describe{
#'       \item{\code{status}}{Character string — \code{"success"} or failure reason.}
#'       \item{\code{model}}{Fitted \code{lmerModLmerTest} object or \code{NULL}.}
#'       \item{\code{fixed_effects}}{Data frame of fixed effects.}
#'       \item{\code{random_effects}}{From \code{\link[lme4]{VarCorr}}.}
#'       \item{\code{n_subjects}}{Integer.}
#'       \item{\code{n_obs}}{Integer.}
#'       \item{\code{convergence_warnings}}{Character vector or \code{NULL}.}
#'       \item{\code{plot_data}}{
#'         A list of plotting data:
#'         \describe{
#'           \item{\code{resid_vs_fitted}}{Columns: \code{fitted}, \code{residuals}.}
#'           \item{\code{qq}}{Columns: \code{theoretical}, \code{sample}.}
#'           \item{\code{obs_vs_pred}}{Columns: \code{observed}, \code{predicted}.}
#'           \item{\code{effects}}{Effect estimates with CIs.}
#'           \item{\code{ranef_caterpillar}}{Random effects summaries.}
#'         }
#'       }
#'     }
#'   }
#'   \item{\code{combined_fixed_effects}}{
#'     Data frame combining fixed effects from all successful models.
#'   }
#'   \item{\code{parameters}}{
#'     List of input parameters for reproducibility.
#'   }
#'   \item{\code{processing_summary}}{
#'     Summary with \code{n_features}, \code{n_success}, 
#'     \code{n_failed}, and \code{success_rate}.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n_subj  <- 20
#' n_times <- 2
#' n       <- n_subj * n_times
#'
#' metadata <- data.frame(
#'   subject_id = factor(rep(seq_len(n_subj), each = n_times)),
#'   time       = rep(c("pre", "post"), n_subj),
#'   age        = rep(round(runif(n_subj, 25, 60)), each = n_times),
#'   stringsAsFactors = FALSE
#' )
#'
#' x <- data.frame(
#'   glucose   = rnorm(n, mean = ifelse(metadata$time == "post", 5.5, 5.0), sd = 0.5),
#'   insulin   = rnorm(n, mean = ifelse(metadata$time == "post", 90,  80 ), sd = 10),
#'   `PA O-28:1` = rnorm(n, mean = 1.2, sd = 0.3),
#'   check.names = FALSE
#' )
#'
#' result <- run_mixedmodel(
#'   x        = x,
#'   metadata = metadata,
#'   formula  = ".feature ~ time + age + (1 | subject_id)",
#'   group    = "time"
#' )
#'
#' print(result)
#' summary(result)
#' head(result$combined_fixed_effects)
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom lme4 VarCorr ranef
#' @importFrom lmerTest lmer
#' @importFrom stats as.formula p.adjust fitted resid ppoints confint vcov
#'   df.residual qt
#'
#' @export
run_mixedmodel <- function(
    x,
    metadata,
    formula,
    subject_col       = NULL,
    group             = NULL,
    p_adjust_method   = "BH",
    exclude           = NULL,
    verbose           = FALSE
) {

  # ---- soft dependencies ---------------------------------------------------
  if (!requireNamespace("lme4",     quietly = TRUE))
    stop("Package 'lme4' is required. Install it with install.packages('lme4').")
  if (!requireNamespace("lmerTest", quietly = TRUE))
    stop("Package 'lmerTest' is required. Install it with install.packages('lmerTest').")

  # ---- input validation: x -------------------------------------------------
  if (!(is.data.frame(x) || is.matrix(x)))
    stop("'x' must be a data frame, tibble, or matrix.")
  if (is.matrix(x)) x <- as.data.frame(x)
  if (ncol(x) == 0L)
    stop("'x' has no columns. Provide at least one feature column.")
  if (nrow(x) == 0L)
    stop("'x' has no rows.")

  # ---- input validation: metadata ------------------------------------------
  if (!is.data.frame(metadata))
    stop("'metadata' must be a data frame or tibble.")
  if (nrow(metadata) != nrow(x))
    stop(sprintf(
      "'metadata' must have the same number of rows as 'x' (%d), but has %d.",
      nrow(x), nrow(metadata)
    ))

  # ---- input validation: formula -------------------------------------------
  if (!is.character(formula) || length(formula) != 1L || !nzchar(trimws(formula)))
    stop("'formula' must be a single non-empty character string.")
  if (!grepl("\\.feature", formula))
    stop(
      "'formula' must contain the placeholder '.feature' as the response variable.\n",
      "  Example: \".feature ~ time + (1 | subject_id)\""
    )

  # ---- input validation: p_adjust_method -----------------------------------
  valid_adj <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!p_adjust_method %in% valid_adj)
    stop(sprintf(
      "'p_adjust_method' must be one of: %s.",
      paste(valid_adj, collapse = ", ")
    ))

  # ---- input validation: verbose -------------------------------------------
  if (!is.logical(verbose) || length(verbose) != 1L)
    stop("'verbose' must be a single logical value (TRUE or FALSE).")

  # ---- exclude columns -----------------------------------------------------
  if (!is.null(exclude)) {
    if (!is.character(exclude))
      stop("'exclude' must be a character vector of column names.")
    bad_excl <- setdiff(exclude, c(names(x), names(metadata)))
    if (length(bad_excl))
      warning(sprintf(
        "The following names in 'exclude' were not found in 'x' or 'metadata' and will be ignored: %s.",
        paste(bad_excl, collapse = ", ")
      ))
    x        <- x[,        !names(x)        %in% exclude, drop = FALSE]
    metadata <- metadata[, !names(metadata) %in% exclude, drop = FALSE]
    if (ncol(x) == 0L)
      stop("All columns in 'x' were removed by 'exclude'. Provide at least one feature.")
  }

  # ---- subject_col auto-detection ------------------------------------------
  if (is.null(subject_col)) {
    candidates <- c("ID", "id", "SubjectID", "subject_id")
    found      <- candidates[candidates %in% names(metadata)]
    if (length(found)) {
      subject_col <- found[[1L]]
      message(sprintf("'subject_col' not specified; using detected column: '%s'.", subject_col))
    } else {
      warning(
        "'subject_col' is NULL and no default column (ID, id, SubjectID, subject_id) ",
        "was found in 'metadata'. Random effects depending on subject identity may ",
        "fail unless the formula specifies a valid grouping factor explicitly."
      )
    }
  } else {
    if (!is.character(subject_col) || length(subject_col) != 1L)
      stop("'subject_col' must be a single character string.")
    if (!subject_col %in% names(metadata))
      stop(sprintf(
        "'subject_col' (\"%s\") was not found in 'metadata'. Available columns: %s.",
        subject_col, paste(names(metadata), collapse = ", ")
      ))
  }

  # ---- group validation ----------------------------------------------------
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L)
      stop("'group' must be a single character string.")
    if (!group %in% names(metadata))
      stop(sprintf(
        "'group' (\"%s\") was not found in 'metadata'. Available columns: %s.",
        group, paste(names(metadata), collapse = ", ")
      ))
  }

  # ---- parse formula variables (excluding .feature) -----------------------
  formula_vars <- all.vars(as.formula(sub("\\.feature", "RESPONSE_PLACEHOLDER", formula)))
  formula_vars <- setdiff(formula_vars, "RESPONSE_PLACEHOLDER")
  missing_from_meta <- setdiff(formula_vars, names(metadata))
  if (length(missing_from_meta))
    stop(sprintf(
      "The following variables appear in 'formula' but are absent from 'metadata': %s.",
      paste(missing_from_meta, collapse = ", ")
    ))

  # ---- feature names -------------------------------------------------------
  feature_names <- names(x)

  # ---- internal helpers ----------------------------------------------------

  # Safely wrap a column name in backticks for formula construction
  .bt <- function(nm) paste0("`", nm, "`")

  # Extract caterpillar data for all random grouping factors
  .ranef_caterpillar <- function(model) {
    re      <- lme4::ranef(model, condVar = TRUE)
    result  <- vector("list", length(re))
    names(result) <- names(re)

    for (grp in names(re)) {
      re_df   <- re[[grp]]
      cvar    <- attr(re[[grp]], "postVar")  # 3-D array: [q, q, n_levels]
      n_terms <- ncol(re_df)
      n_lvls  <- nrow(re_df)

      rows <- vector("list", n_terms * n_lvls)
      idx  <- 0L
      for (j in seq_len(n_terms)) {
        term_nm <- colnames(re_df)[j]
        for (i in seq_len(n_lvls)) {
          se_val <- if (!is.null(cvar)) sqrt(cvar[j, j, i]) else NA_real_
          idx    <- idx + 1L
          rows[[idx]] <- data.frame(
            level    = rownames(re_df)[i],
            term     = term_nm,
            estimate = re_df[i, j],
            se       = se_val,
            ci_lower = re_df[i, j] - 1.96 * se_val,
            ci_upper = re_df[i, j] + 1.96 * se_val,
            stringsAsFactors = FALSE
          )
        }
      }
      result[[grp]] <- do.call(rbind, rows)
    }
    result
  }

  # Extract and standardise the fixed-effects coefficient table
  .fixed_effects_df <- function(model_summary) {
    fe  <- as.data.frame(model_summary$coefficients)
    # Standardise column names defensively
    col_map <- c(
      "Estimate"   = "estimate",
      "Std. Error" = "std_error",
      "df"         = "df",
      "t value"    = "t_value",
      "Pr(>|t|)"   = "p_value"
    )
    for (orig in names(col_map)) {
      if (orig %in% names(fe)) names(fe)[names(fe) == orig] <- col_map[[orig]]
    }
    fe$term <- rownames(fe)
    rownames(fe) <- NULL
    fe[, c("term", intersect(c("estimate","std_error","df","t_value","p_value"), names(fe)))]
  }

  # ---- main loop -----------------------------------------------------------
  all_results   <- vector("list", length(feature_names))
  names(all_results) <- feature_names

  n_success <- 0L
  n_failed  <- 0L

  for (feat in feature_names) {

    if (verbose) message(sprintf("  Processing: %s", feat))

    result_stub <- list(
      status              = NULL,
      model               = NULL,
      fixed_effects       = NULL,
      random_effects      = NULL,
      plot_data           = NULL,
      n_subjects          = NA_integer_,
      n_obs               = NA_integer_,
      convergence_warnings = NULL
    )

    # -- build analysis data frame -------------------------------------------
    feat_vec <- x[[feat]]

    if (!is.numeric(feat_vec)) {
      result_stub$status <- "skipped: feature is not numeric"
      all_results[[feat]] <- result_stub
      n_failed <- n_failed + 1L
      if (verbose) message(sprintf("    Skipped (non-numeric): %s", feat))
      next
    }

    df_model <- cbind(
      data.frame(.feature = feat_vec, stringsAsFactors = FALSE),
      metadata
    )

    # Remove incomplete cases restricted to formula variables
    vars_needed  <- c(".feature", formula_vars)
    vars_present <- intersect(vars_needed, names(df_model))
    complete_idx <- stats::complete.cases(df_model[, vars_present, drop = FALSE])
    df_model     <- df_model[complete_idx, , drop = FALSE]

    if (nrow(df_model) == 0L) {
      result_stub$status <- "skipped: no complete cases"
      all_results[[feat]] <- result_stub
      n_failed <- n_failed + 1L
      next
    }

    n_obs  <- nrow(df_model)
    n_subj <- if (!is.null(subject_col)) length(unique(df_model[[subject_col]])) else NA_integer_

    if (!is.null(subject_col) && !is.na(n_subj) && n_subj < 3L) {
      result_stub$status    <- sprintf("skipped: too few subjects (%d; need >= 3)", n_subj)
      result_stub$n_subjects <- n_subj
      result_stub$n_obs      <- n_obs
      all_results[[feat]]    <- result_stub
      n_failed <- n_failed + 1L
      next
    }

    # -- build formula with backtick-safe response ---------------------------
    feat_formula_str <- sub("\\.feature", ".feature", formula)  # .feature already in df
    feat_formula_obj <- tryCatch(
      as.formula(feat_formula_str),
      error = function(e)
        stop(sprintf("Could not parse formula for feature '%s': %s", feat, e$message))
    )

    # -- fit model -----------------------------------------------------------
    fit <- tryCatch({
      if (verbose) {
        lmerTest::lmer(feat_formula_obj, data = df_model, REML = TRUE)
      } else {
        suppressMessages(suppressWarnings(
          lmerTest::lmer(feat_formula_obj, data = df_model, REML = TRUE)
        ))
      }
    }, error = function(e) {
      structure(list(message = e$message), class = "lmer_error")
    })

    if (inherits(fit, "lmer_error")) {
      result_stub$status    <- paste("failed: model fitting error -", fit$message)
      result_stub$n_subjects <- n_subj
      result_stub$n_obs      <- n_obs
      all_results[[feat]]    <- result_stub
      n_failed <- n_failed + 1L
      next
    }

    # -- convergence warnings ------------------------------------------------
    conv_msgs <- fit@optinfo$conv$lme4$messages
    if (length(conv_msgs) == 0L) conv_msgs <- NULL

    # -- summary object ------------------------------------------------------
    fit_summary <- summary(fit)

    # -- fixed effects -------------------------------------------------------
    fe_df <- tryCatch(.fixed_effects_df(fit_summary), error = function(e) NULL)

    # -- confidence intervals (Wald) -----------------------------------------
    ci_df <- tryCatch({
      ci <- suppressMessages(confint(fit, level = 0.95, method = "Wald"))
      ci <- as.data.frame(ci)
      ci$term <- rownames(ci)
      rownames(ci) <- NULL
      names(ci)[1:2] <- c("ci_lower", "ci_upper")
      ci
    }, error = function(e) {
      if (!is.null(fe_df))
        data.frame(term = fe_df$term, ci_lower = NA_real_, ci_upper = NA_real_)
      else NULL
    })

    if (!is.null(fe_df) && !is.null(ci_df)) {
      fe_df <- merge(fe_df, ci_df, by = "term", all.x = TRUE)
    }

    # -- random effects ------------------------------------------------------
    re_df <- tryCatch(
      as.data.frame(lme4::VarCorr(fit)),
      error = function(e) NULL
    )

    # -- plot_data assembly --------------------------------------------------

    # 1. residuals vs fitted
    raw_resid   <- stats::resid(fit)
    fitted_vals <- stats::fitted(fit)

    resid_vs_fitted <- data.frame(
      fitted    = fitted_vals,
      residuals = raw_resid
    )

    # 2. Q-Q (standardised residuals)
    std_resid <- tryCatch(
      stats::resid(fit, scaled = TRUE),
      error = function(e) as.vector(scale(raw_resid))
    )
    dof_resid <- tryCatch(stats::df.residual(fit), error = function(e) length(std_resid) - 1L)
    qq_data <- data.frame(
      theoretical = stats::qt(stats::ppoints(length(std_resid)), df = dof_resid),
      sample      = sort(std_resid)
    )

    # 3. observed vs predicted
    obs_vs_pred <- data.frame(
      observed  = df_model[[".feature"]],
      predicted = fitted_vals
    )

    # 4. effect plot data (fixed effect CIs, excluding intercept)
    effects_data <- if (!is.null(fe_df)) {
      fe_excl <- fe_df[fe_df$term != "(Intercept)", , drop = FALSE]
      fe_excl[, intersect(c("term","estimate","ci_lower","ci_upper"), names(fe_excl)), drop = FALSE]
    } else NULL

    # 5. random effects caterpillar
    ranef_data <- tryCatch(.ranef_caterpillar(fit), error = function(e) NULL)

    plot_data <- list(
      resid_vs_fitted  = resid_vs_fitted,
      qq               = qq_data,
      obs_vs_pred      = obs_vs_pred,
      effects          = effects_data,
      ranef_caterpillar = ranef_data
    )

    # -- store ---------------------------------------------------------------
    all_results[[feat]] <- list(
      status               = "success",
      model                = fit,
      fixed_effects        = fe_df,
      random_effects       = re_df,
      plot_data            = plot_data,
      n_subjects           = n_subj,
      n_obs                = n_obs,
      convergence_warnings = conv_msgs
    )
    n_success <- n_success + 1L
  }

  # ---- combined fixed effects summary --------------------------------------
  successful <- Filter(function(r) isTRUE(r$status == "success"), all_results)

  combined_fe <- if (length(successful)) {
    rows <- lapply(names(successful), function(feat) {
      fe <- successful[[feat]]$fixed_effects
      if (is.null(fe) || nrow(fe) == 0L) return(NULL)
      fe$feature <- feat
      fe
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) {
      out <- do.call(rbind, rows)
      rownames(out) <- NULL
      # column order
      core_cols  <- c("feature","term","estimate","std_error","df","t_value","p_value","ci_lower","ci_upper")
      out <- out[, c(intersect(core_cols, names(out)),
                     setdiff(names(out), core_cols)), drop = FALSE]

      # p-value adjustment per term
      out$p_adj <- out$p_value
      if (p_adjust_method != "none" && "p_value" %in% names(out)) {
        for (trm in unique(out$term)) {
          idx <- out$term == trm & !is.na(out$p_value)
          if (any(idx)) out$p_adj[idx] <- stats::p.adjust(out$p_value[idx], method = p_adjust_method)
        }
      }
      # significance stars
      out$stars <- ifelse(
        is.na(out$p_adj), "",
        ifelse(out$p_adj < 0.001, "***",
               ifelse(out$p_adj < 0.01, "**",
                      ifelse(out$p_adj < 0.05, "*", "")))
      )
      out
    } else data.frame()
  } else {
    warning("No successful models were produced. 'combined_fixed_effects' is empty.")
    data.frame()
  }

  # ---- assemble return object ----------------------------------------------
  out_list <- c(
    all_results,
    list(
      combined_fixed_effects = combined_fe,
      parameters = list(
        formula         = formula,
        subject_col     = subject_col,
        group           = group,
        p_adjust_method = p_adjust_method,
        exclude         = exclude,
        n_features_in   = length(feature_names)
      ),
      processing_summary = list(
        n_features   = length(feature_names),
        n_success    = n_success,
        n_failed     = n_failed,
        success_rate = round(n_success / length(feature_names) * 100, 2)
      )
    )
  )

  class(out_list) <- c("run_mixedmodel", "list")

  message(sprintf(
    "\nrun_mixedmodel complete: %d/%d features fitted successfully (%.1f%%).",
    n_success, length(feature_names),
    n_success / length(feature_names) * 100
  ))

  invisible(out_list)
}


# ---- S3 methods -------------------------------------------------------------

#' @export
print.run_mixedmodel <- function(x, ...) {
  ps <- x$processing_summary
  cat("=== run_mixedmodel ===\n")
  cat(sprintf("  Features analysed : %d\n", ps$n_features))
  cat(sprintf("  Successful models : %d\n", ps$n_success))
  cat(sprintf("  Failed / skipped  : %d\n", ps$n_failed))
  cat(sprintf("  Success rate      : %.1f%%\n", ps$success_rate))
  invisible(x)
}

#' @export
summary.run_mixedmodel <- function(object, n = 10L, ...) {
  cfe <- object$combined_fixed_effects
  if (!is.null(cfe) && nrow(cfe) > 0L && "p_adj" %in% names(cfe)) {
    top <- head(cfe[order(cfe$p_adj, na.last = TRUE), ], n)
  } else {
    top <- data.frame()
  }
  ans <- list(
    processing_summary = object$processing_summary,
    parameters         = object$parameters,
    top_effects        = top
  )
  class(ans) <- "summary.run_mixedmodel"
  ans
}

#' @export
print.summary.run_mixedmodel <- function(x, ...) {
  ps <- x$processing_summary
  cat("---------------------------------------\n")
  cat("run_mixedmodel — Analysis Summary\n")
  cat("---------------------------------------\n")
  cat(sprintf("  Formula       : %s\n", x$parameters$formula))
  cat(sprintf("  Subject col   : %s\n", x$parameters$subject_col %||% "(not specified)"))
  cat(sprintf("  P-adj method  : %s\n", x$parameters$p_adjust_method))
  cat(sprintf("  Features in   : %d\n", ps$n_features))
  cat(sprintf("  Converged     : %d\n", ps$n_success))
  cat(sprintf("  Failed        : %d\n", ps$n_failed))
  cat("\n-- Top Features by Adjusted P-value --\n")
  if (nrow(x$top_effects) > 0L) {
    print(x$top_effects[, intersect(
      c("feature","term","estimate","p_value","p_adj","stars"), names(x$top_effects)
    )], row.names = FALSE)
  } else {
    cat("  No significant effects found.\n")
  }
  invisible(x)
}

# Lightweight null-coalescing operator (internal only; not exported)
`%||%` <- function(a, b) if (is.null(a)) b else a