# =============================================================================
#  SECTION 1 — Internal helpers (unexported)
# =============================================================================

# -----------------------------------------------------------------------------
#' Check normality of a numeric vector
#'
#' @param data_vector Numeric vector to test.
#' @param group_label Optional label for the group (used in the returned label).
#' @param method One of \code{"auto"}, \code{"shapiro"}, \code{"lilliefors"}.
#' @param alpha_normality Significance level.
#' @param warnings_list Character vector of accumulated warnings; updated
#'   in-place via explicit return (no superassignment).
#'
#' @return A named list with elements \code{violated}, \code{p_value},
#'   \code{method}, \code{label}, and \code{warnings_list}.
#' @keywords internal
.check_normality <- function(data_vector,
                             group_label    = NULL,
                             method         = "auto",
                             alpha_normality = 0.05,
                             warnings_list  = character(0)) {

  label <- if (!is.null(group_label))
    sprintf("Normality: %s", group_label) else "Normality"

  n <- length(data_vector)
  if (n < 3)
    return(list(violated = TRUE, p_value = NA_real_,
                method = "Insufficient data", label = label,
                warnings_list = warnings_list))

  if (length(unique(data_vector)) < 2) {
    msg <- if (!is.null(group_label))
      sprintf("All values in group '%s' are identical; normality assumption violated.",
              group_label)
    else
      "All values are identical; normality assumption violated."
    warnings_list <- c(warnings_list, msg)
    return(list(violated = TRUE, p_value = NA_real_,
                method = "Skipped (zero variance)", label = label,
                warnings_list = warnings_list))
  }

  use_method <- if (method == "auto") {
    if (n < 50) "shapiro" else "lilliefors"
  } else method

  if (use_method == "shapiro") {
    res   <- stats::shapiro.test(data_vector)
    p_val <- res$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method = "Shapiro-Wilk test", label = label,
                warnings_list = warnings_list))
  }

  # lilliefors
  if (!requireNamespace("nortest", quietly = TRUE)) {
    warnings_list <- c(warnings_list,
      "Package 'nortest' not available. Falling back to Shapiro-Wilk test.")
    res   <- stats::shapiro.test(data_vector)
    p_val <- res$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method = "Shapiro-Wilk test (fallback)", label = label,
                warnings_list = warnings_list))
  }
  res   <- nortest::lillie.test(data_vector)
  p_val <- res$p.value
  list(violated = p_val < alpha_normality, p_value = p_val,
       method = "Lilliefors (K-S) test", label = label,
       warnings_list = warnings_list)
}


# -----------------------------------------------------------------------------
#' Run normality checks in parallel (or sequentially)
#'
#' Encapsulates the Windows / POSIX branching so it is not duplicated between
#' \code{run_diff} and \code{.run_diff_single}.
#'
#' @param groups Character vector of group labels.
#' @param y Numeric outcome vector.
#' @param grp Factor grouping vector (same length as \code{y}).
#' @param normality_method Passed to \code{.check_normality}.
#' @param alpha_normality Passed to \code{.check_normality}.
#' @param num_cores Number of cores to use.
#'
#' @return A list of normality-check result lists, one per group.
#'   Each element carries its own \code{warnings_list} slot; callers
#'   must harvest these.
#' @keywords internal
.parallel_normality <- function(groups, y, grp,
                                normality_method, alpha_normality,
                                num_cores) {

  run_one <- function(g)
    .check_normality(y[grp == g], g, normality_method, alpha_normality)

  if (num_cores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(num_cores)
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      parallel::clusterExport(
        cl,
        varlist = c("y", "grp", ".check_normality",
                    "normality_method", "alpha_normality"),
        envir   = environment()
      )
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("nortest", quietly = TRUE)) library(nortest)
      })
      parallel::parLapply(cl, groups, run_one)
    } else {
      parallel::mclapply(groups, run_one, mc.cores = num_cores)
    }
  } else {
    lapply(groups, run_one)
  }
}


# -----------------------------------------------------------------------------
#' Run variance-homogeneity test (Levene or Bartlett fallback)
#'
#' @param y Numeric outcome vector.
#' @param grp Factor grouping vector.
#' @param alpha_variance Significance level.
#' @param verbose Logical.
#' @param warnings_list Character vector of accumulated warnings.
#'
#' @return A list with \code{p_value}, \code{test_name}, \code{violated},
#'   and \code{warnings_list}.
#' @keywords internal
.check_variance <- function(y, grp, alpha_variance,
                            verbose, warnings_list) {

  result <- tryCatch({
    if (requireNamespace("car", quietly = TRUE)) {
      list(p     = car::leveneTest(y, grp)$`Pr(>F)`[1],
           name  = "Levene's test")
    } else {
      if (verbose)
        message("Package 'car' not found. Using Bartlett's test.")
      list(p     = stats::bartlett.test(x = y, g = grp)$p.value,
           name  = "Bartlett's test")
    }
  }, error = function(e) {
    list(p = 1, name = "Error (assumed homogeneous)")
  })

  if (result$name == "Error (assumed homogeneous)") {
    warnings_list <- c(warnings_list,
                       "Error in variance test. Assuming homogeneity.")
    if (verbose) message("Error in variance test. Assuming homogeneity.")
  }

  list(p_value       = result$p,
       test_name     = result$name,
       violated      = !is.na(result$p) && result$p < alpha_variance,
       warnings_list = warnings_list)
}


# -----------------------------------------------------------------------------
#' Build the standard assumptions data frame row from a normality result
#' @keywords internal
.assumption_row <- function(label, method_name, p_value, violated) {
  data.frame(
    check    = label,
    result   = method_name,
    p_value  = p_value,
    decision = ifelse(violated, "VIOLATED", "OK"),
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
#' Format a degrees-of-freedom string consistently across test types
#'
#' For t-tests / Wilcoxon the \code{parameter} slot is a single named number;
#' for ANOVA-like tests it is a named vector of length 2 or more.
#'
#' @param parameter Named numeric vector from an \code{htest} or equivalent.
#' @return A single character string, e.g.
#'   \code{"df = 18"} or \code{"DFn = 2, DFd = 27"}.
#' @keywords internal
.format_df <- function(parameter) {
  if (is.null(parameter) || length(parameter) == 0)
    return(NA_character_)
  paste(names(parameter), round(as.numeric(parameter), 4),
        sep = " = ", collapse = ", ")
}


# =============================================================================
#  SECTION 2 — Test-execution helpers (unexported)
#  Each returns list(test_used, test_result, fitted_model, parametric,
#                    warnings_list)
# =============================================================================

# -----------------------------------------------------------------------------
#' Execute an independent two-sample test
#' @keywords internal
.run_independent_two_sample <- function(y, groups, grp,
                                        normality_violated,
                                        variance_violated,
                                        test_type,
                                        warnings_list) {

  parametric    <- TRUE
  fitted_model  <- NULL
  test_used     <- NULL
  test_result   <- NULL

  if (test_type %in% c("auto", "nonparametric")) {
    do_nonpar <- (test_type == "nonparametric") ||
                 (test_type == "auto" && normality_violated)
    if (do_nonpar) {
      label <- if (test_type == "nonparametric")
        "Mann-Whitney U test (forced)" else "Mann-Whitney U test (Wilcoxon rank-sum)"
      res <- withCallingHandlers(
        stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
        warning = function(w) {
          if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
            warnings_list <<- c(warnings_list, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      return(list(test_used = label, test_result = res,
                  fitted_model = NULL, parametric = FALSE,
                  warnings_list = warnings_list))
    }
  }

  # parametric branch
  if (test_type == "auto" && variance_violated) {
    return(list(test_used    = "Welch's t-test",
                test_result  = stats::t.test(x = y[grp == groups[1]],
                                             y = y[grp == groups[2]],
                                             var.equal = FALSE),
                fitted_model = NULL, parametric = TRUE,
                warnings_list = warnings_list))
  }

  label <- if (test_type == "parametric") "Student's t-test (forced)" else "Student's t-test"
  list(test_used    = label,
       test_result  = stats::t.test(x = y[grp == groups[1]],
                                    y = y[grp == groups[2]],
                                    var.equal = TRUE),
       fitted_model = NULL, parametric = TRUE,
       warnings_list = warnings_list)
}


# -----------------------------------------------------------------------------
#' Execute a paired two-sample test
#' @keywords internal
.run_paired_two_sample <- function(y, groups, grp,
                                   normality_violated,
                                   test_type,
                                   warnings_list) {

  g1 <- y[grp == groups[1]]
  g2 <- y[grp == groups[2]]

  if (test_type %in% c("auto", "nonparametric")) {
    do_nonpar <- (test_type == "nonparametric") ||
                 (test_type == "auto" && normality_violated)
    if (do_nonpar) {
      label <- if (test_type == "nonparametric")
        "Wilcoxon signed-rank test (forced)" else "Wilcoxon signed-rank test"
      res <- withCallingHandlers(
        stats::wilcox.test(g1, g2, paired = TRUE),
        warning = function(w) {
          if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
            warnings_list <<- c(warnings_list, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      return(list(test_used = label, test_result = res,
                  fitted_model = NULL, parametric = FALSE,
                  warnings_list = warnings_list))
    }
  }

  label <- if (test_type == "parametric") "Paired t-test (forced)" else "Paired t-test"
  list(test_used    = label,
       test_result  = stats::t.test(g1, g2, paired = TRUE),
       fitted_model = NULL, parametric = TRUE,
       warnings_list = warnings_list)
}


# -----------------------------------------------------------------------------
#' Execute an independent ANOVA-family test
#' @keywords internal
.run_independent_anova <- function(y, grp,
                                   normality_violated,
                                   variance_violated,
                                   test_type,
                                   warnings_list) {

  if (test_type %in% c("auto", "nonparametric")) {
    do_nonpar <- (test_type == "nonparametric") ||
                 (test_type == "auto" && normality_violated)
    if (do_nonpar) {
      label <- if (test_type == "nonparametric")
        "Kruskal-Wallis test (forced)" else "Kruskal-Wallis test"
      return(list(test_used = label,
                  test_result  = stats::kruskal.test(x = y, g = grp),
                  fitted_model = NULL, parametric = FALSE,
                  warnings_list = warnings_list))
    }
  }

  if (test_type == "auto" && variance_violated) {
    return(list(test_used    = "Welch's ANOVA",
                test_result  = stats::oneway.test(stats::formula(y ~ grp),
                                                  var.equal = FALSE),
                fitted_model = NULL, parametric = TRUE,
                warnings_list = warnings_list))
  }

  label        <- if (test_type == "parametric") "One-way ANOVA (forced)" else "One-way ANOVA"
  fitted_model <- stats::aov(y ~ grp)
  sa           <- summary(fitted_model)[[1]]
  test_result  <- list(
    statistic = stats::setNames(sa$`F value`[1], "F value"),
    parameter = stats::setNames(c(sa$Df[1], sa$Df[2]), c("df1", "df2")),
    p.value   = sa$`Pr(>F)`[1],
    method    = label,
    data.name = deparse(substitute(y))
  )
  list(test_used = label, test_result = test_result,
       fitted_model = fitted_model, parametric = TRUE,
       warnings_list = warnings_list)
}


# =============================================================================
#  SECTION 3 — Post-hoc construction helper (unexported)
# =============================================================================

# -----------------------------------------------------------------------------
#' Build a standardised post-hoc result data frame
#'
#' Accepts the raw output from any of \code{rstatix::tukey_hsd},
#' \code{rstatix::games_howell_test}, or \code{rstatix::dunn_test} and adds
#' \code{n1}, \code{n2}, \code{mean1}, \code{mean2}, formatted p columns,
#' \code{significant}, \code{p.adj.signif}, and \code{interpretation}.
#'
#' @param ph Raw data frame from rstatix.
#' @param posthoc_label Character string naming the post-hoc test used.
#' @param y Numeric outcome vector.
#' @param grp Factor grouping vector.
#' @param group_ns Named integer table of per-group n.
#' @param groups Character vector of group names (levels).
#' @param parametric Logical — use mean-based interpretation when \code{TRUE}.
#' @param test_alpha Significance threshold.
#' @param p_adjust_method p-adjustment method (used only for Dunn's; others
#'   already carry p.adj from rstatix).
#'
#' @return A data frame with a standardised column set.
#' @keywords internal
.build_posthoc_df <- function(ph, posthoc_label, y, grp,
                              group_ns, groups, parametric,
                              test_alpha, p_adjust_method) {

  ph <- as.data.frame(ph)
  ph$posthoc_test <- posthoc_label

  # Drop rstatix's hardcoded p.adj.signif so we can recompute with test_alpha
  ph$p.adj.signif <- NULL

  # ---- ensure n1, n2 ----
  if (!all(c("n1", "n2") %in% names(ph))) {
    n_lookup <- stats::setNames(as.numeric(group_ns), names(group_ns))
    ph$n1    <- n_lookup[match(ph$group1, names(n_lookup))]
    ph$n2    <- n_lookup[match(ph$group2, names(n_lookup))]
  }

  # ---- raw p ----
  if (!"p" %in% names(ph))  ph$p <- NA_real_
  ph$p.format <- ifelse(
    is.na(ph$p), NA_character_,
    ifelse(ph$p < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p))
  )

  # ---- p.adj ----
  if (!"p.adj" %in% names(ph)) ph$p.adj <- NA_real_
  ph$p.adj.format <- ifelse(
    is.na(ph$p.adj), NA_character_,
    ifelse(ph$p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p.adj))
  )

  ph$significant <- !is.na(ph$p.adj) & ph$p.adj < test_alpha

  ph$p.adj.signif <- ifelse(
    ph$significant,
    ifelse(ph$p.adj < 0.001, "***",
    ifelse(ph$p.adj < 0.01,  "**",
    ifelse(ph$p.adj < 0.05,  "*",
           sprintf("* (p<%.2f)", test_alpha)))),
    "ns"
  )

  # ---- means ----
  summary_stats <- data.frame(
    group = groups,
    mean  = as.numeric(tapply(y, grp, mean)),
    stringsAsFactors = FALSE
  )
  ph$mean1 <- summary_stats$mean[match(ph$group1, summary_stats$group)]
  ph$mean2 <- summary_stats$mean[match(ph$group2, summary_stats$group)]

  # ---- interpretation ----
  if (parametric) {
    ph$interpretation <- ifelse(
      !ph$significant,
      "No significant difference between groups",
      ifelse(ph$mean1 < ph$mean2,
             paste0(ph$group1, " has lower mean than ",  ph$group2),
             paste0(ph$group1, " has higher mean than ", ph$group2))
    )
  } else {
    rks        <- rank(y)
    mean_ranks <- tapply(rks, grp, mean)
    ph$interpretation <- mapply(
      function(g1, g2, sig) {
        if (!isTRUE(sig)) return("No significant difference between groups")
        if (mean_ranks[g1] < mean_ranks[g2])
          paste0(g1, " tends to have smaller values than ", g2)
        else
          paste0(g1, " tends to have larger values than ", g2)
      },
      ph$group1, ph$group2, ph$significant,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )
  }

  # ---- canonical column order ----
  final_cols <- c("group1", "group2", "n1", "n2", "statistic",
                  "p", "p.adj", "posthoc_test",
                  "p.format", "p.adj.format", "p.adj.signif",
                  "significant", "mean1", "mean2", "interpretation")
  for (col in final_cols) if (!col %in% names(ph)) ph[[col]] <- NA
  ph[, final_cols]
}


# =============================================================================
#  SECTION 4 — RM / Mixed-ANOVA helper (unexported)
# =============================================================================

# -----------------------------------------------------------------------------
#' Execute RM or mixed ANOVA and return a full run_diff-compatible result list
#'
#' Called exclusively from the single-outcome path of \code{run_diff} when
#' \code{within} is non-NULL.
#'
#' @return A list of class \code{"run_diff"}.
#' @keywords internal
.run_rm_anova <- function(x, outcome, group, within, subject_id,
                          test_family,
                          normality_method, alpha_normality,
                          alpha_sphericity,
                          type_sosquares, rm_effect_size, correction,
                          test_alpha, perform_posthoc, p_adjust_method,
                          verbose, warnings_list, parameters) {

  # ---- 0.  Build analysis data frame ----------------------------------------
  key_cols <- unique(c(subject_id, within, outcome,
                       if (!is.null(group)) group))
  df_rm    <- x[, key_cols, drop = FALSE]

  for (wf in within) df_rm[[wf]] <- as.factor(df_rm[[wf]])
  if (!is.null(group)) df_rm[[group]] <- as.factor(df_rm[[group]])
  df_rm[[subject_id]] <- as.factor(df_rm[[subject_id]])

  cc_rm <- stats::complete.cases(df_rm)
  if (sum(!cc_rm) > 0) {
    msg <- sprintf("Removing %d rows with missing values (RM/Mixed ANOVA).",
                   sum(!cc_rm))
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }
  df_rm <- df_rm[cc_rm, , drop = FALSE]

  # Safe internal column name to avoid formula-parsing issues
  safe_outcome         <- ".outcome_internal_"
  df_rm[[safe_outcome]] <- df_rm[[outcome]]

  # ---- 1.  Normality of residuals -------------------------------------------
  assumptions_rm <- data.frame(
    check = character(), result = character(),
    p_value = numeric(), decision = character(),
    stringsAsFactors = FALSE
  )

  rhs        <- paste(paste0("`", c(if (!is.null(group)) group, within), "`"),
                      collapse = " * ")
  lm_formula <- stats::as.formula(paste(safe_outcome, "~", rhs))
  lm_fit     <- tryCatch(stats::lm(lm_formula, data = df_rm), error = function(e) NULL)

  if (!is.null(lm_fit)) {
    nc <- .check_normality(lm_fit$residuals, NULL, normality_method,
                           alpha_normality, warnings_list)
    warnings_list <- nc$warnings_list
    nc$label      <- "Normality of residuals"
    assumptions_rm <- rbind(assumptions_rm,
                             .assumption_row(nc$label, nc$method,
                                             nc$p_value, nc$violated))
    if (nc$violated) {
      msg <- "Normality of residuals violated; interpret RM/mixed ANOVA results with caution."
      warnings_list <- c(warnings_list, msg)
      if (verbose) message("Warning: ", msg)
    }
  }

  # ---- 2.  rstatix::anova_test ----------------------------------------------
  anova_args <- list(
    data        = df_rm,
    dv          = safe_outcome,
    wid         = subject_id,
    within      = within,
    type        = type_sosquares,
    effect.size = rm_effect_size[1],
    detailed    = TRUE
  )
  if (test_family == "mixed_anova") anova_args[["between"]] <- group

  aov_res <- tryCatch(
    do.call(rstatix::anova_test, anova_args),
    error = function(e) stop("rstatix::anova_test failed: ", e$message)
  )

  aov_tbl <- as.data.frame(
    rstatix::get_anova_table(aov_res, correction = correction)
  )

  # Second effect-size pass when both "ges" and "pes" requested
  if (length(rm_effect_size) == 2) {
    second_es        <- rm_effect_size[2]
    anova_args2      <- anova_args
    anova_args2$effect.size <- second_es
    aov_res2 <- tryCatch(
      do.call(rstatix::anova_test, anova_args2),
      error = function(e) {
        msg <- paste0("Second effect size ('", second_es, "') pass failed: ", e$message)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message("Warning: ", msg)
        NULL
      }
    )
    if (!is.null(aov_res2)) {
      aov_tbl2 <- as.data.frame(
        rstatix::get_anova_table(aov_res2, correction = correction)
      )
      if (second_es %in% names(aov_tbl2) && !second_es %in% names(aov_tbl)) {
        idx              <- match(aov_tbl$Effect, aov_tbl2$Effect)
        aov_tbl[[second_es]] <- aov_tbl2[[second_es]][idx]
      }
    }
  }

  # ---- 3.  Sphericity (Mauchly) ---------------------------------------------
  if (!is.null(aov_res$`Mauchly's Test for Sphericity`)) {
    mauchly_tbl <- as.data.frame(aov_res$`Mauchly's Test for Sphericity`)
    for (i in seq_len(nrow(mauchly_tbl))) {
      sph_p   <- mauchly_tbl$p[i]
      sph_eff <- as.character(mauchly_tbl$Effect[i])
      sph_dec <- ifelse(!is.na(sph_p) & sph_p < alpha_sphericity, "VIOLATED", "OK")
      assumptions_rm <- rbind(assumptions_rm, data.frame(
        check    = paste0("Sphericity (Mauchly): ", sph_eff),
        result   = sprintf("W = %.4f", mauchly_tbl$W[i]),
        p_value  = sph_p,
        decision = sph_dec,
        stringsAsFactors = FALSE
      ))
      if (sph_dec == "VIOLATED") {
        msg <- sprintf("Sphericity violated for '%s' (Mauchly p = %.4f). GG correction applied.",
                       sph_eff, sph_p)
        warnings_list <- c(warnings_list, msg)
        if (verbose) message("Warning: ", msg)
      }
    }
    if (!is.null(aov_res$`Sphericity Corrections`)) {
      sph_cor <- as.data.frame(aov_res$`Sphericity Corrections`)
      for (i in seq_len(nrow(sph_cor))) {
        assumptions_rm <- rbind(assumptions_rm, data.frame(
          check    = paste0("GG epsilon: ", sph_cor$Effect[i]),
          result   = sprintf("eps = %.4f", sph_cor$GGe[i]),
          p_value  = sph_cor$`p[GG]`[i],
          decision = "INFO",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # ---- 4.  Effect size (primary / last effect) ------------------------------
  primary_row    <- aov_tbl[nrow(aov_tbl), ]
  primary_es_col <- rm_effect_size[1]

  es_rm <- if (primary_es_col %in% names(aov_tbl) &&
               !is.na(primary_row[[primary_es_col]])) {
    es_val   <- as.numeric(primary_row[[primary_es_col]])
    es_label <- if (primary_es_col == "ges") "Generalized eta-squared (ges)"
                else                         "Partial eta-squared (pes)"
    mag      <- if (es_val >= 0.14) "large" else if (es_val >= 0.06) "medium" else "small"
    list(
      estimate  = es_val,
      ci_low    = NA_real_,
      ci_high   = NA_real_,
      magnitude = mag,
      metric    = es_label,
      interpretation = sprintf("Effect size: %s (%s = %.3f)", mag, primary_es_col, es_val)
    )
  } else NULL

  # ---- 5.  Post-hoc ----------------------------------------------------------
  posthoc_rm  <- NULL
  sig_effects <- aov_tbl$Effect[!is.na(aov_tbl$p) & aov_tbl$p < test_alpha]

  if (perform_posthoc && requireNamespace("rstatix", quietly = TRUE) &&
      length(sig_effects) > 0) {

    ph_list <- list()
    for (eff in sig_effects) {
      eff_factors    <- trimws(strsplit(eff, ":")[[1]])
      within_in_eff  <- intersect(eff_factors, within)
      between_in_eff <- if (!is.null(group)) intersect(eff_factors, group) else character(0)
      if (length(within_in_eff) == 0 && length(between_in_eff) == 0) next

      ph_one <- tryCatch({
        ph_formula <- stats::as.formula(
          paste(safe_outcome, "~",
                paste0("`", if (length(within_in_eff) > 0)
                               within_in_eff[1] else between_in_eff[1], "`"))
        )
        rstatix::pairwise_t_test(
          data            = df_rm,
          formula         = ph_formula,
          paired          = length(within_in_eff) > 0,
          p.adjust.method = p_adjust_method,
          detailed        = TRUE
        ) |> as.data.frame()
      }, error = function(e) {
        msg <- paste0("Post-hoc failed for effect '", eff, "': ", e$message)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        NULL
      })

      if (!is.null(ph_one) && nrow(ph_one) > 0) {
        ph_one$effect       <- eff
        ph_one$posthoc_test <- if (length(within_in_eff) > 0)
          "Pairwise t-test (paired)" else "Pairwise t-test (independent)"
        ph_list[[eff]] <- ph_one
      }
    }

    if (length(ph_list) > 0) {
      all_cols <- unique(unlist(lapply(ph_list, names)))
      ph_list  <- lapply(ph_list, function(df) {
        for (cn in setdiff(all_cols, names(df))) df[[cn]] <- NA
        df[, all_cols, drop = FALSE]
      })
      posthoc_rm           <- do.call(rbind, ph_list)
      rownames(posthoc_rm) <- NULL

      if ("p.adj" %in% names(posthoc_rm)) {
        posthoc_rm$significant <- posthoc_rm$p.adj < test_alpha
        posthoc_rm$p.adj.signif <- ifelse(
          posthoc_rm$p.adj < test_alpha,
          ifelse(posthoc_rm$p.adj < 0.001, "***",
          ifelse(posthoc_rm$p.adj < 0.01,  "**",
          ifelse(posthoc_rm$p.adj < 0.05,  "*",
                 sprintf("* (p<%.2f)", test_alpha)))),
          "ns"
        )
        posthoc_rm$p.adj.format <- ifelse(
          posthoc_rm$p.adj < 0.001, "p < 0.001",
          sprintf("p = %.3f", posthoc_rm$p.adj)
        )
      }
      if ("p" %in% names(posthoc_rm)) {
        posthoc_rm$p.format <- ifelse(
          is.na(posthoc_rm$p), NA_character_,
          ifelse(posthoc_rm$p < 0.001, "p < 0.001",
                 sprintf("p = %.3f", posthoc_rm$p))
        )
      }

      # marginal means for interpretation
      mean_tbl <- stats::setNames(
        as.numeric(tapply(df_rm[[safe_outcome]],
                          df_rm[[within[1]]], mean, na.rm = TRUE)),
        levels(as.factor(df_rm[[within[1]]]))
      )
      posthoc_rm$mean1 <- mean_tbl[match(posthoc_rm$group1, names(mean_tbl))]
      posthoc_rm$mean2 <- mean_tbl[match(posthoc_rm$group2, names(mean_tbl))]

      posthoc_rm$interpretation <- ifelse(
        !posthoc_rm$significant %in% TRUE,
        "No significant difference",
        ifelse(!is.na(posthoc_rm$mean1) & !is.na(posthoc_rm$mean2) &
                 posthoc_rm$mean1 < posthoc_rm$mean2,
               paste0(posthoc_rm$group1, " has lower mean than ",  posthoc_rm$group2),
               paste0(posthoc_rm$group1, " has higher mean than ", posthoc_rm$group2))
      )
    }
  }

  # ---- 6.  Data summary ------------------------------------------------------
  summary_factors <- if (test_family == "mixed_anova") c(group, within) else within

  ds_rm <- tryCatch({
    agg <- stats::aggregate(
      stats::as.formula(
        paste(safe_outcome, "~",
              paste(paste0("`", summary_factors, "`"), collapse = " + "))
      ),
      data = df_rm,
      FUN  = function(v) c(n      = length(v),
                            mean   = mean(v),
                            sd     = stats::sd(v),
                            se     = stats::sd(v) / sqrt(length(v)),
                            median = stats::median(v),
                            q25    = stats::quantile(v, 0.25, names = FALSE),
                            q75    = stats::quantile(v, 0.75, names = FALSE),
                            min    = min(v),
                            max    = max(v))
    )
    stats_mat        <- as.data.frame(agg[[safe_outcome]])
    names(stats_mat) <- c("n", "mean", "sd", "se", "median",
                          "q25", "q75", "min", "max")
    cbind(agg[, summary_factors, drop = FALSE], stats_mat)
  }, error = function(e) {
    warning("Could not build data_summary for RM/mixed ANOVA: ", e$message)
    NULL
  })

  # ---- 7.  Verbose output ----------------------------------------------------
  if (verbose) {
    message("\n--- RM/Mixed ANOVA Results ---")
    message(sprintf("Design: %s", test_family))
    print(aov_tbl)
    if (!is.null(posthoc_rm) && nrow(posthoc_rm) > 0) {
      message(sprintf("\n--- Post-Hoc (Pairwise t-test, adj: %s) ---", p_adjust_method))
      sig_ph <- posthoc_rm[posthoc_rm$significant %in% TRUE, ]
      if (nrow(sig_ph) > 0) print(sig_ph) else message("No significant pairwise differences.")
    }
  }

  # ---- 8.  test_result stub (compatibility with print/summary) ---------------
  test_result_rm <- list(
    statistic   = stats::setNames(as.numeric(primary_row$F), "F value"),
    parameter   = stats::setNames(
      c(as.numeric(primary_row$DFn), as.numeric(primary_row$DFd)),
      c("DFn", "DFd")
    ),
    p.value     = as.numeric(primary_row$p),
    method      = test_family,
    anova_table = aov_tbl
  )

  # ---- 9.  Return ------------------------------------------------------------
  result_rm <- list(
    test_used      = switch(test_family,
      one_way_rm_anova = "One-way Repeated Measures ANOVA",
      two_way_rm_anova = "Two-way Repeated Measures ANOVA",
      mixed_anova      = "Two-way Mixed ANOVA"
    ),
    test_result    = test_result_rm,
    effect_size    = es_rm,
    posthoc_result = posthoc_rm,
    assumptions    = assumptions_rm,
    parametric     = TRUE,
    data_summary   = ds_rm,
    anova_table    = aov_tbl,
    sphericity     = aov_res$`Mauchly's Test for Sphericity`,
    sphericity_corrections = aov_res$`Sphericity Corrections`,
    outcome        = outcome,
    group_var      = group,
    within_vars    = within,
    subject_id_var = subject_id,
    paired         = FALSE,
    test_type      = parameters$test_type,
    test_alpha     = test_alpha,
    subgroup_analysis = NULL,
    warnings       = warnings_list,
    parameters     = parameters,
    raw_data       = df_rm
  )
  class(result_rm) <- "run_diff"
  result_rm
}


# =============================================================================
#  SECTION 5 — Effect-size helper (unexported)
# =============================================================================

# -----------------------------------------------------------------------------
#' Compute the appropriate effect size for a classical run_diff result
#'
#' @return A list with \code{estimate}, \code{ci_low}, \code{ci_high},
#'   \code{magnitude}, \code{metric}, and \code{interpretation}.
#' @keywords internal
.compute_effect_size <- function(y, grp, groups,
                                 test_family, parametric_used,
                                 variance_violated,
                                 fitted_model,
                                 warnings_list) {

  if (!requireNamespace("effectsize", quietly = TRUE))
    return(list(result = NULL, warnings_list = warnings_list))

  es_result <- tryCatch({
    if (test_family %in% c("independent_two_sample", "paired_two_sample")) {
      is_paired_fam <- (test_family == "paired_two_sample")
      if (parametric_used) {
        es     <- effectsize::cohens_d(
                    x         = y[grp == groups[1]],
                    y         = y[grp == groups[2]],
                    pooled_sd = !variance_violated,
                    paired    = is_paired_fam,
                    verbose   = FALSE
                  )
        interp <- effectsize::interpret(es, rules = "cohen1988")
        list(estimate  = es$Cohens_d,
             ci_low    = es$CI_low,
             ci_high   = es$CI_high,
             magnitude = as.character(interp$Interpretation),
             metric    = if (is_paired_fam) "Cohen's d (paired)" else "Cohen's d",
             interpretation = sprintf("Effect size: %s (%.3f)",
                                      as.character(interp$Interpretation), es$Cohens_d))
      } else {
        es     <- effectsize::rank_biserial(
                    x       = y[grp == groups[1]],
                    y       = y[grp == groups[2]],
                    paired  = is_paired_fam,
                    verbose = FALSE
                  )
        interp <- effectsize::interpret(es, rules = "funder2019")
        list(estimate  = es$r_rank_biserial,
             ci_low    = es$CI_low,
             ci_high   = es$CI_high,
             magnitude = as.character(interp$Interpretation),
             metric    = if (is_paired_fam) "Rank-biserial correlation (paired)"
                         else "Rank-biserial correlation",
             interpretation = sprintf("Effect size: %s (%.3f)",
                                      as.character(interp$Interpretation),
                                      es$r_rank_biserial))
      }
    } else {
      # independent_anova
      if (parametric_used) {
        es_model <- if (is.null(fitted_model)) {
          warnings_list <<- c(warnings_list,
            "Eta-squared calculated from standard ANOVA (Welch's ANOVA used for test).")
          stats::aov(y ~ grp)
        } else fitted_model
        es     <- effectsize::eta_squared(es_model, ci = 0.95, verbose = FALSE)
        interp <- effectsize::interpret(es, rules = "field2013")
        list(estimate  = es$Eta2,
             ci_low    = es$CI_low,
             ci_high   = es$CI_high,
             magnitude = as.character(interp$Interpretation),
             metric    = "Eta-squared",
             interpretation = sprintf("Effect size: %s (%.3f)",
                                      as.character(interp$Interpretation), es$Eta2))
      } else {
        es     <- effectsize::rank_epsilon_squared(x = y, groups = grp,
                                                   ci = 0.95, verbose = FALSE)
        interp <- effectsize::interpret(es, rules = "field2013")
        list(estimate  = es$rank_epsilon_squared,
             ci_low    = if (!is.null(es$CI_low))  es$CI_low  else NA_real_,
             ci_high   = if (!is.null(es$CI_high)) es$CI_high else NA_real_,
             magnitude = as.character(interp$Interpretation),
             metric    = "Rank epsilon-squared",
             interpretation = sprintf("Effect size: %s (%.3f)",
                                      as.character(interp$Interpretation),
                                      es$rank_epsilon_squared))
      }
    }
  }, error = function(e) {
    msg <- paste("Could not calculate effect size:", e$message)
    warnings_list <<- c(warnings_list, msg)
    NULL
  })

  list(result = es_result, warnings_list = warnings_list)
}


# =============================================================================
#  SECTION 6 — Public API
# =============================================================================

#' Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @title Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @description
#' Performs an automatic comparison of a numeric outcome variable across two or
#' more groups. Intelligently selects the appropriate statistical test
#' (parametric or non-parametric) based on checks for normality (of data or
#' residuals) and homogeneity of variances, calculates effect sizes, and
#' performs post-hoc tests when appropriate.
#'
#' When \code{outcome} is a character vector of length > 1 and
#' \code{summary_table = TRUE}, an additional \code{$summary_table} element is
#' returned alongside the per-outcome list.
#'
#' When \code{subgroup} is also specified and \code{summary_table = TRUE},
#' additional named elements \code{$summary_table_<subgroup_var>} (one per
#' subgroup variable) are appended; each is a named list of per-level summary
#' tables, e.g. \code{$summary_table_Gender$Male}.
#'
#' @details
#' \strong{Test Family Selection:}
#' If \code{paired = FALSE}, independent groups are compared. With 2 groups an
#' independent or paired two-sample test is chosen; with > 2 groups an
#' ANOVA-family test is used.
#'
#' \strong{Assumption Checking (\code{test_type = "auto"}):}
#' \itemize{
#'   \item \strong{Normality:} Shapiro-Wilk (n < 50) or Lilliefors (n ≥ 50)
#'     under \code{"auto"}. Per-group for 2-group comparisons; on model
#'     residuals for > 2 groups.
#'   \item \strong{Homogeneity of Variance:} Levene's test (\code{car});
#'     falls back to Bartlett's test if \code{car} is unavailable.
#' }
#'
#' \strong{Test Selection Logic (\code{test_type = "auto"}):}
#' \itemize{
#'   \item Independent two-sample: Mann-Whitney U (normality violated) →
#'     Welch's t (variance violated) → Student's t.
#'   \item Paired two-sample: Wilcoxon signed-rank (normality violated) →
#'     Paired t.
#'   \item Independent > 2 groups: Kruskal-Wallis (normality violated) →
#'     Welch's ANOVA (variance violated) → One-way ANOVA.
#' }
#'
#' \strong{Post-Hoc Tests (> 2 groups, significant omnibus):}
#' Tukey HSD (ANOVA), Games-Howell (Welch's ANOVA), Dunn's test
#' (Kruskal-Wallis), Pairwise t-test (RM / Mixed ANOVA).
#'
#' \strong{Effect Sizes:}
#' Cohen's d / paired Cohen's d (t-tests), rank-biserial correlation
#' (Wilcoxon/Mann-Whitney), eta-squared (ANOVA), epsilon-squared
#' (Kruskal-Wallis), generalized or partial eta-squared (RM / Mixed ANOVA).
#'
#' \strong{Non-Parametric Interpretation (Stochastic Dominance):}
#' For non-parametric tests the \code{interpretation} column reads
#' "X tends to have larger/smaller values than Y".
#'
#' @param x A data frame, or an object of class \code{"run_DIpreprocess"}.
#'   If a data frame, it should contain only the numeric features/metabolites
#'   to be analysed.
#' @param metadata A data frame containing sample-level metadata (e.g.,
#'   grouping, subgroups, subject IDs). Required if \code{x} is a data frame.
#' @param outcome A character string or character vector naming the numeric
#'   outcome variable(s) found in \code{x}.
#' @param group A character string naming the grouping variable found in
#'   \code{metadata}. May be \code{NULL} for pure within-subjects designs.
#' @param filter Optional character vector of levels in \code{group} to exclude
#'   from the analysis. If \code{x} is a \code{run_DIpreprocess} object,
#'   defaults to the identified QC types.
#' @param within Character vector of within-subject factor column name(s)
#'   found in \code{metadata}.
#'   When supplied, \code{subject_id} must also be provided. Triggers
#'   repeated-measures or mixed-ANOVA logic.
#' @param subject_id Character string naming the subject identifier column
#'   found in \code{metadata}.
#'   Required when \code{within} is non-\code{NULL}.
#' @param paired Logical; whether observations are paired. Default \code{FALSE}.
#' @param test_type One of \code{"auto"} (default), \code{"parametric"}, or
#'   \code{"nonparametric"}.
#' @param normality_method One of \code{"auto"} (default), \code{"shapiro"},
#'   or \code{"lilliefors"}.
#' @param alpha_normality Significance level for normality tests. Default
#'   \code{0.05}.
#' @param alpha_variance Significance level for the variance-homogeneity test.
#'   Default \code{0.05}.
#' @param alpha_sphericity Significance level for Mauchly's test. Default
#'   \code{0.05}.
#' @param type_sosquares Integer (1, 2, or 3) specifying the sums-of-squares
#'   type for RM / mixed ANOVA. Default \code{2}.
#' @param rm_effect_size Character vector; one or both of \code{"ges"}
#'   (generalized eta-squared, default) and \code{"pes"} (partial
#'   eta-squared). Generalized eta-squared is the default because it is
#'   directly comparable across designs that differ in the number of within-
#'   and between-subject factors. Use \code{"pes"} when matching output from
#'   software such as SPSS. Specify \code{c("ges", "pes")} to compute both;
#'   the first listed metric is treated as primary. Applies only to RM / mixed
#'   ANOVA designs.
#' @param correction Sphericity correction: \code{"auto"} (default; applies
#'   Greenhouse-Geisser when sphericity is violated), \code{"GG"},
#'   \code{"HF"}, or \code{"none"}. Applies only to RM / mixed ANOVA.
#' @param test_alpha Significance level for the final test. Default \code{0.05}.
#' @param min_n_threshold Minimum per-group sample size. Default \code{3}.
#' @param calculate_effect_size Logical. Default \code{TRUE}.
#' @param perform_posthoc Logical. Default \code{TRUE}.
#' @param p_adjust_method P-value adjustment method for post-hoc comparisons.
#'   Default \code{"BH"}.
#' @param group_order Optional character vector specifying group level order.
#' @param subgroup Optional character vector of column names in \code{metadata}
#'   for stratified analysis.
#' @param summary_table Logical; when \code{outcome} has length > 1 and
#'   \code{TRUE}, appends \code{$summary_table} (and subgroup summary tables)
#'   to the returned list. Default \code{FALSE}.
#' @param verbose Logical. Default \code{TRUE}.
#' @param num_cores Integer or \code{"max"}. Default \code{1}.
#'
#' @return When \code{outcome} is a single string, a list of class
#'   \code{"run_diff"} containing: \code{test_used}, \code{test_result},
#'   \code{effect_size}, \code{posthoc_result}, \code{assumptions},
#'   \code{parametric}, \code{data_summary}, \code{outcome},
#'   \code{group_var}, \code{within_vars} (RM only), \code{subject_id_var}
#'   (RM only), \code{anova_table} (RM only), \code{sphericity} (RM only),
#'   \code{sphericity_corrections} (RM only), \code{paired}, \code{test_type},
#'   \code{test_alpha}, \code{subgroup_analysis}, \code{warnings},
#'   \code{parameters}, \code{raw_data}.
#'
#'   When \code{outcome} is length > 1, a named list of \code{"run_diff"}
#'   objects, optionally with \code{$summary_table} and
#'   \code{$summary_table_<sg>} elements.
#'
#' @importFrom car leveneTest
#' @importFrom effectsize cohens_d rank_biserial eta_squared epsilon_squared interpret
#' @importFrom rstatix tukey_hsd games_howell_test dunn_test
#' @importFrom moments skewness kurtosis
#' @importFrom nortest lillie.test
#' @importFrom tibble tibble
#' @importFrom stats aov as.formula bartlett.test complete.cases kruskal.test lm
#' @importFrom stats median oneway.test quantile sd shapiro.test t.test wilcox.test
#' @importFrom utils install.packages
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{run_DIpreprocess}}, \code{\link{run_summarytable}}, \code{\link{get_volcanodata}}
#' 
#' @examples
#' \dontrun{
#' # --- 3+ Groups (auto) ---
#' res_iris <- run_diff(iris, "Sepal.Length", "Species")
#' print(res_iris); summary(res_iris)
#'
#' # --- Two-sample independent ---
#' mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
#' run_diff(mtcars, "mpg", "am")
#'
#' # --- Paired ---
#' run_diff(sleep, "extra", "group", paired = TRUE)
#'
#' # --- Multi-outcome with summary table ---
#' outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' res_multi <- run_diff(iris, outcomes, "Species",
#'                       summary_table = TRUE, verbose = FALSE)
#' res_multi$summary_table
#'
#' # --- One-way RM ANOVA ---
#' run_diff(df_long, "HeartRate", within = "Time", subject_id = "ID")
#'
#' # --- Mixed ANOVA ---
#' run_diff(df_long, "HeartRate", group = "Group",
#'          within = "Time", subject_id = "ID")
#' }
#'
#' @export
run_diff <- function(
    x,
    metadata              = NULL,
    outcome,
    group                 = NULL,
    filter                = NULL,
    within                = NULL,
    subject_id            = NULL,
    paired                = FALSE,
    test_type             = c("auto", "parametric", "nonparametric"),
    normality_method      = c("auto", "shapiro", "lilliefors"),
    alpha_normality       = 0.05,
    alpha_variance        = 0.05,
    alpha_sphericity      = 0.05,
    type_sosquares        = 2,
    rm_effect_size        = "ges",
    correction            = "auto",
    test_alpha            = 0.05,
    min_n_threshold       = 3,
    calculate_effect_size = TRUE,
    perform_posthoc       = TRUE,
    p_adjust_method       = "BH",
    group_order           = NULL,
    subgroup              = NULL,
    summary_table         = FALSE,
    verbose               = TRUE,
    num_cores             = 1
) {

  # --- run_DIpreprocess integration ------------------------------------------
  if (inherits(x, "run_DIpreprocess")) {
    if (is.null(metadata)) {
      metadata <- if (!is.null(x$metadata_merged)) x$metadata_merged else x$metadata
    }
    x_data  <- if (!is.null(x$data_nonpls_merged)) x$data_nonpls_merged else x$data_nonpls
    
    if (missing(outcome) || is.null(outcome)) outcome <- names(x_data)
    if (is.null(group)) group <- x$parameters$group_col

    if (is.null(filter)) {
      filter <- eval(x$parameters$qc_types)
    }
    x <- x_data
  }

  # --- Default outcome if not specified ---
  if (missing(outcome) || is.null(outcome)) {
    outcome <- names(x)
  }

  # --- Basic Alignment Validation --------------------------------------------
  if (!is.data.frame(x)) stop("'x' must be a data frame or matrix of features.")
  if (is.null(metadata)) stop("'metadata' must be provided.")
  if (nrow(x) != nrow(metadata)) {
    stop(sprintf("Row mismatch: 'x' has %d rows, but 'metadata' has %d rows.", 
                 nrow(x), nrow(metadata)))
  }

  # --- Validate x contains only numeric features -----------------------------
  non_numeric_cols <- names(x)[!vapply(x, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    stop(sprintf("All columns in 'x' must be numeric. Non-numeric column(s) found: %s",
                 paste(non_numeric_cols, collapse = ", ")))
  }

  # --- Filtering -------------------------------------------------------------
  if (!is.null(group) && !is.null(filter)) {
    if (!is.character(filter)) stop("'filter' must be a character vector.")
    if (!group %in% names(metadata)) {
      stop(paste0("Group variable '", group, "' not found in 'metadata'."))
    }

    grp_raw <- as.character(metadata[[group]])
    keep_rows <- !grp_raw %in% filter

    if (sum(keep_rows) == 0L) stop("After applying 'filter', no rows remain.")
    if (sum(!keep_rows) > 0L) {
      x <- x[keep_rows, , drop = FALSE]
      metadata <- metadata[keep_rows, , drop = FALSE]
      if (verbose) message(sprintf("Filtered out %d rows based on 'filter' levels.", sum(!keep_rows)))
    }
  }

  # ===========================================================================
  #  Multi-outcome dispatch
  # ===========================================================================
  if (length(outcome) > 1) {
    if (!is.logical(summary_table) || length(summary_table) != 1L)
      stop("'summary_table' must be a single logical value (TRUE or FALSE).")

    results <- lapply(outcome, function(oc) {
      tryCatch(
        run_diff(
          x = x, metadata = metadata, outcome = oc, group = group, 
          filter = filter, within = within,
          subject_id = subject_id, paired = paired,
          test_type = test_type, normality_method = normality_method,
          alpha_normality = alpha_normality, alpha_variance = alpha_variance,
          alpha_sphericity = alpha_sphericity,
          type_sosquares = type_sosquares, rm_effect_size = rm_effect_size,
          correction = correction, test_alpha = test_alpha,
          min_n_threshold = min_n_threshold,
          calculate_effect_size = calculate_effect_size,
          perform_posthoc = perform_posthoc, p_adjust_method = p_adjust_method,
          group_order = group_order, subgroup = subgroup,
          summary_table = FALSE, verbose = verbose, num_cores = num_cores
        ),
        error = function(e) {
          warning(paste0("run_diff failed for outcome '", oc, "': ", e$message),
                  call. = FALSE)
          NULL
        }
      )
    })
    names(results) <- outcome

    if (summary_table) {
      results[["summary_table"]] <- .build_summary_table(results, outcome, test_alpha)

      if (!is.null(subgroup)) {
        for (sg in subgroup) {
          sg_levels <- if (is.factor(metadata[[sg]])) {
            levels(droplevels(as.factor(metadata[[sg]])))
          } else {
            sort(unique(na.omit(as.character(metadata[[sg]]))))
          }
          sg_tables <- lapply(
            stats::setNames(sg_levels, sg_levels),
            function(lvl) {
              level_results <- stats::setNames(
                lapply(outcome, function(oc) {
                  res <- results[[oc]]
                  if (is.null(res)) return(NULL)
                  res$subgroup_analysis[[sg]][[lvl]]
                }),
                outcome
              )
              .build_summary_table(level_results, outcome, test_alpha)
            }
          )
          results[[paste0("summary_table_", sg)]] <- sg_tables
        }
      }
    }
    class(results) <- c("run_diff", "list")
    return(results)
  }

  # ===========================================================================
  #  Single-outcome path
  # ===========================================================================

  # --- 1. Input validation & preparation -------------------------------------
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  if (!is.logical(summary_table) || length(summary_table) != 1L)
    stop("'summary_table' must be a single logical value (TRUE or FALSE).")

  # Validate subgroup
  if (!is.null(subgroup)) {
    if (!is.character(subgroup))
      stop("'subgroup' must be a character vector of column names.")
    missing_sg <- setdiff(subgroup, names(metadata))
    if (length(missing_sg) > 0)
      stop(paste("Subgroup column(s) not found in 'metadata':",
                 paste(missing_sg, collapse = ", ")))
    non_cat <- subgroup[sapply(subgroup, function(s)
      is.numeric(metadata[[s]]) && !is.factor(metadata[[s]]))]
    if (length(non_cat) > 0)
      stop(paste("Subgroup column(s) must be categorical:",
                 paste(non_cat, collapse = ", ")))
  }

  # Validate RM parameters
  if (!type_sosquares %in% c(1L, 2L, 3L))
    stop("'type_sosquares' must be 1, 2, or 3.")
  type_sosquares <- as.integer(type_sosquares)

  if (!all(rm_effect_size %in% c("ges", "pes")))
    stop("'rm_effect_size' must be \"ges\", \"pes\", or c(\"ges\", \"pes\").")
  rm_effect_size <- intersect(c("ges", "pes"), rm_effect_size)  # canonical order

  correction <- match.arg(correction, c("auto", "GG", "HF", "none"))

  is_rm_design    <- !is.null(within)
  is_mixed_design <- is_rm_design && !is.null(group)   # kept for potential future use

  if (is_rm_design) {
    if (!requireNamespace("rstatix", quietly = TRUE))
      stop("Package 'rstatix' is required for RM/mixed ANOVA. ",
           "Install with: install.packages('rstatix')")
    if (is.null(subject_id))
      stop("'subject_id' must be specified for RM/mixed ANOVA designs.")
    if (!is.character(within) || length(within) < 1)
      stop("'within' must be a character vector of within-subject column name(s).")
    missing_w <- setdiff(within, names(metadata))
    if (length(missing_w) > 0)
      stop(paste("Within-subject column(s) not found in 'metadata':",
                 paste(missing_w, collapse = ", ")))
    if (!subject_id %in% names(metadata))
      stop(paste0("subject_id column '", subject_id, "' not found in 'metadata'."))
  }

  # Resolve num_cores
  avail_cores <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.na(avail_cores)) avail_cores <- 1L
  if (identical(num_cores, "max")) {
    num_cores <- max(1L, avail_cores - 2L)
    if (verbose)
      message(sprintf("Parallel: using 'max' setting (%d cores, 2 reserved).", num_cores))
  } else if (is.numeric(num_cores)) {
    if (num_cores < 1 || num_cores > avail_cores || num_cores %% 1 != 0)
      stop(sprintf("'num_cores' must be an integer between 1 and %d, or \"max\".",
                   avail_cores))
    num_cores <- as.integer(num_cores)
    if (num_cores > 1 && verbose)
      message(sprintf("Parallel: using %d cores.", num_cores))
  } else {
    stop("Invalid 'num_cores'. Must be an integer or \"max\".")
  }

  # Validate outcome / group
  if (!is.data.frame(x)) stop("'x' must be a data frame.")
  if (!outcome %in% names(x))
    stop(paste0("Outcome variable '", outcome, "' not found in 'x'."))
  if (!is.null(group) && !group %in% names(metadata))
    stop(paste0("Group variable '", group, "' not found in 'metadata'."))

  y <- x[[outcome]]
  if (!is.numeric(y)) stop("Outcome variable must be numeric.")

  if (!is.null(group)) {
    grp <- as.factor(metadata[[group]])
    if (!is.null(group_order)) {
      if (!all(group_order %in% levels(grp))) {
        stop("'group_order' contains levels not present in the grouping column.")
      }
      grp <- factor(grp, levels = group_order)
    }
    cc <- stats::complete.cases(y, grp)
    if (sum(!cc) > 0) {
      msg           <- sprintf("Removing %d rows with missing values.", sum(!cc))
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
    y   <- y[cc]
    grp <- if (length(cc) == length(grp)) droplevels(grp[cc]) else grp # handle subsetting
    
    # Sync metadata if needed for RM
    metadata_cc <- metadata[cc, , drop = FALSE]
    
    if (nlevels(grp) < 2)
      stop("Grouping factor must have at least 2 levels after removing missing values.")
    groups   <- levels(grp)
    n_groups <- length(groups)
    group_ns <- table(grp)
    min_n    <- min(group_ns)
  } else {
    if (!is_rm_design) {
       stop("A grouping variable must be provided via 'group' for independent comparisons.")
    }
    # Pure RM design — dummy grouping for internal consistency
    grp      <- factor(rep("all", length(y)))
    groups   <- levels(grp)
    n_groups <- 1L
    group_ns <- table(grp)
    min_n    <- length(y)
    metadata_cc <- metadata
  }

  if (verbose && !is_rm_design) {
    message("\n=== Auto Compare Analysis ===")
    message(sprintf("Outcome: %s | Group: %s | Paired: %s", outcome, group, paired))
    message(sprintf("Number of groups: %d", n_groups))
    message(sprintf("Sample sizes: %s",
                    paste(names(group_ns), "=", group_ns, collapse = ", ")))
    message(sprintf("Test Type: %s", test_type))
  }

  # --- 2. Test-family determination ------------------------------------------
  if (is_rm_design) {
    has_between <- !is.null(group) && 
                   group %in% names(metadata_cc) && 
                   length(unique(na.omit(as.character(metadata_cc[[group]])))) > 1
                   
    test_family <- if (length(within) >= 1 && has_between) {
      "mixed_anova"
    } else if (length(within) >= 2) {
      "two_way_rm_anova"
    } else {
      "one_way_rm_anova"
    }
  } else if (n_groups == 2) {
    test_family <- if (paired) "paired_two_sample" else "independent_two_sample"
  } else {
    test_family <- if (paired) "repeated_measures_anova" else "independent_anova"
  }

  # Build the parameters list used across RM/non-RM branches
  parameters <- list(
    test_type             = test_type,
    normality_method      = normality_method,
    alpha_normality       = alpha_normality,
    alpha_variance        = alpha_variance,
    alpha_sphericity      = alpha_sphericity,
    type_sosquares        = type_sosquares,
    rm_effect_size        = rm_effect_size,
    correction            = correction,
    test_alpha            = test_alpha,
    min_n_threshold       = min_n_threshold,
    calculate_effect_size = calculate_effect_size,
    perform_posthoc       = perform_posthoc,
    p_adjust_method       = p_adjust_method,
    group_order           = group_order,
    filter                = filter,
    subgroup              = subgroup,
    within                = within,
    subject_id            = subject_id,
    paired                = paired
  )

  # --- RM / Mixed ANOVA early-return -----------------------------------------
  if (test_family %in% c("one_way_rm_anova", "two_way_rm_anova", "mixed_anova")) {
    return(.run_rm_anova(
      x = x, metadata = metadata, outcome = outcome, group = group, within = within,
      subject_id = subject_id, test_family = test_family,
      normality_method = normality_method, alpha_normality = alpha_normality,
      alpha_sphericity = alpha_sphericity,
      type_sosquares = type_sosquares, rm_effect_size = rm_effect_size,
      correction = correction, test_alpha = test_alpha,
      perform_posthoc = perform_posthoc, p_adjust_method = p_adjust_method,
      verbose = verbose, warnings_list = warnings_list,
      parameters = parameters
    ))
  }

  # --- 3. Assumption checks (classical designs) ------------------------------
  assumptions <- data.frame(
    check = character(), result = character(),
    p_value = numeric(), decision = character(),
    stringsAsFactors = FALSE
  )

  assumptions <- rbind(assumptions, data.frame(
    check    = "Minimum sample size",
    result   = sprintf("Min n = %d", min_n),
    p_value  = NA_real_,
    decision = ifelse(min_n >= min_n_threshold, "PASS", "WARNING"),
    stringsAsFactors = FALSE
  ))
  if (min_n < min_n_threshold) {
    msg           <- sprintf("Warning: Min sample size (%d) below threshold (%d).",
                             min_n, min_n_threshold)
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }

  normality_violated <- FALSE
  variance_violated  <- FALSE

  if (test_type == "auto") {
    # Normality
    if (n_groups == 2) {
      norm_results <- .parallel_normality(
        groups, y, grp, normality_method, alpha_normality, num_cores
      )
      for (nr in norm_results) {
        if (inherits(nr, "try-error")) { if (verbose) warning("A parallel worker failed."); next }
        warnings_list <- c(warnings_list, nr$warnings_list)
        assumptions   <- rbind(assumptions,
                               .assumption_row(nr$label, nr$method, nr$p_value, nr$violated))
        if (nr$violated) normality_violated <- TRUE
      }
    } else {
      fit <- stats::lm(y ~ grp)
      nr  <- .check_normality(fit$residuals, NULL, normality_method,
                              alpha_normality, warnings_list)
      warnings_list <- nr$warnings_list
      nr$label      <- "Normality of residuals"
      assumptions   <- rbind(assumptions,
                             .assumption_row(nr$label, nr$method, nr$p_value, nr$violated))
      if (nr$violated) normality_violated <- TRUE
    }

    # Variance homogeneity
    if (!paired) {
      vc            <- .check_variance(y, grp, alpha_variance, verbose, warnings_list)
      warnings_list <- vc$warnings_list
      variance_violated <- vc$violated
      assumptions   <- rbind(assumptions, data.frame(
        check    = "Homogeneity of variance",
        result   = vc$test_name,
        p_value  = vc$p_value,
        decision = if (vc$violated) "VIOLATED" else "OK",
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- 4. Select and execute test --------------------------------------------
  test_exec <- switch(
    test_family,
    independent_two_sample = .run_independent_two_sample(
      y, groups, grp, normality_violated, variance_violated,
      test_type, warnings_list
    ),
    paired_two_sample = .run_paired_two_sample(
      y, groups, grp, normality_violated, test_type, warnings_list
    ),
    independent_anova = .run_independent_anova(
      y, grp, normality_violated, variance_violated, test_type, warnings_list
    ),
    stop(sprintf("Unsupported test family: '%s'.", test_family))
  )

  test_used     <- test_exec$test_used
  test_result   <- test_exec$test_result
  fitted_model  <- test_exec$fitted_model
  parametric_used <- test_exec$parametric
  warnings_list <- test_exec$warnings_list

  # --- 5. Effect size --------------------------------------------------------
  es_out        <- .compute_effect_size(
    y, grp, groups, test_family, parametric_used,
    variance_violated, fitted_model, warnings_list
  )
  effect_size_result <- es_out$result
  warnings_list      <- es_out$warnings_list

  # --- 6. Post-hoc -----------------------------------------------------------
  posthoc_result <- NULL

  if (perform_posthoc && n_groups > 2 &&
      !is.null(test_result) && test_result$p.value < test_alpha &&
      requireNamespace("rstatix", quietly = TRUE)) {

    posthoc_result <- tryCatch({
      td <- tibble::tibble(y = y, grp = grp)
      ph <- NULL

      if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
        ph <- rstatix::tukey_hsd(td, y ~ grp)
        ph <- .build_posthoc_df(ph, "Tukey HSD", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      } else if (test_used == "Welch's ANOVA") {
        ph <- rstatix::games_howell_test(td, y ~ grp)
        ph <- .build_posthoc_df(ph, "Games-Howell", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      } else if (test_used %in% c("Kruskal-Wallis test", "Kruskal-Wallis test (forced)")) {
        ph <- rstatix::dunn_test(td, y ~ grp, p.adjust.method = p_adjust_method)
        ph <- .build_posthoc_df(ph, "Dunn's test", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      }
      ph
    }, error = function(e) {
      msg           <- paste("Could not perform post-hoc test:", e$message)
      warnings_list <<- c(warnings_list, msg)
      if (verbose) message(msg)
      NULL
    })
  }

  # --- 7. Data summary -------------------------------------------------------
  dg <- droplevels(grp)
  data_summary <- data.frame(
    group  = levels(dg),
    n      = as.numeric(table(dg)),
    mean   = as.numeric(tapply(y, dg, mean)),
    sd     = as.numeric(tapply(y, dg, stats::sd)),
    se     = as.numeric(tapply(y, dg,
               function(v) stats::sd(v) / sqrt(length(v)))),
    median = as.numeric(tapply(y, dg, stats::median)),
    q25    = as.numeric(tapply(y, dg, function(v) stats::quantile(v, 0.25))),
    q75    = as.numeric(tapply(y, dg, function(v) stats::quantile(v, 0.75))),
    min    = as.numeric(tapply(y, dg, min)),
    max    = as.numeric(tapply(y, dg, max)),
    row.names = NULL
  )
  if (requireNamespace("moments", quietly = TRUE)) {
    data_summary$skewness <- as.numeric(tapply(y, dg, moments::skewness))
    data_summary$kurtosis <- as.numeric(tapply(y, dg, moments::kurtosis))
  }

  # --- 8. Verbose output -----------------------------------------------------
  if (verbose) {
    message("\n--- Test Results ---")
    message(sprintf("Test Selected: %s", test_used))
    message(sprintf("Test statistic: %.4f", test_result$statistic[[1]]))
    if (test_result$p.value < 0.001) message("P-value: < 0.001") else
      message(sprintf("P-value: %.4f", test_result$p.value))
    if (!is.null(effect_size_result))
      message(sprintf("%s", effect_size_result$interpretation))
    if (test_result$p.value < test_alpha)
      message(sprintf("Significant difference detected (p < %.2f) ***", test_alpha))
    else
      message(sprintf("No significant difference (p >= %.2f)", test_alpha))
    if (!is.null(posthoc_result)) {
      message(sprintf("\n--- Post-Hoc (%s) ---",
                      unique(posthoc_result$posthoc_test)))
      message(sprintf("Significant pairwise differences: %d / %d",
                      sum(posthoc_result$significant, na.rm = TRUE),
                      nrow(posthoc_result)))
    }
  }

  # --- 9. Subgroup analysis --------------------------------------------------
  subgroup_analysis <- NULL

  if (!is.null(subgroup)) {
    subgroup_analysis <- list()
    for (sg_var in subgroup) {
      sg_levels <- if (is.factor(metadata[[sg_var]])) {
        levels(droplevels(as.factor(metadata[[sg_var]])))
      } else {
        sort(unique(na.omit(as.character(metadata[[sg_var]]))))
      }
      sg_var_results <- list()
      for (lvl in sg_levels) {
        x_sub   <- x[as.character(metadata[[sg_var]]) == lvl, , drop = FALSE]
        meta_sub <- metadata[as.character(metadata[[sg_var]]) == lvl, , drop = FALSE]
        if (nrow(x_sub) == 0) {
          warning(paste0("Subgroup '", sg_var, "' level '", lvl,
                         "' has no observations. Skipping."))
          next
        }
        grp_sub <- droplevels(as.factor(meta_sub[[group]]))
        if (nlevels(grp_sub) < 2) {
          warning(paste0("Subgroup '", sg_var, "' = '", lvl,
                         "': fewer than 2 group levels. Skipping."))
          next
        }
        sg_result <- tryCatch(
          run_diff(
            x = x_sub, metadata = meta_sub, outcome = outcome, group = group, filter = filter,
            within = within, subject_id = subject_id, paired = paired,
            test_type = test_type, normality_method = normality_method,
            alpha_normality = alpha_normality, alpha_variance = alpha_variance,
            alpha_sphericity = alpha_sphericity,
            type_sosquares = type_sosquares, rm_effect_size = rm_effect_size,
            correction = correction, test_alpha = test_alpha,
            min_n_threshold = min_n_threshold,
            calculate_effect_size = calculate_effect_size,
            perform_posthoc = perform_posthoc, p_adjust_method = p_adjust_method,
            group_order = group_order, verbose = verbose, num_cores = num_cores,
            summary_table = FALSE, subgroup = NULL
          ),
          error = function(e) {
            warning(paste0("run_diff failed for subgroup '", sg_var, "' = '", lvl,
                           "': ", e$message), call. = FALSE)
            NULL
          }
        )
        if (!is.null(sg_result)) sg_var_results[[lvl]] <- sg_result
      }
      subgroup_analysis[[sg_var]] <- sg_var_results
    }
  }

  # --- 10. Return object -----------------------------------------------------
  result <- list(
    test_used         = test_used,
    test_result       = test_result,
    effect_size       = effect_size_result,
    posthoc_result    = posthoc_result,
    assumptions       = assumptions,
    parametric        = parametric_used,
    data_summary      = data_summary,
    outcome           = outcome,
    group_var         = group,
    paired            = paired,
    test_type         = test_type,
    test_alpha        = test_alpha,
    subgroup_analysis = subgroup_analysis,
    warnings          = warnings_list,
    parameters        = parameters,
    raw_data          = tibble::tibble(outcome = y, group = grp)
  )
  class(result) <- c("run_diff", "list")
  result
}


# =============================================================================
#  SECTION 7 — S3 methods
# =============================================================================

#' Print method for run_diff objects
#' @param x An object from \code{run_diff}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.run_diff <- function(x, ...) {
  cat("\n=== Automatic Statistical Comparison ===\n\n")
  cat(sprintf("Test Used: %s\n", x$test_used))
  cat(sprintf("Parametric: %s\n\n", x$parametric))
  cat("Test Results:\n")
  cat(sprintf("  Statistic: %.4f\n", x$test_result$statistic[[1]]))
  if (!is.null(x$test_result$p.value)) {
    if (x$test_result$p.value < 0.001) {
      cat("  P-value: < 0.001 ***\n")
    } else {
      cat(sprintf("  P-value: %.4f", x$test_result$p.value))
      if      (x$test_result$p.value < 0.01) cat(" **")
      else if (x$test_result$p.value < 0.05) cat(" *")
      cat("\n")
    }
  }
  if (!is.null(x$effect_size))
    cat(sprintf("  %s\n", x$effect_size$interpretation))
  cat("\n")

  if (!is.null(x$anova_table)) {
    cat("ANOVA Table (all effects):\n")
    print(x$anova_table, row.names = FALSE)
  } else if (!is.null(x$data_summary)) {
    cat("Data Summary:\n")
    show_cols <- intersect(c("group", "n", "mean", "sd", "median"),
                           names(x$data_summary))
    if (length(show_cols) >= 2)
      print(x$data_summary[, show_cols, drop = FALSE], row.names = FALSE)
    else
      print(x$data_summary, row.names = FALSE)
  }

  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\n\nPost-Hoc Tests (%s):\n",
                paste(unique(x$posthoc_result$posthoc_test), collapse = ", ")))
    sig_only <- x$posthoc_result[x$posthoc_result$significant %in% TRUE, , drop = FALSE]
    ph_show  <- intersect(c("group1", "group2", "effect", "p.adj", "p.adj.signif"),
                          names(sig_only))
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      print(sig_only[, ph_show, drop = FALSE], row.names = FALSE)
    } else {
      cat("No significant pairwise differences found.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat(sprintf("  - %s\n", w))
  }
}


#' Summary method for run_diff objects
#' @param object An object from \code{run_diff}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.run_diff <- function(object, ...) {
  cat("\n=== Detailed Statistical Comparison Summary ===\n\n")
  cat(sprintf("Outcome Variable:  %s\n", object$outcome))
  cat(sprintf("Grouping Variable: %s\n", object$group_var))
  cat(sprintf("Paired Analysis:   %s\n", object$paired))
  cat(sprintf("Test Strategy:     %s (alpha = %.2f)\n",
              object$test_type, object$test_alpha))
  cat(sprintf("Test Selected:     %s\n\n", object$test_used))
  if (!is.null(object$within_vars))
    cat(sprintf("Within-Subject Factor(s): %s\n",
                paste(object$within_vars, collapse = ", ")))
  if (!is.null(object$subject_id_var))
    cat(sprintf("Subject ID Variable:      %s\n", object$subject_id_var))

  cat("--- Descriptive Statistics ---\n")
  print(object$data_summary, row.names = FALSE)
  cat("\n--- Assumption Checks ---\n")
  print(object$assumptions, row.names = FALSE)

  if (!is.null(object$anova_table)) {
    cat("\n--- Full ANOVA Table (all effects) ---\n")
    print(object$anova_table, row.names = FALSE)
  }
  if (!is.null(object$sphericity)) {
    cat("\n--- Mauchly's Test for Sphericity ---\n")
    print(as.data.frame(object$sphericity), row.names = FALSE)
  }
  if (!is.null(object$sphericity_corrections)) {
    cat("\n--- Sphericity Corrections (GG / HF) ---\n")
    print(as.data.frame(object$sphericity_corrections), row.names = FALSE)
  }

  cat("\n--- Statistical Test Results ---\n")
  if (is.null(object$anova_table)) {
    print(object$test_result)
  } else {
    tr <- object$test_result
    cat(sprintf("Primary effect: F(%s, %s) = %.4f, p = %s\n",
                tr$parameter["DFn"], tr$parameter["DFd"],
                tr$statistic[[1]],
                if (tr$p.value < 0.001) "< 0.001"
                else sprintf("%.4f", tr$p.value)))
  }

  if (!is.null(object$effect_size)) {
    cat("\n--- Effect Size ---\n")
    cat(sprintf("Metric:    %s\n", object$effect_size$metric))
    cat(sprintf("Estimate:  %.3f", object$effect_size$estimate))
    if (!is.na(object$effect_size$ci_low))
      cat(sprintf(" [95%% CI: %.3f, %.3f]",
                  object$effect_size$ci_low, object$effect_size$ci_high))
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$effect_size$magnitude))
  }

  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons (%s) ---\n",
                unique(object$posthoc_result$posthoc_test)))
    cat(sprintf("P-value adjustment method: %s\n\n", object$parameters$p_adjust_method))
    cols_to_show <- intersect(
      c("group1", "group2", "n1", "n2", "statistic", "p",
        "p.adj", "p.adj.signif", "interpretation"),
      names(object$posthoc_result)
    )
    print(object$posthoc_result[, cols_to_show], row.names = FALSE)
  }

  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT difference detected (p < %.2f).\n", object$test_alpha))
    if (!is.null(object$effect_size))
      cat(sprintf("Effect size is %s (%s = %.3f).\n",
                  tolower(object$effect_size$magnitude),
                  object$effect_size$metric,
                  object$effect_size$estimate))
  } else {
    cat(sprintf("Result: NO significant difference (p >= %.2f).\n", object$test_alpha))
  }

  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) cat(sprintf("  - %s\n", w))
  }
  invisible(object)
}


# =============================================================================
#  SECTION 8 — .run_diff_single (internal; used by plot_diff)
# =============================================================================

#' Internal single-outcome run_diff for plot_diff
#'
#' A lighter version of \code{run_diff} that returns a rich list suitable for
#' downstream use by \code{plot_diff}. Not intended for direct use.
#'
#' @keywords internal
.run_diff_single <- function(
    x,
    outcome,
    group,
    paired                = FALSE,
    test_type             = c("auto", "parametric", "nonparametric"),
    normality_method      = c("auto", "shapiro", "lilliefors"),
    alpha_normality       = 0.05,
    alpha_variance        = 0.05,
    test_alpha            = 0.05,
    min_n_threshold       = 3,
    calculate_effect_size = TRUE,
    perform_posthoc       = TRUE,
    p_adjust_method       = "BH",
    group_order           = NULL,
    verbose               = FALSE,
    num_cores             = 1L
) {
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  avail_cores <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1L
  num_cores <- if (identical(num_cores, "max")) max(1L, avail_cores - 2L) else as.integer(num_cores)

  if (!is.data.frame(x))     stop("'x' must be a data frame.")
  if (!outcome %in% names(x)) stop(sprintf("Outcome '%s' not found.", outcome))
  if (!group   %in% names(x)) stop(sprintf("Group '%s' not found.", group))

  y   <- x[[outcome]]
  grp <- as.factor(x[[group]])
  if (!is.numeric(y)) stop("Outcome variable must be numeric.")

  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(grp)))
      stop("'group_order' contains levels not present in the data.")
    grp <- factor(grp, levels = group_order)
  }

  cc <- stats::complete.cases(y, grp)
  if (sum(!cc) > 0) {
    msg           <- sprintf("Removing %d rows with missing values.", sum(!cc))
    warnings_list <- c(warnings_list, msg)
  }
  y   <- y[cc]
  grp <- droplevels(grp[cc])

  if (nlevels(grp) < 2)
    stop("Grouping factor must have at least 2 levels after removing missing values.")

  groups   <- levels(grp)
  n_groups <- length(groups)
  group_ns <- table(grp)
  min_n    <- min(group_ns)

  if (min_n < min_n_threshold)
    warnings_list <- c(warnings_list,
                       sprintf("Min sample size (%d) below threshold (%d).",
                               min_n, min_n_threshold))

  if (n_groups == 2) {
    test_family <- if (paired) "paired_two_sample" else "independent_two_sample"
  } else {
    if (paired) stop("Repeated measures ANOVA is not yet implemented in .run_diff_single.")
    test_family <- "independent_anova"
  }
  if (paired && n_groups == 2 && group_ns[1] != group_ns[2])
    stop("Paired tests require equal group sizes after removing missing values.")

  # --- Assumption checks -----------------------------------------------------
  normality_violated <- FALSE
  variance_violated  <- FALSE
  norm_results_raw   <- list()
  variance_p         <- NA_real_
  variance_method    <- NA_character_

  if (test_type == "auto") {
    if (n_groups == 2) {
      norm_results_raw <- .parallel_normality(
        groups, y, grp, normality_method, alpha_normality, num_cores
      )
    } else {
      fit              <- stats::lm(y ~ grp)
      nr               <- .check_normality(fit$residuals, "residuals",
                                           normality_method, alpha_normality)
      norm_results_raw <- list(nr)
    }
    for (nr in norm_results_raw) {
      if (inherits(nr, "try-error")) next
      warnings_list <- c(warnings_list, nr$warnings_list)
      if (isTRUE(nr$violated)) normality_violated <- TRUE
    }

    if (!paired) {
      vc            <- .check_variance(y, grp, alpha_variance, verbose, warnings_list)
      warnings_list <- vc$warnings_list
      variance_p    <- vc$p_value
      variance_method <- vc$test_name
      variance_violated <- vc$violated
    }
  }

  # --- Test execution --------------------------------------------------------
  test_exec <- switch(
    test_family,
    independent_two_sample = .run_independent_two_sample(
      y, groups, grp, normality_violated, variance_violated, test_type, warnings_list
    ),
    paired_two_sample = .run_paired_two_sample(
      y, groups, grp, normality_violated, test_type, warnings_list
    ),
    independent_anova = .run_independent_anova(
      y, grp, normality_violated, variance_violated, test_type, warnings_list
    )
  )
  test_used       <- test_exec$test_used
  test_result     <- test_exec$test_result
  fitted_model    <- test_exec$fitted_model
  parametric_used <- test_exec$parametric
  warnings_list   <- test_exec$warnings_list

  # --- Effect size -----------------------------------------------------------
  es_out <- .compute_effect_size(
    y, grp, groups, test_family, parametric_used,
    variance_violated, fitted_model, warnings_list
  )
  es_result     <- es_out$result
  warnings_list <- es_out$warnings_list

  # --- Post-hoc --------------------------------------------------------------
  posthoc_result <- NULL
  if (perform_posthoc && n_groups > 2 &&
      !is.null(test_result) && test_result$p.value < test_alpha &&
      requireNamespace("rstatix", quietly = TRUE)) {
    posthoc_result <- tryCatch({
      td <- tibble::tibble(y = y, grp = grp)
      ph <- NULL
      if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
        ph <- rstatix::tukey_hsd(td, y ~ grp)
        ph <- .build_posthoc_df(ph, "Tukey HSD", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      } else if (test_used == "Welch's ANOVA") {
        ph <- rstatix::games_howell_test(td, y ~ grp)
        ph <- .build_posthoc_df(ph, "Games-Howell", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      } else if (test_used %in% c("Kruskal-Wallis test", "Kruskal-Wallis test (forced)")) {
        ph <- rstatix::dunn_test(td, y ~ grp, p.adjust.method = p_adjust_method)
        ph <- .build_posthoc_df(ph, "Dunn's test", y, grp, group_ns, groups,
                                parametric_used, test_alpha, p_adjust_method)
      }
      ph
    }, error = function(e) {
      warnings_list <<- c(warnings_list, paste("Post-hoc failed:", e$message))
      NULL
    })
  }

  # --- Data summary ----------------------------------------------------------
  dg <- droplevels(grp)
  data_summary <- data.frame(
    group     = levels(dg),
    n         = as.numeric(table(dg)),
    mean      = as.numeric(tapply(y, dg, mean)),
    sd        = as.numeric(tapply(y, dg, stats::sd)),
    se        = as.numeric(tapply(y, dg,
                  function(v) stats::sd(v) / sqrt(length(v)))),
    median    = as.numeric(tapply(y, dg, stats::median)),
    q25       = as.numeric(tapply(y, dg, function(v) stats::quantile(v, 0.25))),
    q75       = as.numeric(tapply(y, dg, function(v) stats::quantile(v, 0.75))),
    min       = as.numeric(tapply(y, dg, min)),
    max       = as.numeric(tapply(y, dg, max)),
    row.names = NULL
  )

  list(
    test_used          = test_used,
    test_result        = test_result,
    effect_size        = es_result,
    posthoc_result     = posthoc_result,
    norm_results       = norm_results_raw,
    variance_p         = variance_p,
    variance_method    = variance_method,
    normality_violated = normality_violated,
    variance_violated  = variance_violated,
    parametric         = parametric_used,
    data_summary       = data_summary,
    outcome            = outcome,
    group_var          = group,
    groups             = groups,
    n_groups           = n_groups,
    group_ns           = group_ns,
    paired             = paired,
    test_type          = test_type,
    test_alpha         = test_alpha,
    warnings           = warnings_list,
    raw_data           = tibble::tibble(outcome = y, group = grp)
  )
}


# =============================================================================
#  SECTION 9 — .build_summary_table (internal)
# =============================================================================

#' Build a summary_table data frame from a list of run_diff objects
#'
#' Used by the multi-outcome dispatch in \code{run_diff} to construct the
#' \code{$summary_table} element and subgroup-level summary tables.
#'
#' @param results_list Named list of \code{"run_diff"} objects (one per outcome).
#' @param outcome_names Character vector of outcome names (defines row order).
#' @param test_alpha Significance threshold for the \code{significant} column.
#'
#' @return A data frame with one row per outcome.
#' @keywords internal
.build_summary_table <- function(results_list, outcome_names, test_alpha) {

  # Collect the union of all group names across outcomes for average_* columns
  all_groups <- character(0)
  for (oc in outcome_names) {
    res <- results_list[[oc]]
    if (!is.null(res) && !is.null(res$data_summary))
      all_groups <- union(all_groups, as.character(res$data_summary$group))
  }
  avg_col_names <- if (length(all_groups) > 0)
    paste0("average_", all_groups) else character(0)

  tbl_rows <- lapply(outcome_names, function(oc) {
    res <- results_list[[oc]]

    # Skeleton row for failed/NULL outcomes
    if (is.null(res)) {
      base <- data.frame(
        outcome               = oc,
        test_used             = NA_character_,
        statistic             = NA_real_,
        df                    = NA_character_,
        p_value               = NA_real_,
        effect_size_metric    = NA_character_,
        effect_size_estimate  = NA_real_,
        effect_size_ci_low    = NA_real_,
        effect_size_ci_high   = NA_real_,
        effect_size_magnitude = NA_character_,
        significant           = NA,
        n_significant_posthoc = NA_integer_,
        posthoc_pairs         = NA_integer_,
        interpretation        = NA_character_,
        stringsAsFactors      = FALSE
      )
      for (cn in avg_col_names) base[[cn]] <- NA_real_
      return(base)
    }

    tr <- res$test_result
    es <- res$effect_size
    ph <- res$posthoc_result
    ds <- res$data_summary

    # Consistent df formatting via .format_df
    df_val <- .format_df(tr$parameter)

    n_sig_ph <- if (!is.null(ph) && nrow(ph) > 0)
      sum(ph$significant, na.rm = TRUE) else NA_integer_
    total_ph <- if (!is.null(ph) && nrow(ph) > 0) nrow(ph) else NA_integer_

    # Interpretation
    interp_val <- "No significant difference"
    if (!is.null(tr$p.value) && tr$p.value < test_alpha) {
      if (!is.null(ph) && nrow(ph) > 0) {
        sig_interps <- ph$interpretation[
          ph$significant %in% TRUE &
          ph$interpretation != "No significant difference between groups"
        ]
        interp_val <- if (length(sig_interps) > 0)
          paste(sig_interps, collapse = "; ")
        else
          "No significant difference in post hoc"
      } else if (!is.null(ds) && nrow(ds) == 2 && "group" %in% names(ds)) {
        gn <- as.character(ds$group)
        if (isTRUE(res$parametric)) {
          interp_val <- if (ds$mean[1] < ds$mean[2])
            paste0(gn[1], " has lower mean than ",  gn[2])
          else
            paste0(gn[1], " has higher mean than ", gn[2])
        } else {
          y_raw <- res$raw_data$outcome
          g_raw <- res$raw_data$group
          rks   <- rank(y_raw)
          mr    <- tapply(rks, g_raw, mean)
          interp_val <- if (mr[gn[1]] < mr[gn[2]])
            paste0(gn[1], " tends to have smaller values than ", gn[2])
          else
            paste0(gn[1], " tends to have larger values than ", gn[2])
        }
      }
    }

    # Per-group averages
    avg_vals <- stats::setNames(
      vector("list", length(avg_col_names)), avg_col_names
    )
    for (cn in avg_col_names) {
      grp_label    <- sub("^average_", "", cn)
      avg_vals[[cn]] <- if (!is.null(ds) && grp_label %in% ds$group)
        ds$mean[ds$group == grp_label] else NA_real_
    }

    row <- data.frame(
      outcome          = oc,
      test_used        = res$test_used,
      statistic        = if (!is.null(tr$statistic)) as.numeric(tr$statistic[[1]]) else NA_real_,
      df               = df_val,
      p_value          = if (!is.null(tr$p.value)) tr$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
    for (cn in avg_col_names) row[[cn]] <- as.numeric(avg_vals[[cn]])
    row$effect_size_metric    <- if (!is.null(es)) es$metric    else NA_character_
    row$effect_size_estimate  <- if (!is.null(es)) es$estimate  else NA_real_
    row$effect_size_ci_low    <- if (!is.null(es)) es$ci_low    else NA_real_
    row$effect_size_ci_high   <- if (!is.null(es)) es$ci_high   else NA_real_
    row$effect_size_magnitude <- if (!is.null(es)) es$magnitude else NA_character_
    row$significant           <- if (!is.null(tr$p.value)) tr$p.value < test_alpha else NA
    row$n_significant_posthoc <- n_sig_ph
    row$posthoc_pairs         <- total_ph
    row$interpretation        <- interp_val
    row
  })

  tbl           <- do.call(rbind, tbl_rows)
  rownames(tbl) <- NULL
  tbl
}