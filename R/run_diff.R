# =============================================================================
#  .run_diff_single  (unchanged — omitted for brevity)
# =============================================================================

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

  # ---- num_cores ----
  avail_cores <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1L
  if (identical(num_cores, "max")) {
    num_cores <- max(1L, avail_cores - 2L)
  } else {
    num_cores <- as.integer(num_cores)
  }

  # ---- basic checks ----
  if (!is.data.frame(x))
    stop("'x' must be a data frame.")
  if (!outcome %in% names(x))
    stop(sprintf("Outcome '%s' not found in data.", outcome))
  if (!group %in% names(x))
    stop(sprintf("Group '%s' not found in data.", group))

  y   <- x[[outcome]]
  grp <- as.factor(x[[group]])
  if (!is.numeric(y))
    stop("Outcome variable must be numeric.")

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

  # ---- test family ----
  if (n_groups == 2) {
    test_family <- if (paired) "paired_two_sample" else "independent_two_sample"
  } else {
    test_family <- if (paired) "repeated_measures_anova" else "independent_anova"
  }
  if (test_family == "repeated_measures_anova")
    stop("Repeated measures ANOVA is not yet implemented.")
  if (paired && n_groups == 2 && group_ns[1] != group_ns[2])
    stop("Paired tests require equal group sizes after removing missing values.")

  # ---- normality helper ----
  .check_norm <- function(data_vector, group_label = NULL) {
    n <- length(data_vector)
    if (n < 3)
      return(list(violated = TRUE, p_value = NA_real_,
                  method = "Insufficient data", label = group_label))
    use_m <- if (normality_method == "auto") {
      if (n < 50) "shapiro" else "lilliefors"
    } else normality_method

    if (use_m == "shapiro") {
      res <- stats::shapiro.test(data_vector)
      return(list(violated = res$p.value < alpha_normality,
                  p_value  = res$p.value,
                  method   = "Shapiro-Wilk",
                  label    = group_label))
    } else {
      if (!requireNamespace("nortest", quietly = TRUE)) {
        res <- stats::shapiro.test(data_vector)
        warnings_list <<- c(warnings_list,
                            "nortest unavailable; fell back to Shapiro-Wilk.")
        return(list(violated = res$p.value < alpha_normality,
                    p_value  = res$p.value,
                    method   = "Shapiro-Wilk (fallback)",
                    label    = group_label))
      }
      res <- nortest::lillie.test(data_vector)
      return(list(violated = res$p.value < alpha_normality,
                  p_value  = res$p.value,
                  method   = "Lilliefors (K-S)",
                  label    = group_label))
    }
  }

  # ---- assumption checks ----
  normality_violated <- FALSE
  variance_violated  <- FALSE
  norm_results_raw   <- list()
  variance_p         <- NA_real_
  variance_method    <- NA_character_

  if (test_type == "auto") {
    if (n_groups == 2) {
      run_nc <- function(g) .check_norm(y[grp == g], g)
      norm_results_raw <- if (num_cores > 1 && .Platform$OS.type != "windows") {
        parallel::mclapply(groups, run_nc, mc.cores = num_cores)
      } else if (num_cores > 1) {
        cl <- parallel::makeCluster(num_cores)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        parallel::clusterExport(
          cl,
          varlist = c("y", "grp", ".check_norm", "normality_method", "alpha_normality"),
          envir   = environment()
        )
        parallel::clusterEvalQ(cl, {
          if (requireNamespace("nortest", quietly = TRUE)) library(nortest)
        })
        parallel::parLapply(cl, groups, run_nc)
      } else {
        lapply(groups, run_nc)
      }
    } else {
      fit              <- stats::lm(y ~ grp)
      nr               <- .check_norm(fit$residuals, "residuals")
      norm_results_raw <- list(nr)
    }

    for (nr in norm_results_raw) {
      if (inherits(nr, "try-error")) next
      if (isTRUE(nr$violated)) normality_violated <- TRUE
    }

    if (!paired) {
      variance_p <- tryCatch({
        if (requireNamespace("car", quietly = TRUE)) {
          variance_method <- "Levene's test"
          car::leveneTest(y, grp)$`Pr(>F)`[1]
        } else {
          variance_method <- "Bartlett's test"
          stats::bartlett.test(x = y, g = grp)$p.value
        }
      }, error = function(e) {
        warnings_list <<- c(warnings_list,
                            "Variance test failed; assumed homogeneous.")
        1
      })
      if (!is.na(variance_p) && variance_p < alpha_variance)
        variance_violated <- TRUE
    }
  }

  # ---- select & run test ----
  test_used       <- NULL
  test_result     <- NULL
  parametric_used <- TRUE
  fitted_model    <- NULL

  if (test_type == "auto") {
    if (test_family == "independent_two_sample") {
      if (normality_violated) {
        test_used       <- "Mann-Whitney U test (Wilcoxon rank-sum)"
        test_result     <- stats::wilcox.test(x = y[grp == groups[1]],
                                              y = y[grp == groups[2]])
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's t-test"
        test_result <- stats::t.test(x = y[grp == groups[1]],
                                     y = y[grp == groups[2]], var.equal = FALSE)
      } else {
        test_used   <- "Student's t-test"
        test_result <- stats::t.test(x = y[grp == groups[1]],
                                     y = y[grp == groups[2]], var.equal = TRUE)
      }
    } else if (test_family == "paired_two_sample") {
      g1 <- y[grp == groups[1]]
      g2 <- y[grp == groups[2]]
      if (normality_violated) {
        test_used       <- "Wilcoxon signed-rank test"
        test_result     <- stats::wilcox.test(g1, g2, paired = TRUE)
        parametric_used <- FALSE
      } else {
        test_used   <- "Paired t-test"
        test_result <- stats::t.test(g1, g2, paired = TRUE)
      }
    } else {
      # independent_anova
      if (normality_violated) {
        test_used       <- "Kruskal-Wallis test"
        test_result     <- stats::kruskal.test(x = y, g = grp)
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's ANOVA"
        test_result <- stats::oneway.test(stats::formula(y ~ grp), var.equal = FALSE)
      } else {
        test_used    <- "One-way ANOVA"
        fitted_model <- stats::aov(y ~ grp)
        sa           <- summary(fitted_model)[[1]]
        test_result  <- list(
          statistic = c("F value" = sa$`F value`[1]),
          parameter = c(df1 = sa$Df[1], df2 = sa$Df[2]),
          p.value   = sa$`Pr(>F)`[1],
          method    = test_used
        )
      }
    }
  } else if (test_type == "parametric") {
    if (test_family == "independent_two_sample") {
      test_used   <- "Student's t-test (forced)"
      test_result <- stats::t.test(x = y[grp == groups[1]],
                                   y = y[grp == groups[2]], var.equal = TRUE)
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Paired t-test (forced)"
      test_result <- stats::t.test(y[grp == groups[1]], y[grp == groups[2]],
                                   paired = TRUE)
    } else {
      test_used    <- "One-way ANOVA (forced)"
      fitted_model <- stats::aov(y ~ grp)
      sa           <- summary(fitted_model)[[1]]
      test_result  <- list(
        statistic = c("F value" = sa$`F value`[1]),
        parameter = c(df1 = sa$Df[1], df2 = sa$Df[2]),
        p.value   = sa$`Pr(>F)`[1],
        method    = test_used
      )
    }
  } else {
    # nonparametric
    parametric_used <- FALSE
    if (test_family == "independent_two_sample") {
      test_used   <- "Mann-Whitney U test (forced)"
      test_result <- stats::wilcox.test(x = y[grp == groups[1]],
                                        y = y[grp == groups[2]])
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Wilcoxon signed-rank test (forced)"
      # test_result <- stats::wilcox.test(y[grp == groups[1]], y[grp == groups[2]],
      #                                   paired = TRUE)
      test_result <- withCallingHandlers(
        stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
        warning = function(w) {
          if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
            warnings_list <<- c(warnings_list, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    } else {
      test_used   <- "Kruskal-Wallis test (forced)"
      test_result <- stats::kruskal.test(x = y, g = grp)
    }
  }

  # ---- effect size ----
  es_result <- NULL
  if (calculate_effect_size && !is.null(test_result) &&
      requireNamespace("effectsize", quietly = TRUE)) {
    es_result <- tryCatch({
      if (test_family %in% c("independent_two_sample", "paired_two_sample")) {
        is_paired_fam <- (test_family == "paired_two_sample")
        if (parametric_used) {
          es     <- effectsize::cohens_d(
                      x         = y[grp == groups[1]],
                      y         = y[grp == groups[2]],
                      pooled_sd = !variance_violated,
                      verbose   = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "cohen1988")
          list(estimate  = es$Cohens_d,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = if (is_paired_fam) "Cohen's d (paired)" else "Cohen's d")
        } else {
          es     <- effectsize::rank_biserial(
                      x       = y[grp == groups[1]],
                      y       = y[grp == groups[2]],
                      verbose = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "funder2019")
          list(estimate  = es$r_rank_biserial,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = if (is_paired_fam) "Rank-biserial correlation (paired)"
                           else "Rank-biserial correlation")
        }
      } else {
        # independent_anova
        if (parametric_used) {
          es_model <- if (is.null(fitted_model)) {
            warnings_list <<- c(warnings_list,
              "Eta-squared calculated from standard ANOVA (Welch's ANOVA used for test).")
            stats::aov(y ~ grp)
          } else fitted_model
          es     <- effectsize::eta_squared(es_model, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate  = es$Eta2,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = "Eta-squared")
        } else {
          es     <- effectsize::rank_epsilon_squared(x = y, groups = grp,
                                                     verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate  = es$rank_epsilon_squared,
               ci_low    = if (!is.null(es$CI_low))  es$CI_low  else NA_real_,
               ci_high   = if (!is.null(es$CI_high)) es$CI_high else NA_real_,
               magnitude = as.character(interp$Interpretation),
               metric    = "Rank epsilon-squared")
        }
      }
    }, error = function(e) {
      warnings_list <<- c(warnings_list, paste("Effect size failed:", e$message))
      NULL
    })
  }

  # ---- post-hoc ----
  posthoc_result <- NULL
  if (perform_posthoc && n_groups > 2 &&
      !is.null(test_result) && test_result$p.value < test_alpha &&
      requireNamespace("rstatix", quietly = TRUE)) {
    posthoc_result <- tryCatch({
      td <- tibble::tibble(y = y, grp = grp)
      ph <- NULL

      if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
        ph              <- rstatix::tukey_hsd(td, y ~ grp)
        ph$posthoc_test <- "Tukey HSD"
        ph$statistic    <- NA_real_
        ph$p            <- NA_real_
        nl              <- stats::setNames(as.numeric(group_ns), names(group_ns))
        ph$n1           <- nl[ph$group1]
        ph$n2           <- nl[ph$group2]
      } else if (test_used == "Welch's ANOVA") {
        ph              <- rstatix::games_howell_test(td, y ~ grp)
        ph$posthoc_test <- "Games-Howell"
        ph$p            <- NA_real_
      } else if (test_used %in% c("Kruskal-Wallis test",
                                   "Kruskal-Wallis test (forced)")) {
        ph              <- rstatix::dunn_test(td, y ~ grp,
                                              p.adjust.method = p_adjust_method)
        ph$posthoc_test <- "Dunn's test"
      }

      if (!is.null(ph)) {
        if ("p.adj" %in% names(ph)) {
          ph$p.adj.signif <- ifelse(ph$p.adj < 0.001, "***",
                             ifelse(ph$p.adj < 0.01,  "**",
                             ifelse(ph$p.adj < 0.05,  "*", "ns")))
          ph$significant  <- ph$p.adj < test_alpha
        }

        # Stochastic dominance interpretation via mean ranks
        rks        <- rank(y)
        mean_ranks <- tapply(rks, grp, mean)
        ph$interpretation <- mapply(function(g1, g2, sig) {
          if (is.na(sig) || sig == "ns") return("No significant difference")
          if (mean_ranks[g1] < mean_ranks[g2])
            paste0(g1, " tends to have smaller values than ", g2)
          else
            paste0(g1, " tends to have larger values than ", g2)
        },
        ph$group1, ph$group2,
        if ("p.adj.signif" %in% names(ph)) ph$p.adj.signif else rep("ns", nrow(ph)),
        SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }
      ph
    }, error = function(e) {
      warnings_list <<- c(warnings_list, paste("Post-hoc failed:", e$message))
      NULL
    })
  }

  # ---- data summary ----
  dg <- droplevels(grp)
  data_summary <- data.frame(
    group     = levels(dg),
    n         = as.numeric(table(dg)),
    mean      = as.numeric(tapply(y, dg, mean)),
    sd        = as.numeric(tapply(y, dg, stats::sd)),
    se        = as.numeric(tapply(y, dg,
                  function(v) stats::sd(v) / sqrt(length(v)))),
    median    = as.numeric(tapply(y, dg, stats::median)),
    q25       = as.numeric(tapply(y, dg,
                  function(v) stats::quantile(v, 0.25))),
    q75       = as.numeric(tapply(y, dg,
                  function(v) stats::quantile(v, 0.75))),
    min       = as.numeric(tapply(y, dg, min)),
    max       = as.numeric(tapply(y, dg, max)),
    row.names = NULL
  )

  # ---- return rich list (mirrors old run_diff structure for plot_diff) ----
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
#  Helper: build a summary_table data frame from a named list of run_diff
#          objects (one element per outcome).  Used both for the main table
#          and for each subgroup × level slice.
# =============================================================================
.build_summary_table <- function(results_list, outcome_names, test_alpha) {
  tbl_rows <- lapply(outcome_names, function(oc) {
    res <- results_list[[oc]]

    if (is.null(res)) {
      return(data.frame(
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
        stringsAsFactors      = FALSE
      ))
    }

    tr <- res$test_result
    es <- res$effect_size
    ph <- res$posthoc_result

    df_val <- if (!is.null(tr$parameter) && length(tr$parameter) > 0)
      paste(names(tr$parameter), round(tr$parameter, 4), sep = " = ", collapse = ", ")
    else NA_character_

    n_sig_ph <- if (!is.null(ph) && nrow(ph) > 0) sum(ph$significant, na.rm = TRUE) else NA_integer_
    total_ph <- if (!is.null(ph) && nrow(ph) > 0) nrow(ph) else NA_integer_

    data.frame(
      outcome               = oc,
      test_used             = res$test_used,
      statistic             = if (!is.null(tr$statistic)) as.numeric(tr$statistic[[1]]) else NA_real_,
      df                    = df_val,
      p_value               = if (!is.null(tr$p.value)) tr$p.value else NA_real_,
      effect_size_metric    = if (!is.null(es)) es$metric    else NA_character_,
      effect_size_estimate  = if (!is.null(es)) es$estimate  else NA_real_,
      effect_size_ci_low    = if (!is.null(es)) es$ci_low    else NA_real_,
      effect_size_ci_high   = if (!is.null(es)) es$ci_high   else NA_real_,
      effect_size_magnitude = if (!is.null(es)) es$magnitude else NA_character_,
      significant           = if (!is.null(tr$p.value)) tr$p.value < test_alpha else NA,
      n_significant_posthoc = n_sig_ph,
      posthoc_pairs         = total_ph,
      stringsAsFactors      = FALSE
    )
  })

  tbl <- do.call(rbind, tbl_rows)
  rownames(tbl) <- NULL
  tbl
}


#' Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @title Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @description
#' Performs an automatic comparison of a numeric outcome variable across two or more
#' groups. Intelligently selects the appropriate statistical test (parametric or
#' non-parametric) based on checks for normality (of data or residuals) and homogeneity
#' of variances, calculates effect sizes, and performs post-hoc tests when appropriate.
#'
#' When \code{outcome} is a character vector of length > 1 and \code{summary_table = TRUE},
#' an additional \code{$summary_table} element is returned alongside the per-outcome list,
#' providing a single tidy data frame that consolidates key results across all outcomes —
#' useful for high-throughput screening (e.g., metabolomics, proteomics).
#'
#' When \code{subgroup} is also specified and \code{summary_table = TRUE}, additional
#' named elements are appended to the returned list, one per subgroup variable, each
#' containing a named list of per-level summary tables. For example, with
#' \code{subgroup = "Gender"}, the result will include \code{$summary_table_Gender},
#' and inside it \code{$summary_table_Gender$Male} and \code{$summary_table_Gender$Female}.
#'
#' @details
#' The function follows a decision-making process to determine the most appropriate test:
#'
#' \strong{Test Family Selection:}
#' If \code{paired = FALSE}, it compares independent groups. If there are 2 groups, it
#' selects a two-sample test family. If there are more than 2 groups, it selects
#' an ANOVA-family test.
#'
#' \strong{Assumption Checking (only for \code{test_type = "auto"}):}
#' \itemize{
#'   \item \strong{Normality:} Assessed using the Shapiro-Wilk test for small samples
#'     (n < 50). For larger samples (n >= 50), uses the Lilliefors (Kolmogorov-Smirnov)
#'     test from the \code{nortest} package. For two groups, each group is tested.
#'     For more than two groups, the normality of model residuals is tested.
#'   \item \strong{Homogeneity of Variance:} For independent tests, assessed using
#'     Levene's test from the \code{car} package. If \code{car} is not installed,
#'     falls back to Bartlett's test.
#' }
#'
#' \strong{Test Selection Logic (\code{test_type = "auto"}):}
#' \itemize{
#'   \item \strong{Independent Two-Sample:}
#'     \itemize{
#'       \item If normality is violated -> Mann-Whitney U test (Wilcoxon rank-sum).
#'       \item If normality is OK but variance is violated -> Welch's t-test.
#'       \item If both are OK -> Student's t-test.
#'     }
#'   \item \strong{Paired Two-Sample:}
#'     \itemize{
#'       \item If normality is violated -> Wilcoxon signed-rank test.
#'       \item If normality is OK -> Paired t-test.
#'     }
#'   \item \strong{Independent ANOVA (>2 groups):}
#'     \itemize{
#'       \item If normality is violated -> Kruskal-Wallis test.
#'       \item If normality is OK but variance violated -> Welch's ANOVA.
#'       \item If both OK -> One-way ANOVA.
#'     }
#' }
#'
#' \strong{Post-Hoc Tests:}
#' For significant results with >2 groups, appropriate post-hoc tests are automatically
#' performed:
#' \itemize{
#'   \item Tukey HSD for one-way ANOVA
#'   \item Games-Howell for Welch's ANOVA
#'   \item Dunn's test for Kruskal-Wallis
#' }
#'
#' \strong{Effect Sizes:}
#' The function automatically calculates appropriate effect sizes:
#' \itemize{
#'   \item Cohen's d for t-tests
#'   \item Rank-biserial correlation for Mann-Whitney U
#'   \item Eta-squared for ANOVA
#'   \item Epsilon-squared for Kruskal-Wallis
#' }
#'
#' \strong{Non-Parametric Interpretation (Stochastic Dominance):}
#' For non-parametric tests, group comparisons are interpreted using stochastic
#' dominance rather than median comparisons. The \code{interpretation} column in
#' \code{posthoc_result} (and the 2-group case) reads "X tends to have larger/smaller
#' values than Y", reflecting that the Wilcoxon / Mann-Whitney family tests whether
#' one group's values tend to be systematically higher or lower than another's.
#'
#' @param x A data frame containing the variables for analysis.
#' @param outcome A character string or character vector specifying the name(s) of the
#'   numeric outcome variable(s). When a vector is supplied, the function loops over each
#'   outcome and returns a named list of \code{"run_diff"} objects. If
#'   \code{summary_table = TRUE}, an additional \code{$summary_table} element is appended
#'   to the returned list (see \code{summary_table}).
#' @param group A character string specifying the name of the grouping factor variable.
#' @param paired A logical indicating whether the observations are paired. Default is
#'   \code{FALSE}.
#' @param test_type A character string specifying the test strategy. Must be one of
#'   \code{"auto"} (default), \code{"parametric"}, or \code{"nonparametric"}.
#' @param normality_method Method for assessing normality: \code{"shapiro"}
#'   (Shapiro-Wilk test), \code{"lilliefors"} (Lilliefors (Kolmogorov-Smirnov) test),
#'   or \code{"auto"} (uses Shapiro-Wilk for n < 50, Lilliefors for n >= 50).
#'   Default is \code{"auto"}.
#' @param alpha_normality Significance level for normality tests. Default is 0.05.
#' @param alpha_variance Significance level for variance homogeneity tests. Default is
#'   0.05.
#' @param test_alpha Significance level for the final statistical test. Default is 0.05.
#' @param min_n_threshold Minimum sample size required per group. Default is 3.
#' @param calculate_effect_size Logical indicating whether to calculate effect sizes.
#'   Default is \code{TRUE}.
#' @param perform_posthoc Logical indicating whether to perform post-hoc tests for >2
#'   groups. Default is \code{TRUE}.
#' @param p_adjust_method The p-value adjustment method for post-hoc comparisons:
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}. Default is \code{"BH"}.
#' @param group_order Character vector specifying the order of groups. If \code{NULL},
#'   uses factor levels in alphabetical order.
#' @param subgroup Character vector or \code{NULL}. Column name(s) in \code{x} to use for
#'   subgroup analysis. Each must be a categorical column (factor or character) in
#'   \code{x}. When specified, the full analysis is repeated independently for each level
#'   of each subgroup column. Results are returned under \code{$subgroup_analysis} in
#'   the output, nested as \code{subgroup_var -> level -> run_diff result}.
#'   Default is \code{NULL}.
#' @param summary_table Logical. Only relevant when \code{outcome} is a character vector
#'   of length > 1. If \code{TRUE}, a tidy summary data frame is constructed across all
#'   outcomes and appended to the returned list as \code{$summary_table}. When
#'   \code{subgroup} is also specified, additional per-subgroup summary table lists are
#'   appended as \code{$summary_table_<subgroup_var>}, each containing one data frame per
#'   level (e.g., \code{$summary_table_Gender$Male}). The data frame contains one row per
#'   outcome with the following columns:
#'   \code{outcome}, \code{test_used}, \code{statistic}, \code{df}, \code{p_value},
#'   \code{effect_size_metric}, \code{effect_size_estimate}, \code{effect_size_ci_low},
#'   \code{effect_size_ci_high}, \code{effect_size_magnitude}, \code{significant},
#'   \code{n_significant_posthoc} (number of significant post-hoc pairs; \code{NA} for
#'   2-group comparisons), and \code{posthoc_pairs} (total post-hoc pairs tested;
#'   \code{NA} for 2-group comparisons). Default is \code{FALSE}.
#' @param verbose A logical indicating whether to print detailed messages. Default is
#'   \code{TRUE}.
#' @param num_cores Integer specifying the number of cores to use for parallel
#'   processing, or the character string \code{"max"}. If \code{"max"}, the function
#'   uses the number of available cores minus 2 (to reserve system resources), with a
#'   minimum of 1. Default is 1. Supports both Windows and POSIX (Mac/Linux) systems.
#'
#' @return When \code{outcome} is a single string, returns a list of class
#'   \code{"run_diff"} containing:
#'   \item{test_used}{Character string of the statistical test performed.}
#'   \item{test_result}{List containing results of the statistical test.}
#'   \item{effect_size}{List containing effect size estimate, confidence interval,
#'     magnitude, and interpretation.}
#'   \item{posthoc_result}{Data frame of post-hoc pairwise comparisons (\code{NULL} if
#'     not applicable). The \code{interpretation} column uses stochastic dominance
#'     framing for non-parametric tests ("X tends to have larger/smaller values than Y")
#'     and mean-based framing for parametric tests.}
#'   \item{assumptions}{Data frame summarizing assumption check results.}
#'   \item{data_summary}{Data frame with comprehensive descriptive statistics per group.}
#'   \item{outcome}{Name of the outcome variable.}
#'   \item{group_var}{Name of the group variable.}
#'   \item{paired}{The paired setting used.}
#'   \item{test_type}{The test_type strategy used.}
#'   \item{test_alpha}{The significance level used.}
#'   \item{subgroup_analysis}{Named list of subgroup results, nested as
#'     \code{subgroup_var -> level -> run_diff result}. \code{NULL} if no subgroup
#'     specified.}
#'   \item{warnings}{Character vector of any warnings generated.}
#'   \item{parameters}{List of all parameters used in the analysis.}
#'   \item{raw_data}{Data frame containing the variables used.}
#'
#'   When \code{outcome} is a character vector of length > 1, returns a named list where
#'   each element is a \code{"run_diff"} object for the corresponding outcome variable.
#'   If \code{summary_table = TRUE}, the list also contains:
#'   \itemize{
#'     \item \code{$summary_table} — overall tidy summary across all outcomes.
#'     \item \code{$summary_table_<sg>} (one per subgroup variable \code{sg}) — a named
#'       list of per-level summary tables, e.g.
#'       \code{$summary_table_Gender$Male}, \code{$summary_table_Gender$Female}.
#'   }
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
#' @references Aaron Schlege, \emph{Games-Howell Post-Hoc Test}. \url{https://rpubs.com/aaronsc32/games-howell-test}.
#' @references Albers, C., & Lakens, D. (2018). \emph{When power analyses based on pilot data are biased: Inaccurate effect size estimators and follow-up bias}. Journal of Experimental Social Psychology. URL \url{https://doi.org/10.1016/j.jesp.2017.09.004}
#' @references Algina, J., Keselman, H. J., & Penfield, R. D. (2006). \emph{Confidence intervals for an effect size when variances are not equal}. Journal of Modern Applied Statistical Methods, 5(1), 2. URL \url{https://jmasm.com/index.php/jmasm/article/view/222} FULL PDF \url{https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1256&context=jmasm}
#' @references Allen, R. (2017). \emph{Statistics and Experimental Design for Psychologists: A Model Comparison Approach}. World Scientific Publishing Company.
#' @references Bache S, Wickham H (2022). \emph{_magrittr: A Forward-Pipe Operator for R_}. doi:10.32614/CRAN.package.magrittr <https://doi.org/10.32614/CRAN.package.magrittr>, R package version 2.0.3, <https://CRAN.R-project.org/package=magrittr>.
#' @references Bartlett, M. S. (1937). \emph{Properties of sufficiency and statistical tests}. Proceedings of the Royal Society of London Series A 160, 268-282. doi:10.1098/rspa.1937.0109.
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#' @references Benjamini, Y., and Hochberg, Y. (1995). \emph{Controlling the false discovery rate: a practical and powerful approach to multiple testing}. Journal of the Royal Statistical Society Series B, 57, 289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
#' @references Benjamini, Y., and Yekutieli, D. (2001). \emph{The control of the false discovery rate in multiple testing under dependency}. Annals of Statistics, 29, 1165-1188. doi:10.1214/aos/1013699998.
#' @references Ben-Shachar M, Ludecke D, Makowski D (2020). \emph{effectsize: Estimation of Effect Size Indices and Standardized Parameters}. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815
#' @references Chambers, J. M., Freeny, A and Heiberger, R. M. (1992) \emph{Analysis of variance; designed experiments}. Chapter 5 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' @references Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.
#' @references Cliff, N. (1993). \emph{Dominance statistics: Ordinal analyses to answer ordinal questions}. Psychological bulletin, 114(3), 494.
#' @references Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd Ed.)}. New York: Routledge.
#' @references Cureton, E. E. (1956). \emph{Rank-biserial correlation}. Psychometrika, 21(3), 287-290.
#' @references Dallal, G.E. and Wilkinson, L. (1986): \emph{An analytic approximation to the distribution of Lilliefors' test for normality}. The American Statistician, 40, 294-296.
#' @references David F. Bauer (1972). \emph{Constructing confidence sets using rank statistics}. Journal of the American Statistical Association 67, 687-690. doi:10.1080/01621459.1972.10481279.
#' @references Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021, May 7). \emph{Why Hedges' g*s based on the non-pooled standard deviation should be reported with Welch's t-test}. doi:10.31234/osf.io/tu6mp
#' @references Dunn, O. J. (1964) \emph{Multiple comparisons using rank sums Technometrics}, 6(3):241-252.
#' @references Eric Langford (2006) \emph{Quartiles in Elementary Statistics}, Journal of Statistics Education 14 3; doi:10.1080/10691898.2006.11910589
#' @references Fox, J. (2016) \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.
#' @references Fox J, Weisberg S (2019). \emph{_An R Companion to Applied Regression_}, Third edition. Sage, Thousand Oaks CA. <https://www.john-fox.ca/Companion/>.
#' @references Glass, G. V. (1965). \emph{A ranking variable analogue of biserial correlation: Implications for short-cut item analysis}. Journal of Educational Measurement, 2(1), 91-95.
#' @references Gross J, Ligges U (2015). \emph{_nortest: Tests for Normality_}. doi:10.32614/CRAN.package.nortest <https://doi.org/10.32614/CRAN.package.nortest>, R package version 1.0-4, <https://CRAN.R-project.org/package=nortest>.
#' @references Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. Orlando, FL: Academic Press.
#' @references Hochberg, Y. (1988). \emph{A sharper Bonferroni procedure for multiple tests of significance}. Biometrika, 75, 800-803. doi:10.2307/2336325
#' @references Holm, S. (1979). \emph{A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics}, 6, 65-70. \url{https://www.jstor.org/stable/4615733}.
#' @references Hommel, G. (1988). \emph{A stagewise rejective multiple test procedure based on a modified Bonferroni test}. Biometrika, 75, 383-386. doi:10.2307/2336190.
#' @references Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis: Correcting error and bias in research findings}. Sage.
#' @references Hyndman, R. J. and Fan, Y. (1996) \emph{Sample quantiles in statistical packages}, American Statistician 50, 361-365. doi:10.2307/2684934.
#' @references Kassambara A (2023). \emph{_rstatix: Pipe-Friendly Framework for Basic Statistical Tests_}. doi:10.32614/CRAN.package.rstatix <https://doi.org/10.32614/CRAN.package.rstatix>, R package version 0.7.2, <https://CRAN.R-project.org/package=rstatix>.
#' @references Kelley, T. (1935) \emph{An unbiased correlation ratio measure}. Proceedings of the National Academy of Sciences. 21(9). 554-559.
#' @references Kerby, D. S. (2014). \emph{The simple difference formula: An approach to teaching nonparametric correlation}. Comprehensive Psychology, 3, 11-IT.
#' @references King, B. M., & Minium, E. W. (2008). \emph{Statistical reasoning in the behavioral sciences}. John Wiley & Sons Inc.
#' @references Komsta L, Novomestky F (2022). \emph{_moments: Moments, Cumulants, Skewness, Kurtosis and Related Tests_}. doi:10.32614/CRAN.package.moments <https://doi.org/10.32614/CRAN.package.moments>, R package version 0.14.1, <https://CRAN.R-project.org/package=moments>.
#' @references Makkonen, L. and Pajari, M. (2014) \emph{Defining Sample Quantiles by the True Rank Probability}, Journal of Probability and Statistics; Hindawi Publ.Corp. doi:10.1155/2014/326579
#' @references Myles Hollander and Douglas A. Wolfe (1973). \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 27-33 (one-sample), 68-75 (two-sample). Or second edition (1999).
#' @references Myles Hollander and Douglas A. Wolfe (1973), \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 115-120.
#' @references Muller K, Wickham H (2025). \emph{_tibble: Simple Data Frames_}. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' @references Olejnik, S., & Algina, J. (2003). \emph{Generalized eta and omega squared statistics: measures of effect size for some common research designs}. Psychological methods, 8(4), 434.
#' @references Patrick Royston (1982). \emph{Algorithm AS 181: The W test for Normality}. Applied Statistics, 31, 176-180. doi:10.2307/2347986.
#' @references Patrick Royston (1982). \emph{An extension of Shapiro and Wilk's W test for normality to large samples}. Applied Statistics, 31, 115-124. doi:10.2307/2347973.
#' @references Patrick Royston (1995). \emph{Remark AS R94: A remark on Algorithm AS 181: The W test for normality}. Applied Statistics, 44, 547-551. doi:10.2307/2986146.
#' @references Ruxton, G.D., and Beauchamp, G. (2008) \emph{'Time for some a priori thinking about post hoc testing'}, Behavioral Ecology, 19(3), pp. 690-693. doi: 10.1093/beheco/arn020.
#' @references Sangseok Lee, Dong Kyu Lee. \emph{What is the proper way to apply the multiple comparison test?}. Korean J Anesthesiol. 2018;71(5):353-360.
#' @references Sarkar, S. (1998). \emph{Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture}. Annals of Statistics, 26, 494-504. doi:10.1214/aos/1013699998.
#' @references Sarkar, S., and Chang, C. K. (1997). \emph{The Simes method for multiple hypothesis testing with positively dependent test statistics}. Journal of the American Statistical Association, 92, 1601-1608. doi:10.2307/2965431.
#' @references Shaffer, J. P. (1995). \emph{Multiple hypothesis testing. Annual Review of Psychology}, 46, 561-584. doi:10.1146/annurev.ps.46.020195.003021.
#' @references Steiger, J. H. (2004). \emph{Beyond the F test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis}. Psychological Methods, 9, 164-182.
#' @references Stephens, M.A. (1974): \emph{EDF statistics for goodness of fit and some comparisons}. Journal of the American Statistical Association, 69, 730-737.
#' @references Thode Jr., H.C. (2002): \emph{Testing for Normality}. Marcel Dekker, New York.
#' @references Tomczak, M., & Tomczak, E. (2014). \emph{The need to report effect size estimates revisited}. An overview of some recommended measures of effect size.
#' @references Wickham H, Averick M, Bryan J, Chang W, McGowan LD, Francois R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Muller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). \emph{"Welcome to the tidyverse."} _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
#' @references Wicklin, R. (2017) \emph{Sample quantiles: A comparison of 9 definitions}; SAS Blog. https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.htm
#' @references Wright, S. P. (1992). \emph{Adjusted P-values for simultaneous inference}. Biometrics, 48, 1005-1013. doi:10.2307/2532694.
#' @references R Core Team (2025). \emph{_R: A Language and Environment for Statistical Computing_}. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # --- Basic Usage with 3+ Groups ---
#' data(iris)
#' res_iris <- run_diff(x = iris, outcome = "Sepal.Length", group = "Species")
#' print(res_iris)
#' summary(res_iris)
#'
#' # --- Two-Sample Independent Comparison ---
#' data(mtcars)
#' mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
#' res_mtcars <- run_diff(x = mtcars, outcome = "mpg", group = "am")
#' print(res_mtcars)
#'
#' # --- Paired-Sample Comparison ---
#' data(sleep)
#' res_sleep <- run_diff(x = sleep, outcome = "extra", group = "group", paired = TRUE)
#' summary(res_sleep)
#'
#' # --- Multi-outcome with summary table (e.g., metabolomics screening) ---
#' data(iris)
#' outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' res_multi <- run_diff(
#'   x             = iris,
#'   outcome       = outcomes,
#'   group         = "Species",
#'   summary_table = TRUE,
#'   verbose       = FALSE
#' )
#' # Access per-outcome results as usual
#' print(res_multi$Sepal.Length)
#' # Filter features with significant omnibus p-value
#' sig_features <- res_multi$summary_table[res_multi$summary_table$significant, ]
#' print(sig_features)
#'
#' # --- Multi-outcome + subgroup summary tables ---
#' # Suppose df has columns: metabolites..., "Disease", "Gender"
#' res_sg <- run_diff(
#'   x             = df,
#'   outcome       = metabolite_cols,
#'   group         = "Disease",
#'   subgroup      = "Gender",
#'   summary_table = TRUE,
#'   verbose       = FALSE
#' )
#' # Overall summary table
#' res_sg$summary_table
#' # Per-subgroup summary tables
#' res_sg$summary_table_Gender$Male
#' res_sg$summary_table_Gender$Female
#'
#' # --- Forcing a Non-Parametric Test ---
#' data(iris)
#' res_iris_np <- run_diff(
#'   x         = iris,
#'   outcome   = "Sepal.Length",
#'   group     = "Species",
#'   test_type = "nonparametric"
#' )
#' print(res_iris_np)
#'
#' # --- Forcing a Parametric Test & Custom P-Value Adjustment ---
#' data(ToothGrowth)
#' ToothGrowth$dose <- as.factor(ToothGrowth$dose)
#' res_tooth <- run_diff(
#'   x               = ToothGrowth,
#'   outcome         = "len",
#'   group           = "dose",
#'   test_type       = "parametric",
#'   p_adjust_method = "bonferroni"
#' )
#' print(res_tooth$posthoc_result)
#'
#' # --- Custom Group Ordering ---
#' data(PlantGrowth)
#' res_plant <- run_diff(
#'   x           = PlantGrowth,
#'   outcome     = "weight",
#'   group       = "group",
#'   group_order = c("trt1", "ctrl", "trt2"),
#'   verbose     = FALSE
#' )
#' print(res_plant$data_summary)
#' }
#'
#' @export
run_diff <- function(
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
    subgroup              = NULL,
    summary_table         = FALSE,
    verbose               = TRUE,
    num_cores             = 1
) {

  # ---------------------------------------------------------------------------
  # Multi-outcome dispatch
  # ---------------------------------------------------------------------------
  if (length(outcome) > 1) {

    if (!is.logical(summary_table) || length(summary_table) != 1L)
      stop("'summary_table' must be a single logical value (TRUE or FALSE).")

    results <- lapply(outcome, function(oc) {
      tryCatch(
        run_diff(
          x                     = x,
          outcome               = oc,
          group                 = group,
          paired                = paired,
          test_type             = test_type,
          normality_method      = normality_method,
          alpha_normality       = alpha_normality,
          alpha_variance        = alpha_variance,
          test_alpha            = test_alpha,
          min_n_threshold       = min_n_threshold,
          calculate_effect_size = calculate_effect_size,
          perform_posthoc       = perform_posthoc,
          p_adjust_method       = p_adjust_method,
          group_order           = group_order,
          subgroup              = subgroup,
          summary_table         = FALSE,   # single-outcome call; no recursion
          verbose               = verbose,
          num_cores             = num_cores
        ),
        error = function(e) {
          warning(paste0("run_diff failed for outcome '", oc, "': ", e$message),
                  call. = FALSE)
          NULL
        }
      )
    })
    names(results) <- outcome

    # -------------------------------------------------------------------------
    # summary_table block
    # -------------------------------------------------------------------------
    if (summary_table) {

      # 1. Overall summary table (across all outcomes, no subgroup filter)
      results[["summary_table"]] <- .build_summary_table(results, outcome, test_alpha)

      # 2. Per-subgroup summary tables
      #    For each subgroup variable sg and each level lvl, extract the
      #    per-outcome subgroup result and build a separate summary table.
      #    Result is stored as results[[ "summary_table_<sg>" ]][[ lvl ]].
      if (!is.null(subgroup)) {
        for (sg in subgroup) {

          # Determine the levels that were analysed (same logic used inside
          # the single-outcome subgroup block)
          sg_levels <- if (is.factor(x[[sg]])) {
            levels(droplevels(as.factor(x[[sg]])))
          } else {
            sort(unique(na.omit(as.character(x[[sg]]))))
          }

          sg_tables <- list()   # will become results[[ "summary_table_<sg>" ]]

          for (lvl in sg_levels) {

            # For each outcome, pull the run_diff result that lives at
            #   results[[ outcome ]]$subgroup_analysis[[ sg ]][[ lvl ]]
            # and build a flat list keyed by outcome name so
            # .build_summary_table() can consume it directly.
            level_results <- stats::setNames(
              lapply(outcome, function(oc) {
                res <- results[[oc]]
                if (is.null(res)) return(NULL)
                res$subgroup_analysis[[sg]][[lvl]]   # NULL if missing
              }),
              outcome
            )

            sg_tables[[lvl]] <- .build_summary_table(level_results, outcome, test_alpha)
          }

          # Attach to the top-level results list under the name
          # "summary_table_<sg>", e.g. "summary_table_Gender"
          results[[ paste0("summary_table_", sg) ]] <- sg_tables
        }
      }
    }

    return(results)
  }

  # ===========================================================================
  # Single-outcome path  (unchanged from your original — paste your full
  # existing single-outcome code here as-is, starting from step 1)
  # ===========================================================================

  # --- 1. Input Validation and Preparation ---
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  if (!is.logical(summary_table) || length(summary_table) != 1L)
    stop("'summary_table' must be a single logical value (TRUE or FALSE).")

  # --- Subgroup Validation ---
  if (!is.null(subgroup)) {
    if (!is.character(subgroup))
      stop("'subgroup' must be a character vector of column names.")
    missing_sg <- setdiff(subgroup, names(x))
    if (length(missing_sg) > 0)
      stop(paste("Subgroup column(s) not found in data:", paste(missing_sg, collapse = ", ")))
    non_cat <- subgroup[sapply(subgroup, function(s) is.numeric(x[[s]]) && !is.factor(x[[s]]))]
    if (length(non_cat) > 0)
      stop(paste("Subgroup column(s) must be categorical (factor or character):",
                 paste(non_cat, collapse = ", ")))
  }

  avail_cores                         <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1

  if (identical(num_cores, "max")) {
    num_cores <- max(1L, avail_cores - 2L)
    if (verbose)
      message(sprintf("Parallel processing: Using 'max' setting (%d cores active, 2 reserved).",
                      num_cores))
  } else if (is.numeric(num_cores)) {
    if (num_cores < 1 || num_cores > avail_cores || num_cores %% 1 != 0)
      stop(sprintf("'num_cores' must be an integer between 1 and %d, or the string 'max'.",
                   avail_cores))
    num_cores <- as.integer(num_cores)
    if (num_cores > 1 && verbose)
      message(sprintf("Parallel processing: Using %d cores.", num_cores))
  } else {
    stop("Invalid 'num_cores' argument. Must be an integer or 'max'.")
  }

  if (!is.data.frame(x))
    stop("'x' must be a data frame.")
  if (!all(c(outcome, group) %in% names(x)))
    stop(paste0("Outcome variable '", outcome, "' or group variable '", group,
                "' not found in data."))

  y   <- x[[outcome]]
  grp <- as.factor(x[[group]])

  if (!is.numeric(y))
    stop("Outcome variable must be numeric.")

  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(grp)))
      stop("group_order contains levels not present in the data.")
    grp <- factor(grp, levels = group_order)
  }

  complete_cases <- complete.cases(y, grp)
  if (sum(!complete_cases) > 0) {
    msg           <- sprintf("Removing %d rows with missing values.", sum(!complete_cases))
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }
  y   <- y[complete_cases]
  grp <- droplevels(grp[complete_cases])

  if (nlevels(grp) < 2)
    stop("Grouping factor must have at least 2 levels after removing missing values.")

  groups   <- levels(grp)
  n_groups <- length(groups)
  group_ns <- table(grp)
  min_n    <- min(group_ns)

  if (verbose) {
    message("\n=== Auto Compare Analysis ===")
    message(sprintf("Outcome: %s | Group: %s | Paired: %s", outcome, group, paired))
    message(sprintf("Number of groups: %d", n_groups))
    message(sprintf("Sample sizes: %s", paste(names(group_ns), "=", group_ns, collapse = ", ")))
    message(sprintf("Test Type: %s", test_type))
  }

  # --- 2. Initialize Assumption Checks and Determine Test Family ---
  assumptions <- data.frame(
    check            = character(),
    result           = character(),
    p_value          = numeric(),
    decision         = character(),
    stringsAsFactors = FALSE
  )

  assumptions <- rbind(assumptions, data.frame(
    check            = "Minimum sample size",
    result           = sprintf("Min n = %d", min_n),
    p_value          = NA_real_,
    decision         = ifelse(min_n >= min_n_threshold, "PASS", "WARNING"),
    stringsAsFactors = FALSE
  ))

  if (min_n < min_n_threshold) {
    msg           <- sprintf("Warning: Minimum sample size (%d) is below threshold (%d)",
                             min_n, min_n_threshold)
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }

  if (n_groups == 2) {
    test_family <- ifelse(paired, "paired_two_sample", "independent_two_sample")
  } else if (n_groups > 2) {
    test_family <- ifelse(paired, "repeated_measures_anova", "independent_anova")
  }

  if (paired && n_groups == 2) {
    if (group_ns[1] != group_ns[2])
      stop("For paired tests, groups must have equal sample sizes after removing missing values.")
  }

  # --- 3. Helper Function: Check Normality ---
  check_normality <- function(data_vector, group_label = NULL, method = normality_method) {
    n <- length(data_vector)
    if (n < 3)
      return(list(violated = TRUE, p_value = NA_real_, method = "Insufficient data"))
    use_method <- if (method == "auto") { if (n < 50) "shapiro" else "lilliefors" } else method
    if (use_method == "shapiro") {
      shapiro_result <- shapiro.test(data_vector)
      p_val          <- shapiro_result$p.value
      return(list(violated = p_val < alpha_normality, p_value = p_val,
                  method = "Shapiro-Wilk test",
                  label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
    } else {
      if (!requireNamespace("nortest", quietly = TRUE)) {
        warning("Package 'nortest' not available. Falling back to Shapiro-Wilk test.")
        shapiro_result <- shapiro.test(data_vector)
        p_val          <- shapiro_result$p.value
        return(list(violated = p_val < alpha_normality, p_value = p_val,
                    method = "Shapiro-Wilk test (fallback)",
                    label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
      }
      lillie_result <- nortest::lillie.test(data_vector)
      p_val         <- lillie_result$p.value
      return(list(violated = p_val < alpha_normality, p_value = p_val,
                  method = "Lilliefors (K-S) test",
                  label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
    }
  }

  # --- 4. Assumption Checking ---
  normality_violated <- FALSE
  variance_violated  <- FALSE

  if (test_type == "auto") {
    if (n_groups == 2) {
      run_norm_check <- function(g) check_normality(y[grp == g], g, normality_method)
      norm_results   <- list()
      if (num_cores > 1) {
        if (.Platform$OS.type == "windows") {
          cl <- parallel::makeCluster(num_cores)
          on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
          parallel::clusterExport(cl,
                                  varlist = c("y", "grp", "check_normality",
                                              "normality_method", "alpha_normality"),
                                  envir   = environment())
          parallel::clusterEvalQ(cl, { if (requireNamespace("nortest", quietly = TRUE)) library(nortest) })
          norm_results <- parallel::parLapply(cl, groups, run_norm_check)
        } else {
          norm_results <- parallel::mclapply(groups, run_norm_check, mc.cores = num_cores)
        }
      } else {
        norm_results <- lapply(groups, run_norm_check)
      }
      for (norm_check in norm_results) {
        if (inherits(norm_check, "try-error")) { if (verbose) warning("A parallel worker failed."); next }
        assumptions <- rbind(assumptions, data.frame(
          check = norm_check$label, result = norm_check$method,
          p_value = norm_check$p_value,
          decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
          stringsAsFactors = FALSE))
        if (norm_check$violated) normality_violated <- TRUE
      }
    } else {
      fit        <- lm(y ~ grp)
      norm_check <- check_normality(fit$residuals, NULL, normality_method)
      norm_check$label <- "Normality of residuals"
      assumptions <- rbind(assumptions, data.frame(
        check = norm_check$label, result = norm_check$method,
        p_value = norm_check$p_value,
        decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
        stringsAsFactors = FALSE))
      if (norm_check$violated) normality_violated <- TRUE
    }
    if (!paired) {
      variance_check_result <- tryCatch({
        if (requireNamespace("car", quietly = TRUE)) {
          car::leveneTest(y, grp)$`Pr(>F)`[1]
        } else {
          if (verbose) message("Package 'car' not found. Using Bartlett's test.")
          bartlett.test(x = y, g = grp)$p.value
        }
      }, error = function(e) {
        msg <- "Error in variance test. Assuming homogeneity."
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        1
      })
      p_val          <- variance_check_result
      decision       <- ifelse(p_val < alpha_variance, "VIOLATED", "OK")
      test_used_name <- if (requireNamespace("car", quietly = TRUE)) "Levene's test" else "Bartlett's test"
      assumptions <- rbind(assumptions, data.frame(
        check = "Homogeneity of variance", result = test_used_name,
        p_value = p_val, decision = decision, stringsAsFactors = FALSE))
      if (decision == "VIOLATED") variance_violated <- TRUE
    }
  }

  # --- 5. Select and Perform Test ---
  test_used       <- NULL
  test_result     <- NULL
  parametric_used <- TRUE
  fitted_model    <- NULL

  if (test_type == "auto") {
    if (test_family == "independent_two_sample") {
      if (normality_violated) {
        test_used       <- "Mann-Whitney U test (Wilcoxon rank-sum)"
        # test_result     <- wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]])
        test_result <- withCallingHandlers(
            stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
            warning = function(w) {
              if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
                warnings_list <<- c(warnings_list, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's t-test"
        test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = FALSE)
      } else {
        test_used   <- "Student's t-test"
        test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = TRUE)
      }
    } else if (test_family == "paired_two_sample") {
      g1 <- y[grp == groups[1]]; g2 <- y[grp == groups[2]]
      if (normality_violated) {
        test_used <- "Wilcoxon signed-rank test"; test_result <- wilcox.test(g1, g2, paired = TRUE)
        parametric_used <- FALSE
      } else {
        test_used <- "Paired t-test"; test_result <- t.test(g1, g2, paired = TRUE)
      }
    } else if (test_family == "independent_anova") {
      if (normality_violated) {
        test_used       <- "Kruskal-Wallis test"
        test_result     <- kruskal.test(x = y, g = grp)
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's ANOVA"
        test_result <- oneway.test(formula(y ~ grp), var.equal = FALSE)
      } else {
        test_used    <- "One-way ANOVA"
        fitted_model <- aov(y ~ grp)
        summary_aov  <- summary(fitted_model)[[1]]
        test_result  <- list(
          statistic = c("F value" = summary_aov$`F value`[1]),
          parameter = c("df1" = summary_aov$Df[1], "df2" = summary_aov$Df[2]),
          p.value   = summary_aov$`Pr(>F)`[1],
          method    = test_used,
          data.name = deparse(substitute(y))
        )
      }
    }
  } else if (test_type == "parametric") {
    if (test_family == "independent_two_sample") {
      test_used   <- "Student's t-test (forced)"
      test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = TRUE)
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Paired t-test (forced)"
      test_result <- t.test(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE)
    } else if (test_family == "independent_anova") {
      test_used    <- "One-way ANOVA (forced)"
      fitted_model <- aov(y ~ grp)
      summary_aov  <- summary(fitted_model)[[1]]
      test_result  <- list(
        statistic = c("F value" = summary_aov$`F value`[1]),
        parameter = c("df1" = summary_aov$Df[1], "df2" = summary_aov$Df[2]),
        p.value   = summary_aov$`Pr(>F)`[1],
        method    = test_used,
        data.name = deparse(substitute(y))
      )
    }
  } else if (test_type == "nonparametric") {
    parametric_used <- FALSE
    if (test_family == "independent_two_sample") {
      test_used   <- "Mann-Whitney U test (forced)"
      # test_result <- wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]])
      test_result <- withCallingHandlers(
        stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
        warning = function(w) {
          if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
            warnings_list <<- c(warnings_list, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Wilcoxon signed-rank test (forced)"
      test_result <- wilcox.test(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE)
    } else if (test_family == "independent_anova") {
      test_used   <- "Kruskal-Wallis test (forced)"
      test_result <- kruskal.test(x = y, g = grp)
    }
  }

  if (test_family == "repeated_measures_anova")
    stop("Repeated measures ANOVA is not yet implemented.")

  # --- 6. Calculate Effect Size ---
  effect_size_result <- NULL

  if (calculate_effect_size && !is.null(test_result) &&
      requireNamespace("effectsize", quietly = TRUE)) {
    effect_size_result <- tryCatch({
      if (test_family %in% c("independent_two_sample", "paired_two_sample")) {
        is_paired_fam <- (test_family == "paired_two_sample")
        if (parametric_used) {
          es     <- effectsize::cohens_d(x = y[grp == groups[1]], y = y[grp == groups[2]],
                                         pooled_sd = !variance_violated, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "cohen1988")
          list(estimate = es$Cohens_d, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric = if (is_paired_fam) "Cohen's d (paired)" else "Cohen's d",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation), es$Cohens_d))
        } else {
          es     <- effectsize::rank_biserial(x = y[grp == groups[1]], y = y[grp == groups[2]],
                                              verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "funder2019")
          list(estimate = es$r_rank_biserial, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric = if (is_paired_fam) "Rank-biserial correlation (paired)"
                        else "Rank-biserial correlation",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation), es$r_rank_biserial))
        }
      } else if (test_family == "independent_anova") {
        if (parametric_used) {
          es_model <- if (is.null(fitted_model)) {
            warnings_list <<- c(warnings_list,
              "Eta-squared calculated from standard ANOVA (Welch's ANOVA used for test).")
            aov(y ~ grp)
          } else fitted_model
          es     <- effectsize::eta_squared(es_model, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate = es$Eta2, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation), metric = "Eta-squared",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation), es$Eta2))
        } else {
          es     <- effectsize::rank_epsilon_squared(x = y, groups = grp, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate = es$rank_epsilon_squared,
               ci_low  = if (!is.null(es$CI_low))  es$CI_low  else NA_real_,
               ci_high = if (!is.null(es$CI_high)) es$CI_high else NA_real_,
               magnitude = as.character(interp$Interpretation), metric = "Rank epsilon-squared",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation),
                                        es$rank_epsilon_squared))
        }
      }
    }, error = function(e) {
      msg <- paste("Could not calculate effect size:", e$message)
      warnings_list <<- c(warnings_list, msg)
      if (verbose) message(msg)
      NULL
    })
  }

  # --- 7. Perform Post-Hoc Tests ---
  posthoc_result <- NULL

  if (perform_posthoc && n_groups > 2 &&
      !is.null(test_result) && test_result$p.value < test_alpha &&
      requireNamespace("rstatix", quietly = TRUE)) {
    posthoc_result <- tryCatch({
      test_data <- tibble::tibble(y = y, grp = grp)
      ph        <- NULL
      if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
        ph              <- rstatix::tukey_hsd(test_data, y ~ grp)
        ph$posthoc_test <- "Tukey HSD"
        ph$statistic    <- NA_real_
        ph$p            <- NA_real_
        n_lookup        <- as.numeric(group_ns); names(n_lookup) <- names(group_ns)
        ph$n1           <- n_lookup[ph$group1]; ph$n2 <- n_lookup[ph$group2]
      } else if (test_used == "Welch's ANOVA") {
        ph              <- rstatix::games_howell_test(test_data, y ~ grp)
        ph$posthoc_test <- "Games-Howell"
        ph$p            <- NA_real_
      } else if (test_used %in% c("Kruskal-Wallis test", "Kruskal-Wallis test (forced)")) {
        ph              <- rstatix::dunn_test(test_data, y ~ grp, p.adjust.method = p_adjust_method)
        ph$posthoc_test <- "Dunn's test"
      }
      if ("p" %in% names(ph)) {
        ph$p.format <- ifelse(ph$p < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p))
      } else { ph$p.format <- NA_character_ }
      if ("p.adj" %in% names(ph)) {
        ph$p.adj.signif <- ifelse(ph$p.adj < 0.001, "***",
                           ifelse(ph$p.adj < 0.01,  "**",
                           ifelse(ph$p.adj < 0.05,  "*", "ns")))
        ph$p.adj.format <- ifelse(ph$p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p.adj))
        ph$significant  <- ph$p.adj < test_alpha
      } else { ph$p.adj.signif <- NA_character_; ph$p.adj.format <- NA_character_; ph$significant <- NA }
      summary_stats <- data.frame(group = groups, mean = as.numeric(tapply(y, grp, mean)),
                                   stringsAsFactors = FALSE)
      ph <- merge(ph, summary_stats, by.x = "group1", by.y = "group", all.x = TRUE)
      names(ph)[names(ph) == "mean"] <- "mean1"
      ph <- merge(ph, summary_stats, by.x = "group2", by.y = "group", all.x = TRUE)
      names(ph)[names(ph) == "mean"] <- "mean2"
      if ("p.adj.signif" %in% names(ph)) {
        if (parametric_used) {
          ph$interpretation <- ifelse(ph$p.adj.signif == "ns",
            "No significant difference between groups",
            ifelse(ph$mean1 < ph$mean2,
                   paste0(ph$group1, " has lower mean than ",  ph$group2),
                   paste0(ph$group1, " has higher mean than ", ph$group2)))
        } else {
          ranks      <- rank(y); mean_rank1 <- tapply(ranks, grp, mean)
          ph$interpretation <- mapply(function(g1, g2, signif) {
            if (signif == "ns") return("No significant difference between groups")
            if (mean_rank1[g1] < mean_rank1[g2]) paste0(g1, " tends to have smaller values than ", g2)
            else paste0(g1, " tends to have larger values than ", g2)
          }, ph$group1, ph$group2, ph$p.adj.signif, SIMPLIFY = TRUE, USE.NAMES = FALSE)
        }
      } else { ph$interpretation <- NA_character_ }
      final_cols <- c("group1", "group2", "n1", "n2", "statistic", "p", "p.adj",
                      "posthoc_test", "p.format", "p.adj.format", "p.adj.signif",
                      "significant", "mean1", "mean2", "interpretation")
      for (col in final_cols) if (!col %in% names(ph)) ph[[col]] <- NA
      ph[, final_cols]
    }, error = function(e) {
      msg <- paste("Could not perform post-hoc test:", e$message)
      warnings_list <<- c(warnings_list, msg)
      if (verbose) message(msg)
      NULL
    })
  }

  # --- 8. Comprehensive Data Summary ---
  present_groups <- levels(droplevels(grp))
  present_ns     <- table(droplevels(grp))
  data_summary <- data.frame(
    group  = present_groups,
    n      = as.numeric(present_ns),
    mean   = as.numeric(tapply(y, droplevels(grp), mean)),
    sd     = as.numeric(tapply(y, droplevels(grp), sd)),
    se     = as.numeric(tapply(y, droplevels(grp), function(v) sd(v) / sqrt(length(v)))),
    median = as.numeric(tapply(y, droplevels(grp), median)),
    q25    = as.numeric(tapply(y, droplevels(grp), function(v) quantile(v, 0.25))),
    q75    = as.numeric(tapply(y, droplevels(grp), function(v) quantile(v, 0.75))),
    min    = as.numeric(tapply(y, droplevels(grp), min)),
    max    = as.numeric(tapply(y, droplevels(grp), max)),
    row.names = NULL
  )
  if (requireNamespace("moments", quietly = TRUE)) {
    data_summary$skewness <- as.numeric(tapply(y, droplevels(grp), moments::skewness))
    data_summary$kurtosis <- as.numeric(tapply(y, droplevels(grp), moments::kurtosis))
  }

  # --- 9. Print Results ---
  if (verbose) {
    message("\n--- Test Results ---")
    message(sprintf("Test Selected: %s", test_used))
    message(sprintf("Test statistic: %.4f", test_result$statistic))
    if (test_result$p.value < 0.001) message("P-value: < 0.001") else
      message(sprintf("P-value: %.4f", test_result$p.value))
    if (!is.null(effect_size_result)) message(sprintf("%s", effect_size_result$interpretation))
    if (test_result$p.value < test_alpha)
      message(sprintf("Interpretation: Significant difference detected (p < %.2f) ***", test_alpha))
    else
      message(sprintf("Interpretation: No significant difference detected (p >= %.2f)", test_alpha))
    if (!is.null(posthoc_result)) {
      message(sprintf("\n--- Post-Hoc Comparisons (%s) ---", unique(posthoc_result$posthoc_test)))
      sig_comparisons   <- sum(posthoc_result$significant, na.rm = TRUE)
      total_comparisons <- nrow(posthoc_result)
      message(sprintf("Significant pairwise differences: %d out of %d",
                      sig_comparisons, total_comparisons))
    }
  }

  # --- 9b. Subgroup Analysis ---
  subgroup_analysis <- NULL

  if (!is.null(subgroup)) {
    subgroup_analysis <- list()
    for (sg_var in subgroup) {
      sg_levels <- if (is.factor(x[[sg_var]])) {
        levels(droplevels(as.factor(x[[sg_var]])))
      } else {
        sort(unique(na.omit(as.character(x[[sg_var]]))))
      }
      sg_var_results <- list()
      for (lvl in sg_levels) {
        x_sub   <- x[as.character(x[[sg_var]]) == lvl, , drop = FALSE]
        if (nrow(x_sub) == 0) {
          warning(paste0("Subgroup '", sg_var, "' level '", lvl, "' has no observations. Skipping."))
          next
        }
        grp_sub <- droplevels(as.factor(x_sub[[group]]))
        if (nlevels(grp_sub) < 2) {
          warning(paste0("Subgroup '", sg_var, "' = '", lvl,
                         "': fewer than 2 group levels after subsetting. Skipping."))
          next
        }
        sg_result <- tryCatch(
          run_diff(
            x = x_sub, outcome = outcome, group = group, paired = paired,
            test_type = test_type, normality_method = normality_method,
            alpha_normality = alpha_normality, alpha_variance = alpha_variance,
            test_alpha = test_alpha, min_n_threshold = min_n_threshold,
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

  # --- 10. Return Object ---
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
    parameters        = list(
      test_type             = test_type,
      normality_method      = normality_method,
      alpha_normality       = alpha_normality,
      alpha_variance        = alpha_variance,
      test_alpha            = test_alpha,
      min_n_threshold       = min_n_threshold,
      calculate_effect_size = calculate_effect_size,
      perform_posthoc       = perform_posthoc,
      p_adjust_method       = p_adjust_method,
      group_order           = group_order,
      subgroup              = subgroup,
      paired                = paired
    ),
    raw_data = tibble::tibble(outcome = y, group = grp)
  )

  return(result)
}


# =============================================================================
# S3 methods  (unchanged)
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
  cat(sprintf("  Statistic: %.4f\n", x$test_result$statistic))
  if (x$test_result$p.value < 0.001) {
    cat("  P-value: < 0.001 ***\n")
  } else {
    cat(sprintf("  P-value: %.4f", x$test_result$p.value))
    if      (x$test_result$p.value < 0.01) cat(" **")
    else if (x$test_result$p.value < 0.05) cat(" *")
    cat("\n")
  }
  if (!is.null(x$effect_size))
    cat(sprintf("  %s\n", x$effect_size$interpretation))
  cat("\n")
  cat("Data Summary:\n")
  print(x$data_summary[, c("group", "n", "mean", "sd", "median")], row.names = FALSE)
  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\n\nPost-Hoc Tests (%s):\n", unique(x$posthoc_result$posthoc_test)))
    sig_only <- x$posthoc_result[x$posthoc_result$significant %in% TRUE,
                                  c("group1", "group2", "p.adj", "p.adj.signif")]
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      print(sig_only, row.names = FALSE)
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
  cat(sprintf("Test Strategy:     %s (alpha = %.2f)\n", object$test_type, object$test_alpha))
  cat(sprintf("Test Selected:     %s\n\n", object$test_used))
  cat("--- Descriptive Statistics ---\n")
  print(object$data_summary, row.names = FALSE)
  cat("\n--- Assumption Checks ---\n")
  print(object$assumptions, row.names = FALSE)
  cat("\n--- Statistical Test Results ---\n")
  print(object$test_result)
  if (!is.null(object$effect_size)) {
    cat("\n--- Effect Size ---\n")
    cat(sprintf("Metric:    %s\n", object$effect_size$metric))
    cat(sprintf("Estimate:  %.3f", object$effect_size$estimate))
    if (!is.na(object$effect_size$ci_low))
      cat(sprintf(" [95%% CI: %.3f, %.3f]", object$effect_size$ci_low, object$effect_size$ci_high))
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$effect_size$magnitude))
  }
  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons (%s) ---\n",
                unique(object$posthoc_result$posthoc_test)))
    cat(sprintf("P-value adjustment method: %s\n\n", object$parameters$p_adjust_method))
    cols_to_show <- c("group1", "group2", "n1", "n2", "statistic", "p",
                      "p.adj", "p.adj.signif", "interpretation")
    cols_to_show <- cols_to_show[cols_to_show %in% names(object$posthoc_result)]
    print(object$posthoc_result[, cols_to_show], row.names = FALSE)
  }
  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT difference detected between groups (p < %.2f).\n",
                object$test_alpha))
    if (!is.null(object$effect_size))
      cat(sprintf("The effect size is %s (%s = %.3f).\n",
                  tolower(object$effect_size$magnitude),
                  object$effect_size$metric,
                  object$effect_size$estimate))
  } else {
    cat(sprintf("Result: NO significant difference detected between groups (p >= %.2f).\n",
                object$test_alpha))
  }
  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) cat(sprintf("  - %s\n", w))
  }
  invisible(object)
}