#' Automatic Statistical Comparison for Categorical Variables
#'
#' @title Automatic Statistical Comparison for Categorical Variables
#'
#' @description
#' Performs an automatic comparison between two categorical variables. It
#' intelligently selects the appropriate statistical test (Chi-square, Fisher's
#' Exact Test, or McNemar's Test), calculates measures of association, and
#' performs post-hoc tests when appropriate.
#'
#' @details
#' \strong{Test Family Selection:}
#' \itemize{
#'   \item If \code{paired = TRUE}, requires a 2x2 table and performs
#'     \strong{McNemar's Test}.
#'   \item If \code{paired = FALSE}, assumes independent samples and builds a
#'     contingency table.
#' }
#'
#' \strong{Test Selection Logic (\code{test_type = "auto"} for independent samples):}
#' The choice between Chi-square and Fisher's test is based on Cochran's rule.
#' \itemize{
#'   \item The function first calculates expected frequencies via a Chi-square
#'     check.
#'   \item If \emph{any} cell has an expected frequency below
#'     \code{expected_freq_threshold} (default 5), Fisher's Exact Test is used.
#'   \item Otherwise, Pearson's Chi-square Test is used.
#' }
#'
#' \strong{Post-Hoc Tests:}
#' For a significant result on a table where one dimension is binary (2 x C or
#' R x 2), pairwise Fisher's tests are performed via
#' \code{rstatix::pairwise_fisher_test}. P-values are adjusted using
#' \code{p_adjust_method}.
#'
#' \strong{Measures of Association (Effect Sizes):}
#' \itemize{
#'   \item \strong{Phi} and \strong{Odds Ratio} for 2x2 independent tables.
#'   \item \strong{Cramer's V} for RxC independent tables (larger than 2x2).
#'   \item \strong{Odds Ratio (paired)} for 2x2 paired (McNemar's) tables.
#' }
#'
#' @param x A data frame containing the variables for analysis.
#' @param var1 A character vector of one or more column names in \code{x}
#'   specifying the categorical variable(s) to test. Each variable in
#'   \code{var1} will be independently paired with \code{var2} for a test of
#'   association. All listed columns must be categorical (factor or character),
#'   unless \code{force_categorical = TRUE}.
#' @param var2 A character string specifying the name of a single categorical
#'   variable (typically the grouping or outcome column). Each variable in
#'   \code{var1} will be compared against this variable.
#' @param force_categorical Logical. If \code{TRUE}, any column in \code{var1}
#'   or \code{var2} that is not already a factor or character will be coerced
#'   to a factor before analysis. This is useful when numeric codes represent
#'   categorical groups (e.g., \code{0}/\code{1}, \code{1}/\code{2}/\code{3}).
#'   A warning is issued for each coerced column. Default is \code{FALSE}.
#' @param weight An optional character string specifying the name of a weight
#'   column. Use this if \code{x} is in a frequency-aggregated format (e.g.,
#'   from \code{as.data.frame.table()}).
#' @param paired A logical indicating whether the observations are paired (for
#'   McNemar's test). Default is \code{FALSE}.
#' @param test_type A character string specifying the test strategy for
#'   independent samples. One of \code{"auto"} (default), \code{"chisq"}, or
#'   \code{"fisher"}. Ignored if \code{paired = TRUE}.
#' @param expected_freq_threshold The minimum expected cell frequency to justify
#'   using the Chi-square test when \code{test_type = "auto"}. Default is 5.
#' @param continuity_correction A logical indicating whether to apply Yates'
#'   continuity correction for 2x2 Chi-square or McNemar's tests. Default is
#'   \code{TRUE}.
#' @param fisher_simulate A logical indicating whether to use Monte Carlo
#'   simulation for Fisher's Exact Test p-values. Recommended for tables larger
#'   than 2x2. Default is \code{FALSE}.
#' @param test_alpha Significance level for the final statistical test. Default
#'   is 0.05.
#' @param calculate_association Logical. Whether to calculate measures of
#'   association (Cramer's V, Phi, Odds Ratio). Requires \code{effectsize}.
#'   Default is \code{TRUE}.
#' @param perform_posthoc Logical. Whether to perform post-hoc pairwise Fisher's
#'   tests for significant tables with one binary dimension. Requires
#'   \code{rstatix}. Default is \code{TRUE}.
#' @param p_adjust_method P-value adjustment method for post-hoc comparisons.
#'   See \code{\link[stats]{p.adjust}}. Default is \code{"BH"}.
#' @param verbose A logical indicating whether to print detailed messages.
#'   Default is \code{TRUE}.
#'
#' @return A list with class \code{"run_assoc"} containing:
#'   \item{test_used}{Character string of the statistical test performed.}
#'   \item{test_result}{An \code{htest} object from the selected test.}
#'   \item{contingency_table}{Observed contingency table as a data frame.}
#'   \item{expected_table}{Expected frequencies as a data frame (\code{NULL}
#'     for Fisher's/McNemar's).}
#'   \item{residuals}{Pearson residuals as a data frame (\code{NULL} for
#'     Fisher's/McNemar's).}
#'   \item{stdres}{Standardized residuals as a data frame (\code{NULL} for
#'     Fisher's/McNemar's). Values beyond \u00b1 1.96 or \u00b1 2.58 indicate cells
#'     contributing to significance.}
#'   \item{association}{Data frame of association metrics (metric, estimate,
#'     CI, magnitude, interpretation), or \code{NULL}.}
#'   \item{posthoc_result}{Data frame of post-hoc pairwise comparisons, or
#'     \code{NULL}.}
#'   \item{assumptions}{Data frame summarising assumption check results.}
#'   \item{var1}{Name of the row variable.}
#'   \item{var2}{Name of the column variable.}
#'   \item{paired}{The paired setting used.}
#'   \item{test_alpha}{The significance level used.}
#'   \item{warnings}{Character vector of any warnings generated.}
#'   \item{parameters}{List of all parameters used in the analysis.}
#'
#' @importFrom stats chisq.test fisher.test mcnemar.test xtabs na.omit as.formula
#' @importFrom effectsize phi cramers_v oddsratio interpret_phi interpret_cramers_v
#' @importFrom rstatix pairwise_fisher_test
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{run_summarytable}}
#' 
#' @examples
#' \dontrun{
#' # --- Example 1: 2x2 Independent Test (Auto -> Chi-square) ---
#' data(mtcars)
#' res <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
#'                  force_categorical = TRUE)
#' print(res)
#' summary(res)
#'
#' # --- Example 2: RxC Independent Test ---
#' data(HairEyeColor)
#' he <- as.data.frame.table(HairEyeColor)
#' he_male <- he[he$Sex == "Male", ]
#' res_hair <- run_assoc(x = he_male, var1 = "Hair", var2 = "Eye", weight = "Freq")
#' summary(res_hair)
#'
#' # --- Example 3: Paired 2x2 (McNemar's) ---
#' set.seed(42)
#' N <- 100
#' before <- sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.7, 0.3))
#' after  <- before
#' after[before == "No"][1:14]  <- "Yes"
#' after[before == "Yes"][1:3]  <- "No"
#' paired_df <- data.frame(before = before, after = after)
#' res_paired <- run_assoc(x = paired_df, var1 = "before", var2 = "after",
#'                          paired = TRUE)
#' summary(res_paired)
#'
#' # --- Example 4: Force Chi-square ---
#' res_chisq <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
#'                         test_type = "chisq", continuity_correction = TRUE,
#'                         force_categorical = TRUE)
#' print(res_chisq)
#' }
#'
#' @export
run_assoc <- function(x,
                      var1,
                      var2,
                      weight                  = NULL,
                      paired                  = FALSE,
                      force_categorical       = FALSE,
                      test_type               = c("auto", "chisq", "fisher"),
                      expected_freq_threshold = 5,
                      continuity_correction   = TRUE,
                      fisher_simulate         = FALSE,
                      test_alpha              = 0.05,
                      calculate_association   = TRUE,
                      perform_posthoc         = TRUE,
                      p_adjust_method         = "BH",
                      verbose                 = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Input Validation and Preparation
  # ---------------------------------------------------------------------------
  test_type     <- match.arg(test_type)
  warnings_list <- character(0)

  if (!is.data.frame(x))
    stop("'x' must be a data frame.")

  # var1 can now be a vector of column names; var2 must be scalar
  if (!is.character(var1) || length(var1) < 1L)
    stop("'var1' must be a non-empty character vector of column name(s).")
  if (!is.character(var2) || length(var2) != 1L)
    stop("'var2' must be a single character string.")

  vars_to_check <- if (is.null(weight)) c(var1, var2) else c(var1, var2, weight)
  missing_vars  <- setdiff(vars_to_check, names(x))
  if (length(missing_vars) > 0)
    stop(paste("Variable(s) not found in 'x':", paste(missing_vars, collapse = ", ")))

  if (!is.logical(force_categorical) || length(force_categorical) != 1L)
    stop("'force_categorical' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(paired) || length(paired) != 1L)
    stop("'paired' must be a single logical value (TRUE or FALSE).")
  if (!is.numeric(expected_freq_threshold) || length(expected_freq_threshold) != 1L ||
      expected_freq_threshold <= 0)
    stop("'expected_freq_threshold' must be a single positive number.")
  if (!is.logical(continuity_correction) || length(continuity_correction) != 1L)
    stop("'continuity_correction' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(fisher_simulate) || length(fisher_simulate) != 1L)
    stop("'fisher_simulate' must be a single logical value (TRUE or FALSE).")
  if (!is.numeric(test_alpha) || length(test_alpha) != 1L ||
      test_alpha <= 0 || test_alpha >= 1)
    stop("'test_alpha' must be a single numeric value between 0 and 1.")
  if (!is.logical(calculate_association) || length(calculate_association) != 1L)
    stop("'calculate_association' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(perform_posthoc) || length(perform_posthoc) != 1L)
    stop("'perform_posthoc' must be a single logical value (TRUE or FALSE).")
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      !p_adjust_method %in% stats::p.adjust.methods)
    stop(paste0("'p_adjust_method' must be one of: ",
                paste(stats::p.adjust.methods, collapse = ", "), "."))

  # -- Categorical type check / force_categorical coercion --------------------
  cat_cols <- c(var1, var2)
  for (col in cat_cols) {
    is_cat <- is.factor(x[[col]]) || is.character(x[[col]])
    if (!is_cat) {
      if (!force_categorical) {
        stop(sprintf(
          paste0("Column '%s' is not a factor or character (found: %s). ",
                 "Coerce it to a factor first, or set `force_categorical = TRUE`."),
          col, class(x[[col]])[1L]
        ))
      } else {
        msg           <- sprintf(
          "Column '%s' (%s) coerced to factor via `force_categorical = TRUE`.",
          col, class(x[[col]])[1L]
        )
        warnings_list <- c(warnings_list, msg)
        if (verbose) message("Warning: ", msg)
        x[[col]]      <- as.factor(x[[col]])
      }
    }
  }

  # -- If var1 has multiple variables, dispatch a loop and return a list ------
  if (length(var1) > 1L) {
    if (verbose)
      message(sprintf(
        "\n=== run_assoc: running %d pairwise associations with '%s' ===",
        length(var1), var2
      ))
    results <- lapply(var1, function(v) {
      if (verbose) message(sprintf("\n--- Testing: %s vs. %s ---", v, var2))
      tryCatch(
        run_assoc(
          x                       = x,
          var1                    = v,
          var2                    = var2,
          weight                  = weight,
          paired                  = paired,
          force_categorical       = force_categorical,
          test_type               = test_type,
          expected_freq_threshold = expected_freq_threshold,
          continuity_correction   = continuity_correction,
          fisher_simulate         = fisher_simulate,
          test_alpha              = test_alpha,
          calculate_association   = calculate_association,
          perform_posthoc         = perform_posthoc,
          p_adjust_method         = p_adjust_method,
          verbose                 = verbose
        ),
        error = function(e) {
          msg <- sprintf("Error in '%s' vs. '%s': %s", v, var2, conditionMessage(e))
          if (verbose) message("  ERROR: ", msg)
          structure(list(error = msg, var1 = v, var2 = var2), class = "run_assoc_error")
        }
      )
    })
    names(results) <- paste0(var1, "_vs_", var2)
    class(results) <- "run_assoc_multi"
    return(invisible(results))
  }

  # ---------------------------------------------------------------------------
  # 2. Build Contingency Table
  # ---------------------------------------------------------------------------
  assumptions <- data.frame(
    check    = character(),
    result   = character(),
    decision = character(),
    stringsAsFactors = FALSE
  )

  # Drop rows with NA in any analysis column
  cols_for_complete <- if (is.null(weight)) c(var1, var2) else c(var1, var2, weight)
  data_complete     <- x[stats::complete.cases(x[, cols_for_complete, drop = FALSE]), ]

  if (nrow(data_complete) == 0L)
    stop("No complete cases found after removing rows with missing values.")
  if (nrow(data_complete) < nrow(x) && verbose)
    message(sprintf("Note: %d row(s) with missing values removed before analysis.",
                    nrow(x) - nrow(data_complete)))

  if (is.null(weight)) {
    tbl <- table(data_complete[[var1]], data_complete[[var2]])
  } else {
    if (!is.numeric(data_complete[[weight]]))
      stop("'weight' column must be numeric.")
    tbl_formula <- stats::as.formula(paste(weight, "~", var1, "+", var2))
    tbl <- stats::xtabs(tbl_formula, data = data_complete)
  }

  tbl              <- as.table(as.matrix(tbl))
  names(dimnames(tbl)) <- c(var1, var2)

  n_rows <- nrow(tbl)
  n_cols <- ncol(tbl)

  if (n_rows < 2 || n_cols < 2)
    stop("Contingency table must have at least 2 rows and 2 columns.")

  # ---------------------------------------------------------------------------
  # 3. Select and Perform Test
  # ---------------------------------------------------------------------------
  test_used      <- NULL
  test_result    <- NULL
  expected_table <- NULL
  residuals_out  <- NULL
  stdres_out     <- NULL

  # Internal helper: run Chi-square
  .run_chisq <- function(warn_small = FALSE) {
    res        <- suppressWarnings(
      stats::chisq.test(tbl, correct = continuity_correction)
    )
    small      <- any(res$expected < expected_freq_threshold)
    n_small    <- sum(res$expected < expected_freq_threshold)
    n_total    <- length(res$expected)

    if (warn_small && small) {
      msg           <- sprintf(
        "Chi-square used, but %d cell(s) have expected counts < %.1f. Consider Fisher's test.",
        n_small, expected_freq_threshold
      )
      warnings_list <<- c(warnings_list, msg)
      if (verbose) message(paste("Warning:", msg))
    }

    assumptions <<- rbind(assumptions, data.frame(
      check    = sprintf("Expected freq. > %.1f", expected_freq_threshold),
      result   = sprintf("%d/%d cells OK", n_total - n_small, n_total),
      decision = if (small) "VIOLATED" else "OK",
      stringsAsFactors = FALSE
    ))

    expected_table <<- res$expected
    residuals_out  <<- res$residuals
    stdres_out     <<- res$stdres
    res
  }

  # Internal helper: run Fisher's Exact Test
  .run_fisher <- function() {
    res <- tryCatch(
      stats::fisher.test(tbl, simulate.p.value = fisher_simulate),
      error = function(e) {
        if (!fisher_simulate) {
          msg           <- "Exact Fisher's test failed (possibly out of memory). Switching to simulation."
          warnings_list <<- c(warnings_list, msg)
          if (verbose) message(msg)
          stats::fisher.test(tbl, simulate.p.value = TRUE)
        } else {
          stop(e)
        }
      }
    )
    if (fisher_simulate && verbose)
      message("Note: Using simulated p-value for Fisher's test (Monte Carlo).")
    assumptions <<- rbind(assumptions, data.frame(
      check    = "Test Selection",
      result   = "Fisher's Exact Test used",
      decision = "OK",
      stringsAsFactors = FALSE
    ))
    res
  }

  if (paired) {
    if (n_rows != 2 || n_cols != 2)
      stop(sprintf("McNemar's test requires a 2x2 table, but got %dx%d.", n_rows, n_cols))
    test_used   <- "McNemar's test"
    test_result <- stats::mcnemar.test(tbl, correct = continuity_correction)
    assumptions <- rbind(assumptions, data.frame(
      check    = "Data Structure",
      result   = "2x2 paired table",
      decision = "OK",
      stringsAsFactors = FALSE
    ))

  } else {

    if (test_type == "auto") {
      chisq_check <- suppressWarnings(stats::chisq.test(tbl, correct = FALSE))
      if (any(chisq_check$expected < expected_freq_threshold)) {
        if (verbose)
          message(sprintf(
            "Note: Expected counts < %.1f detected. Using Fisher's Exact Test.",
            expected_freq_threshold
          ))
        test_used   <- "Fisher's Exact Test"
        test_result <- .run_fisher()
      } else {
        test_used   <- "Pearson's Chi-square test"
        test_result <- .run_chisq()
      }
    } else if (test_type == "chisq") {
      test_used   <- "Pearson's Chi-square test (forced)"
      test_result <- .run_chisq(warn_small = TRUE)
    } else {
      test_used   <- "Fisher's Exact Test (forced)"
      test_result <- .run_fisher()
    }
  }

  # ---------------------------------------------------------------------------
  # 4. Measures of Association
  # ---------------------------------------------------------------------------
  association_result <- NULL

  if (calculate_association) {
    if (!requireNamespace("effectsize", quietly = TRUE)) {
      msg           <- "Package 'effectsize' not available. Skipping association measures."
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    } else {
      association_result <- tryCatch({
        if (paired) {
          # Odds ratio from discordant cells (McNemar)
          b      <- tbl[1L, 2L]
          cc     <- tbl[2L, 1L]
          or_est <- b / cc
          log_or <- log(or_est)
          se_log <- sqrt(1 / b + 1 / cc)
          ci_lo  <- exp(log_or - 1.96 * se_log)
          ci_hi  <- exp(log_or + 1.96 * se_log)
          mag    <- if (or_est > 1) "positive change" else
                    if (or_est < 1) "negative change" else "no change"
          list(
            metric         = "Odds Ratio (paired)",
            estimate       = or_est,
            ci_low         = ci_lo,
            ci_high        = ci_hi,
            magnitude      = mag,
            interpretation = sprintf(
              "Odds Ratio: %.3f (%s) (95%% CI: %.3f, %.3f)", or_est, mag, ci_lo, ci_hi
            )
          )
        } else if (n_rows == 2L && n_cols == 2L) {
          phi_val <- effectsize::phi(tbl, ci = 0.95)
          or_val  <- effectsize::oddsratio(tbl, ci = 0.95)
          mag     <- as.character(effectsize::interpret_phi(abs(phi_val$phi)))
          list(
            metric         = "Phi / Odds Ratio",
            estimate       = phi_val$phi,
            ci_low         = phi_val$CI_low,
            ci_high        = phi_val$CI_high,
            magnitude      = mag,
            interpretation = sprintf(
              "Phi: %.3f (%s), OR: %.3f (95%% CI: %.3f, %.3f)",
              phi_val$phi, mag,
              or_val$Odds_Ratio, or_val$CI_low, or_val$CI_high
            )
          )
        } else {
          v_val <- effectsize::cramers_v(tbl, ci = 0.95)
          v_est <- v_val$Cramers_V
          mag   <- tryCatch(
            as.character(effectsize::interpret_cramers_v(v_est, nrow = n_rows, ncol = n_cols)),
            error = function(e) "unknown"
          )
          list(
            metric         = "Cramer's V",
            estimate       = v_est,
            ci_low         = v_val$CI_low,
            ci_high        = v_val$CI_high,
            magnitude      = mag,
            interpretation = sprintf("Cramer's V: %.3f (%s)", v_est, mag)
          )
        }
      }, error = function(e) {
        msg           <- paste("Could not calculate association:", e$message)
        warnings_list <<- c(warnings_list, msg)
        NULL
      })
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Post-Hoc Tests
  # ---------------------------------------------------------------------------
  posthoc_result     <- NULL
  posthoc_applicable <- (!paired) && ((n_rows > 2L && n_cols == 2L) ||
                                      (n_cols > 2L && n_rows == 2L))

  if (perform_posthoc && posthoc_applicable &&
      !is.null(test_result) && test_result$p.value < test_alpha) {

    if (!requireNamespace("rstatix", quietly = TRUE)) {
      msg           <- "Package 'rstatix' not available. Skipping post-hoc tests."
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    } else {
      posthoc_result <- tryCatch({
        if (verbose)
          message(sprintf(
            "Performing post-hoc pairwise Fisher's tests (p-adjust: %s)...",
            p_adjust_method
          ))
        # pairwise_fisher_test expects the raw table object directly
        rstatix::pairwise_fisher_test(tbl, p.adjust.method = p_adjust_method)
      }, error = function(e) {
        msg           <- paste("Could not perform post-hoc test:", e$message)
        warnings_list <<- c(warnings_list, msg)
        NULL
      })
    }

  } else if (perform_posthoc && !posthoc_applicable && verbose) {
    message("Post-hoc tests skipped: applicable only for tables where one dimension is binary (2 x C or R x 2).")
  }

  # ---------------------------------------------------------------------------
  # 6. Verbose Output
  # ---------------------------------------------------------------------------
  if (verbose) {
    message("\n--- Test Results ---")
    message(sprintf("Test Selected: %s", test_used))
    if ("statistic" %in% names(test_result)) {
      sn  <- names(test_result$statistic)
      sv  <- as.numeric(test_result$statistic)
      message(sprintf("%s: %.4f", if (is.null(sn) || sn == "") "Statistic" else sn, sv))
    }
    if (test_result$p.value < 0.001)
      message("P-value: < 0.001")
    else
      message(sprintf("P-value: %.4f", test_result$p.value))

    if (test_result$p.value < test_alpha) {
      message(sprintf("Interpretation: Significant association detected (p < %.2f) ***", test_alpha))
      if (!is.null(association_result))
        message(association_result$interpretation)
    } else {
      message(sprintf("Interpretation: No significant association detected (p >= %.2f)", test_alpha))
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Return Object
  # ---------------------------------------------------------------------------
  result <- list(
    test_used         = test_used,
    test_result       = test_result,
    contingency_table = as.data.frame.matrix(tbl),
    expected_table    = if (!is.null(expected_table)) as.data.frame.matrix(expected_table) else NULL,
    residuals         = if (!is.null(residuals_out))  as.data.frame.matrix(residuals_out)  else NULL,
    stdres            = if (!is.null(stdres_out))      as.data.frame.matrix(stdres_out)     else NULL,
    association       = if (!is.null(association_result)) {
      as.data.frame(
        lapply(association_result, function(v) if (length(v) == 0L || is.null(v)) NA else v),
        stringsAsFactors = FALSE
      )
    } else NULL,
    posthoc_result    = posthoc_result,
    assumptions       = assumptions,
    var1              = var1,
    var2              = var2,
    paired            = paired,
    test_alpha        = test_alpha,
    warnings          = unique(warnings_list),
    parameters        = list(
      test_type               = test_type,
      expected_freq_threshold = expected_freq_threshold,
      continuity_correction   = continuity_correction,
      fisher_simulate         = fisher_simulate,
      test_alpha              = test_alpha,
      calculate_association   = calculate_association,
      perform_posthoc         = perform_posthoc,
      p_adjust_method         = p_adjust_method,
      paired                  = paired,
      weight                  = weight
    )
  )

  class(result) <- "run_assoc"
  result
}


# =============================================================================
# S3 Methods
# =============================================================================

#' Print method for run_assoc objects
#' @param x An object from \code{run_assoc}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.run_assoc <- function(x, ...) {
  cat("\n=== Automatic Categorical Comparison ===\n\n")
  cat(sprintf("Test Used: %s\n", x$test_used))
  cat(sprintf("Variables: %s vs. %s\n", x$var1, x$var2))
  cat(sprintf("Paired: %s\n\n", x$paired))

  cat("Test Results:\n")
  if ("statistic" %in% names(x$test_result)) {
    sn <- names(x$test_result$statistic)
    sv <- as.numeric(x$test_result$statistic)
    cat(sprintf("  %s: %.4f\n",
                if (is.null(sn) || sn == "") "Statistic" else sn, sv))
  }
  if (x$test_result$p.value < 0.001) {
    cat("  P-value: < 0.001 ***\n")
  } else {
    pv <- x$test_result$p.value
    cat(sprintf("  P-value: %.4f%s\n", pv,
                if (pv < 0.01) " **" else if (pv < 0.05) " *" else ""))
  }
  cat(sprintf("  Interpretation: %s\n",
              if (x$test_result$p.value < x$test_alpha)
                "Significant association detected."
              else
                "No significant association detected."))

  if (!is.null(x$association))
    cat(sprintf("  %s\n", x$association$interpretation))

  cat("\nObserved Frequencies:\n")
  print(x$contingency_table)

  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\nPost-Hoc Tests (p-adjust: %s):\n", x$parameters$p_adjust_method))
    sig_only <- x$posthoc_result[
      !is.na(x$posthoc_result$p.adj) & x$posthoc_result$p.adj < x$test_alpha, ]
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      cols <- intersect(c("group1", "group2", "p", "p.adj"), names(sig_only))
      print(sig_only[, cols], row.names = FALSE)
    } else {
      cat("No significant pairwise differences found.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat(sprintf("  - %s\n", w))
  }
  invisible(x)
}

#' Summary method for run_assoc objects
#' @param object An object from \code{run_assoc}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.run_assoc <- function(object, ...) {
  cat("\n=== Detailed Categorical Comparison Summary ===\n\n")
  cat(sprintf("Row Variable:    %s\n", object$var1))
  cat(sprintf("Col Variable:    %s\n", object$var2))
  cat(sprintf("Paired:          %s\n", object$paired))
  cat(sprintf("Test Strategy:   %s (alpha = %.2f)\n",
              if (object$paired) "mcnemar" else object$parameters$test_type,
              object$test_alpha))
  cat(sprintf("Test Selected:   %s\n\n", object$test_used))

  cat("--- Observed Frequencies ---\n")
  print(object$contingency_table)

  if (!is.null(object$expected_table)) {
    cat("\n--- Expected Frequencies ---\n")
    print(round(object$expected_table, 2))
  }

  if (!is.null(object$residuals)) {
    cat("\n--- Pearson Residuals ---\n")
    print(round(object$residuals, 2))
  }

  if (!is.null(object$stdres)) {
    cat("\n--- Standardized Residuals ---\n")
    print(round(object$stdres, 2))
    cat("  (Values beyond |1.96| or |2.58| indicate cells driving the association)\n")
  }

  cat("\n--- Assumption Checks ---\n")
  if (nrow(object$assumptions) > 0)
    print(object$assumptions, row.names = FALSE)
  else
    cat("  No specific assumptions checked for this test.\n")

  cat("\n--- Statistical Test Results ---\n")
  cat(sprintf("  Method: %s\n", object$test_result$method))
  if ("statistic" %in% names(object$test_result)) {
    sn <- names(object$test_result$statistic)
    sv <- as.numeric(object$test_result$statistic)
    cat(sprintf("  %s: %.4f\n",
                if (is.null(sn) || sn == "") "Statistic" else sn, sv))
  }
  if ("parameter" %in% names(object$test_result)) {
    pn  <- names(object$test_result$parameter)
    pv2 <- as.numeric(object$test_result$parameter)
    if (!is.na(pv2))
      cat(sprintf("  %s: %d\n",
                  if (is.null(pn) || pn == "") "Parameter" else pn,
                  as.integer(pv2)))
  }
  if (object$test_result$p.value < 0.001)
    cat("  p-value: < 0.001\n")
  else
    cat(sprintf("  p-value: %.4f\n", object$test_result$p.value))

  if ("conf.int" %in% names(object$test_result))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
                object$test_result$conf.int[1L], object$test_result$conf.int[2L]))

  if (!is.null(object$association)) {
    cat("\n--- Measures of Association ---\n")
    cat(sprintf("Metric:    %s\n", object$association$metric))
    cat(sprintf("Estimate:  %.3f", object$association$estimate))
    if (!is.na(object$association$ci_low))
      cat(sprintf(" [95%% CI: %.3f, %.3f]",
                  object$association$ci_low, object$association$ci_high))
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$association$magnitude))
  }

  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons ---\n"))
    cat(sprintf("P-value adjustment: %s\n\n", object$parameters$p_adjust_method))
    print(object$posthoc_result, row.names = FALSE)
  }

  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT association detected (p < %.2f).\n", object$test_alpha))
    if (!is.null(object$association))
      cat(sprintf("Strength: %s (%s = %.3f).\n",
                  tolower(object$association$magnitude),
                  object$association$metric,
                  object$association$estimate))
  } else {
    cat(sprintf("Result: NO significant association detected (p >= %.2f).\n", object$test_alpha))
  }

  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) cat(sprintf("  - %s\n", w))
  }

  invisible(object)
}

#' Print method for run_assoc_multi objects
#' @param x An object of class \code{run_assoc_multi} from \code{run_assoc}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.run_assoc_multi <- function(x, ...) {
  cat(sprintf("\n=== run_assoc: %d pairwise association result(s) ===\n", length(x)))
  for (nm in names(x)) {
    cat(sprintf("\n--- %s ---\n", nm))
    if (inherits(x[[nm]], "run_assoc_error")) {
      cat(sprintf("  ERROR: %s\n", x[[nm]]$error))
    } else {
      print(x[[nm]])
    }
  }
  invisible(x)
}

#' Summary method for run_assoc_multi objects
#' @param object An object of class \code{run_assoc_multi} from \code{run_assoc}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.run_assoc_multi <- function(object, ...) {
  cat(sprintf("\n=== run_assoc: Summary of %d pairwise association(s) ===\n", length(object)))
  for (nm in names(object)) {
    cat(sprintf("\n--- %s ---\n", nm))
    if (inherits(object[[nm]], "run_assoc_error")) {
      cat(sprintf("  ERROR: %s\n", object[[nm]]$error))
    } else {
      summary(object[[nm]])
    }
  }
  invisible(object)
}

#' Print method for run_assoc_error objects
#' @param x An object of class \code{run_assoc_error}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.run_assoc_error <- function(x, ...) {
  cat(sprintf("\n[run_assoc ERROR] %s vs. %s\n  %s\n", x$var1, x$var2, x$error))
  invisible(x)
}