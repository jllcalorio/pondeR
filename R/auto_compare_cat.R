#' Automatic Statistical Comparison for Categorical Variables
#'
#' This function performs an automatic comparison between two categorical variables.
#' It intelligently selects the appropriate statistical test (Chi-square,
#' Fisher's Exact Test, or McNemar's Test), calculates measures of association,
#' and performs post-hoc tests when appropriate.
#'
#' @details
#' The function follows a decision-making process to determine the most appropriate test:
#'
#' \strong{Test Family Selection:}
#' \itemize{
#'   \item If `paired = TRUE`, the function assumes the data is for matched pairs
#'     and requires a 2x2 table. It performs \strong{McNemar's Test}.
#'   \item If `paired = FALSE`, the function assumes independent samples and
#'     builds a contingency table.
#' }
#'
#' \strong{Test Selection Logic (\code{test_type = "auto"} for independent samples):}
#' The choice between Chi-square and Fisher's test is based on Cochran's rule.
#' \itemize{
#'   \item \strong{Assumption Check:} The function first runs a Chi-square test
#'     to calculate the *expected frequencies* for each cell in the contingency table.
#'   \item \strong{Decision:}
#'     \itemize{
#'       \item If *any* cell has an expected frequency less than `expected_freq_threshold`
#'         (default 5), the Chi-square approximation may be inaccurate. The function
#'         discards this result and performs a \strong{Fisher's Exact Test}.
#'       \item If all cells have an expected frequency >= `expected_freq_threshold`,
#'         the function proceeds with the \strong{Pearson's Chi-square Test}.
#'     }
#' }
#'
#' \strong{Post-Hoc Tests:}
#' For a significant result on an independent table larger than 2x2, a post-hoc
#' test is performed using `rstatix::pairwise_fisher_test` to determine which
#' specific pairs of categories are driving the association. P-values are
#' adjusted using the method specified in `p_adjust_method`.
#'
#' \strong{Measures of Association (Effect Sizes):}
#' The function automatically calculates appropriate measures of association:
#' \itemize{
#'   \item \strong{Phi ($\phi$)} and \strong{Odds Ratio} for 2x2 independent tables.
#'   \item \strong{Cramér's V} for RxC independent tables (larger than 2x2).
#'   \item \strong{Odds Ratio (paired)} for 2x2 paired (McNemar's) tables.
#' }
#'
#' @param data A data frame containing the variables for analysis.
#' @param var1 A character string specifying the name of the first categorical
#'   variable (typically rows).
#' @param var2 A character string specifying the name of the second categorical
#'   variable (typically columns).
#' @param weight An optional character string specifying the name of a weight
#'   column. Use this if your `data` is in a "long" format (e.g., from
#'   `as.data.frame.table()`).
#' @param paired A logical indicating whether the observations are paired
#'   (for McNemar's test). Default is `FALSE`.
#' @param test_type A character string specifying the test strategy for
#'   independent samples. Must be one of `auto` (default), `chisq`, or `fisher`.
#'   This parameter is ignored if `paired = TRUE`.
#' @param expected_freq_threshold The minimum expected cell frequency to
#'   justify using the Chi-square test when `test_type = "auto"`. Default is 5.
#' @param continuity_correction A logical indicating whether to apply Yates'
#'   continuity correction for 2x2 Chi-square tests. Default is `TRUE`.
#' @param fisher_simulate A logical indicating whether to use Monte Carlo
#'   simulation for Fisher's test p-value in tables larger than 2x2
#'   (which can be computationally intensive). Default is `TRUE`.
#' @param test_alpha Significance level for the final statistical test. Default is 0.05.
#' @param calculate_association Logical indicating whether to calculate
#'   measures of association (e.g., Cramér's V, Odds Ratio). Default is `TRUE`.
#' @param perform_posthoc Logical indicating whether to perform post-hoc tests
#'   for significant RxC tables. Default is `TRUE`.
#' @param p_adjust_method P-value adjustment method for post-hoc comparisons.
#'   See `stats::p.adjust.methods`. Default is "BH".
#' @param verbose A logical indicating whether to print detailed messages. Default is TRUE.
#'
#' @return An object of class "auto_compare_cat" which is a list containing:
#'   \item{test_used}{Character string of the statistical test performed.}
#'   \item{test_result}{List containing results of the statistical test (a 'htest' object).}
#'   \item{contingency_table}{The observed contingency table (as a data.frame).}
#'   \item{expected_table}{The expected frequencies table (as a data.frame, NULL if Fisher's/McNemar's).}
#'   \item{residuals}{Pearson residuals (as a data.frame, NULL if Fisher's/McNemar's).}
#'   \item{stdres}{Standardized residuals (as a data.frame, NULL if Fisher's/McNemar's).}
#'   \item{association}{Data frame containing association metrics, magnitude, and interpretation.}
#'   \item{posthoc_result}{Data frame of post-hoc pairwise comparisons (NULL if not applicable).}
#'   \item{assumptions}{Data frame summarizing assumption check results.}
#'   \item{var1}{Name of the row variable.}
#'   \item{var2}{Name of the column variable.}
#'   \item{paired}{The paired setting used.}
#'   \item{test_alpha}{The significance level used.}
#'   \item{warnings}{Character vector of any warnings generated.}
#'   \item{parameters}{List of all parameters used in the analysis.}
#'
#' @importFrom stats chisq.test fisher.test mcnemar.test xtabs na.omit as.formula
#' @importFrom effectsize phi cramers_v oddsratio interpret_phi interpret_cramers_v interpret_oddsratio
#' @importFrom rstatix pairwise_fisher_test
#'
#' @examples
#' # --- Example 1: 2x2 Independent Test (Auto -> Chi-square) ---
#' # Using 'mtcars' to compare transmission (am) vs. engine shape (vs)
#' # This will use Chi-square as all expected counts are > 5.
#' data(mtcars)
#' res_mtcars <- auto_compare_cat(data = mtcars, var1 = "am", var2 = "vs")
#' print(res_mtcars)
#' summary(res_mtcars)
#'
#' # --- Example 2: RxC Independent Test (Auto -> Chi-square) ---
#' # Using 'HairEyeColor' data, which is a table.
#' # We must first convert it to a data frame and use the 'Freq' column as a weight.
#' data(HairEyeColor)
#' he_data <- as.data.frame.table(HairEyeColor)
#' # We'll just compare Hair vs. Eye
#' he_data_sub <- he_data[he_data$Sex == "Male", ] # Subset to avoid 3D
#'
#' res_hair <- auto_compare_cat(
#'   data = he_data_sub,
#'   var1 = "Hair",
#'   var2 = "Eye",
#'   weight = "Freq"
#' )
#' print(res_hair)
#' summary(res_hair) # Will include post-hoc results
#'
#' # --- Example 3: Paired 2x2 Test (McNemar's) ---
#' # We must create a dataset for this.
#' # Imagine a campaign to get people to vote "Yes" (1).
#' set.seed(42)
#' N <- 100
#' before_vote <- sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.7, 0.3))
#' after_vote <- before_vote
#' # 20% of "No" voters switch to "Yes"
#' after_vote[before_vote == "No"][1:14] <- "Yes"
#' # 10% of "Yes" voters switch to "No"
#' after_vote[before_vote == "Yes"][1:3] <- "No"
#'
#' paired_data <- data.frame(
#'   person_id = 1:N,
#'   before = before_vote,
#'   after = after_vote
#' )
#'
#' res_paired <- auto_compare_cat(
#'   data = paired_data,
#'   var1 = "before",
#'   var2 = "after",
#'   paired = TRUE
#' )
#' summary(res_paired)
#'
#' # --- Example 4: Force a test_type ---
#' # Force Chi-square on the mtcars data, with continuity correction
#' res_mtcars_chisq <- auto_compare_cat(
#'   data = mtcars,
#'   var1 = "am",
#'   var2 = "vs",
#'   test_type = "chisq",
#'   continuity_correction = TRUE
#' )
#' print(res_mtcars_chisq)
#'
#' @export
auto_compare_cat <- function(data,
                             var1,
                             var2,
                             weight = NULL,
                             paired = FALSE,
                             test_type = c("auto", "chisq", "fisher"),
                             expected_freq_threshold = 5,
                             continuity_correction = TRUE,
                             fisher_simulate = TRUE,
                             test_alpha = 0.05,
                             calculate_association = TRUE,
                             perform_posthoc = TRUE,
                             p_adjust_method = "BH",
                             verbose = TRUE) {

  # --- 1. Input Validation and Preparation ---
  test_type <- match.arg(test_type)
  warnings_list <- character(0)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  vars_to_check <- if (is.null(weight)) c(var1, var2) else c(var1, var2, weight)
  if (!all(vars_to_check %in% names(data))) {
    stop("One or more specified variables not found in data.")
  }

  # Remove missing values
  data_complete <- stats::na.omit(data[, vars_to_check, drop = FALSE])
  if (nrow(data_complete) < nrow(data)) {
    msg <- sprintf("Removing %d rows with missing values.", nrow(data) - nrow(data_complete))
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }

  # Convert to factors AND drop unused levels (Fixes "Ghost Level" bugs)
  data_complete[[var1]] <- droplevels(as.factor(data_complete[[var1]]))
  data_complete[[var2]] <- droplevels(as.factor(data_complete[[var2]]))

  if (verbose) {
    message("\n=== Auto Compare Categorical Analysis ===")
    message(sprintf("Row: %s | Column: %s | Paired: %s", var1, var2, paired))
    message(sprintf("Test Type: %s", ifelse(paired, "mcnemar", test_type)))
  }

  # --- 2. Initialize Assumption Checks and Build Table ---
  assumptions <- data.frame(
    check = character(),
    result = character(),
    decision = character(),
    stringsAsFactors = FALSE
  )

  tbl <- NULL
  if (is.null(weight)) {
    tbl <- table(data_complete[[var1]], data_complete[[var2]])
  } else {
    if (!is.numeric(data_complete[[weight]])) {
      stop("Weight column must be numeric.")
    }
    tbl_formula <- stats::as.formula(paste(weight, "~", var1, "+", var2))
    tbl <- stats::xtabs(tbl_formula, data = data_complete)
  }

  # Ensure the table is a standard 2D matrix/table without 'xtabs' or other class issues
  tbl <- as.table(as.matrix(tbl))

  # Rename dimensions for clarity in output
  names(dimnames(tbl)) <- c(var1, var2)

  n_rows <- nrow(tbl)
  n_cols <- ncol(tbl)

  if (n_rows < 2 || n_cols < 2) {
    stop("Contingency table must have at least 2 rows and 2 columns.")
  }

  # --- 3. Select and Perform Test ---
  test_used <- NULL
  test_result <- NULL
  expected_table <- NULL
  residuals <- NULL
  stdres <- NULL # Added standardized residuals

  if (paired) {
    # --- 3a. Paired Test (McNemar's) ---
    if (n_rows != 2 || n_cols != 2) {
      stop(sprintf("McNemar's test requires a 2x2 table, but got %dx%d.", n_rows, n_cols))
    }
    test_used <- "McNemar's test"
    test_result <- stats::mcnemar.test(tbl, correct = continuity_correction)
    assumptions <- rbind(assumptions, data.frame(
      check = "Data Structure",
      result = "2x2 table",
      decision = "OK"
    ))

  } else {
    # --- 3b. Independent Samples Tests ---
    run_chisq <- function(warn_small = FALSE) {
      # chisq.test can throw a warning, capture it
      test <- tryCatch(
        {
          suppressWarnings(stats::chisq.test(tbl, correct = continuity_correction))
        },
        warning = function(w) {
          # Store warning
          warnings_list <<- c(warnings_list, w$message)
          # Return the test result anyway
          suppressWarnings(stats::chisq.test(tbl, correct = continuity_correction))
        }
      )

      small_counts <- any(test$expected < expected_freq_threshold)
      if (warn_small && small_counts) {
        msg <- sprintf("Chi-square used, but %d cell(s) have expected counts < %.1f. Consider Fisher's test.",
                       sum(test$expected < expected_freq_threshold), expected_freq_threshold)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(paste("Warning:", msg))
      }
      assumptions <<- rbind(assumptions, data.frame(
        check = sprintf("Expected Freq. > %.1f", expected_freq_threshold),
        result = sprintf("%d/%d cells OK", sum(test$expected >= expected_freq_threshold), length(test$expected)),
        decision = ifelse(small_counts, "VIOLATED", "OK")
      ))

      expected_table <<- test$expected
      residuals <<- test$residuals
      stdres <<- test$stdres
      return(test)
    }

    run_fisher <- function() {
      sim_pval <- (n_rows > 2 || n_cols > 2) && fisher_simulate
      if (sim_pval && verbose) {
        message("Note: Using simulated p-value for Fisher's test on RxC table.")
      }
      test <- stats::fisher.test(tbl, simulate.p.value = sim_pval)
      assumptions <<- rbind(assumptions, data.frame(
        check = "Test Selection",
        result = "Fisher's test used",
        decision = "OK"
      ))
      return(test)
    }

    run_fisher <- function() {
      # Use simulation only if requested AND table is larger than 2x2 (or if user forced it)
      # Standard fisher.test does exact for RxC if workspace allows, otherwise needs sim.
      do_sim <- fisher_simulate

      if (do_sim && verbose) {
        message("Note: Using simulated p-value for Fisher's test (Monte Carlo).")
      }

      test <- tryCatch({
        stats::fisher.test(tbl, simulate.p.value = do_sim)
      }, error = function(e) {
        # Fallback if Exact fails due to memory
        if(!do_sim) {
          msg <- "Exact Fisher test failed (out of memory?). Switching to simulation."
          warnings_list <<- c(warnings_list, msg)
          if(verbose) message(msg)
          return(stats::fisher.test(tbl, simulate.p.value = TRUE))
        } else {
          stop(e)
        }
      })

      assumptions <<- rbind(assumptions, data.frame(
        check = "Test Selection",
        result = "Fisher's test used",
        decision = "OK"
      ))
      return(test)
    }

    if (test_type == "auto") {
      # First, run chisq to check assumptions
      chisq_check <- suppressWarnings(stats::chisq.test(tbl, correct = FALSE)) # No correction for check
      small_counts <- any(chisq_check$expected < expected_freq_threshold)

      if (small_counts) {
        if (verbose) message(sprintf("Note: Expected counts < %.1f. Switching to Fisher's Exact Test.", expected_freq_threshold))
        test_used <- "Fisher's Exact Test"
        test_result <- run_fisher()
      } else {
        test_used <- "Pearson's Chi-square test"
        test_result <- run_chisq()
      }
    } else if (test_type == "chisq") {
      test_used <- "Pearson's Chi-square test (forced)"
      test_result <- run_chisq(warn_small = TRUE)
    } else if (test_type == "fisher") {
      test_used <- "Fisher's Exact Test (forced)"
      test_result <- run_fisher()
    }
  }

  # --- 4. Calculate Measures of Association ---
  association_result <- NULL
  if (calculate_association && requireNamespace("effectsize", quietly = TRUE)) {
    association_result <- tryCatch({
      if (paired) {
        # For paired data, calculate odds ratio from discordant cells
        b <- tbl[1, 2]  # before=No, after=Yes
        c <- tbl[2, 1]  # before=Yes, after=No
        or_est <- b / c
        # Simple CI using normal approximation on log scale
        log_or <- log(or_est)
        se_log_or <- sqrt(1/b + 1/c)
        ci_low <- exp(log_or - 1.96 * se_log_or)
        ci_high <- exp(log_or + 1.96 * se_log_or)

        mag <- if (or_est > 1) "positive change" else if (or_est < 1) "negative change" else "no change"

        list(
          metric = "Odds Ratio (paired)",
          estimate = or_est,
          ci_low = ci_low,
          ci_high = ci_high,
          magnitude = mag,
          interpretation = sprintf("Odds Ratio: %.3f (%s) (95%% CI: %.3f, %.3f)",
                                   or_est, mag, ci_low, ci_high)
        )
      } else if (n_rows == 2 && n_cols == 2) {
        phi_val <- effectsize::phi(tbl, ci = 0.95)
        or_val <- effectsize::oddsratio(tbl, ci = 0.95)
        mag <- effectsize::interpret_phi(abs(phi_val$phi))
        if (is.null(mag) || is.na(mag)) mag <- "unknown"
        list(
          metric = "Phi / Odds Ratio",
          estimate = phi_val$phi,
          ci_low = phi_val$CI_low,
          ci_high = phi_val$CI_high,
          magnitude = as.character(mag),
          interpretation = sprintf("Phi: %.3f (%s), OR: %.3f (95%% CI: %.3f, %.3f)",
                                   phi_val$phi, mag,
                                   or_val$Odds_Ratio, or_val$CI_low, or_val$CI_high)
        )
      } else {
        v_val <- effectsize::cramers_v(tbl, ci = 0.95)
        mag <- effectsize::interpret_cramers_v(v_val$Cramers_V, nrow = n_rows, ncol = n_cols)
        if (is.null(mag) || is.na(mag)) mag <- "unknown"
        list(
          metric = "Cramér's V",
          estimate = v_val$Cramers_V,
          ci_low = v_val$CI_low,
          ci_high = v_val$CI_high,
          magnitude = as.character(mag),
          interpretation = sprintf("Cramér's V: %.3f (%s)",
                                   v_val$Cramers_V, mag)
        )
      }
    }, error = function(e) {
      msg <- paste("Could not calculate association:", e$message)
      warnings_list <<- c(warnings_list, msg)
      NULL
    })
  }

  # --- 5. Perform Post-Hoc Tests ---
  posthoc_result <- NULL
  if (perform_posthoc && !paired && (n_rows > 2 || n_cols > 2) &&
      !is.null(test_result) && test_result$p.value < test_alpha) {

    if (requireNamespace("rstatix", quietly = TRUE)) {
      posthoc_result <- tryCatch({
        if (verbose) message(sprintf("Performing post-hoc pairwise Fisher's tests (p-adjust: %s)...", p_adjust_method))
        # Use pairwise_fisher_test as a robust post-hoc for contingency tables
        rstatix::pairwise_fisher_test(tbl, p.adjust.method = p_adjust_method)
      }, error = function(e) {
        msg <- paste("Could not perform post-hoc test:", e$message)
        warnings_list <<- c(warnings_list, msg)
        NULL
      })
    } else {
      msg <- "Package 'rstatix' not available. Skipping post-hoc tests."
      warnings_list <- c(warnings_list, msg)
    }
  }

  # --- 6. Print Results ---
  if (verbose) {
    message(sprintf("\n--- Test Results ---"))
    message(sprintf("Test Selected: %s", test_used))
    if ("statistic" %in% names(test_result)) {
      stat_name <- names(test_result$statistic)
      stat_val <- as.numeric(test_result$statistic)
      message(sprintf("%s: %.4f", ifelse(is.null(stat_name) || stat_name == "", "Statistic", stat_name), stat_val))
    }

    if (test_result$p.value < 0.001) {
      message(sprintf("P-value: < 0.001"))
    } else {
      message(sprintf("P-value: %.4f", test_result$p.value))
    }

    if (test_result$p.value < test_alpha) {
      message(sprintf("Interpretation: Significant association detected (p < %.2f) ***", test_alpha))
      if (!is.null(association_result)) {
        message(sprintf("%s", association_result$interpretation))
      }
    } else {
      message(sprintf("Interpretation: No significant association detected (p >= %.2f)", test_alpha))
    }
  }

  # --- 7. Return Object ---
  result <- list(
    test_used = test_used,
    test_result = test_result,
    contingency_table = as.data.frame.matrix(tbl),
    expected_table = if (!is.null(expected_table)) as.data.frame.matrix(expected_table) else NULL,
    residuals = if (!is.null(residuals)) as.data.frame.matrix(residuals) else NULL,
    stdres = if (!is.null(stdres)) as.data.frame.matrix(stdres) else NULL,
    association = if (!is.null(association_result)) as.data.frame(association_result) else NULL,
    posthoc_result = posthoc_result,
    assumptions = assumptions,
    var1 = var1,
    var2 = var2,
    paired = paired,
    test_alpha = test_alpha,
    warnings = unique(warnings_list), # Use unique to avoid duplicate warnings
    parameters = list(
      test_type = test_type,
      expected_freq_threshold = expected_freq_threshold,
      continuity_correction = continuity_correction,
      fisher_simulate = fisher_simulate,
      test_alpha = test_alpha,
      calculate_association = calculate_association,
      perform_posthoc = perform_posthoc,
      p_adjust_method = p_adjust_method,
      paired = paired,
      weight = weight
    )
  )

  class(result) <- "auto_compare_cat"
  return(result)
}

#' Print method for auto_compare_cat objects
#' @param x An object of class "auto_compare_cat".
#' @param ... Further arguments passed to or from other methods.
#' @export
print.auto_compare_cat <- function(x, ...) {
  cat("\n=== Automatic Categorical Comparison ===\n\n")
  cat(sprintf("Test Used: %s\n", x$test_used))
  cat(sprintf("Variables: %s vs. %s\n", x$var1, x$var2))
  cat(sprintf("Paired: %s\n\n", x$paired))

  cat("Test Results:\n")
  if ("statistic" %in% names(x$test_result)) {
    stat_name <- names(x$test_result$statistic)
    stat_val <- as.numeric(x$test_result$statistic)
    cat(sprintf("  %s: %.4f\n", ifelse(is.null(stat_name) || stat_name == "", "Statistic", stat_name), stat_val))
  }

  if (x$test_result$p.value < 0.001) {
    cat("  P-value: < 0.001 ***\n")
  } else {
    cat(sprintf("  P-value: %.4f", x$test_result$p.value))
    if (x$test_result$p.value < 0.001) cat(" ***")
    else if (x$test_result$p.value < 0.01) cat(" **")
    else if (x$test_result$p.value < 0.05) cat(" *")
    cat("\n")
  }

  if (x$test_result$p.value < x$test_alpha) {
    cat("  Interpretation: Significant association detected.\n")
  } else {
    cat("  Interpretation: No significant association detected.\n")
  }

  if (!is.null(x$association)) {
    cat(sprintf("  %s\n", x$association$interpretation))
  }

  cat("\n")

  cat("Observed Frequencies:\n")
  print(x$contingency_table)

  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\n\nPost-Hoc Tests (p-adjust: %s):\n", x$parameters$p_adjust_method))
    sig_only <- x$posthoc_result[x$posthoc_result$p.adj < x$test_alpha, ]
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      print(sig_only[, c("group1", "group2", "p", "p.adj")], row.names = FALSE)
    } else {
      cat("No significant pairwise differences found.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  }
}

#' Summary method for auto_compare_cat objects
#' @param object An object of class "auto_compare_cat".
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.auto_compare_cat <- function(object, ...) {
  cat("\n=== Detailed Categorical Comparison Summary ===\n\n")
  cat(sprintf("Row Variable: %s\n", object$var1))
  cat(sprintf("Col Variable: %s\n", object$var2))
  cat(sprintf("Paired Analysis: %s\n", object$paired))
  cat(sprintf("Test Strategy: %s (alpha = %.2f)\n",
              ifelse(object$paired, "mcnemar", object$parameters$test_type),
              object$test_alpha))
  cat(sprintf("Test Selected: %s\n\n", object$test_used))

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
    cat("  (Values > |1.96| or |2.58| suggest contribution to significance)\n")
  }

  cat("\n--- Assumption Checks ---\n")
  if (nrow(object$assumptions) > 0) {
    print(object$assumptions, row.names = FALSE)
  } else {
    cat("  No specific assumptions checked for this test.\n")
  }

  cat("\n--- Statistical Test Results ---\n")
  cat(sprintf("  Method: %s\n", object$test_result$method))
  if ("statistic" %in% names(object$test_result)) {
    stat_name <- names(object$test_result$statistic)
    stat_val <- as.numeric(object$test_result$statistic)
    cat(sprintf("  %s: %.4f\n", ifelse(is.null(stat_name) || stat_name == "", "Statistic", stat_name), stat_val))
  }
  if ("parameter" %in% names(object$test_result)) {
    param_name <- names(object$test_result$parameter)
    param_val <- as.numeric(object$test_result$parameter)
    if(!is.na(param_val)) { # Fisher's has no df
      cat(sprintf("  %s: %d\n", ifelse(is.null(param_name) || param_name == "", "Parameter", param_name), param_val))
    }
  }

  if (object$test_result$p.value < 0.001) {
    cat("  p.value: < 0.001\n")
  } else {
    cat(sprintf("  p.value: %.4f\n", object$test_result$p.value))
  }

  if ("conf.int" %in% names(object$test_result)) {
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", object$test_result$conf.int[1], object$test_result$conf.int[2]))
  }

  if (!is.null(object$association)) {
    cat("\n--- Measures of Association ---\n")
    cat(sprintf("Metric: %s\n", object$association$metric))
    cat(sprintf("Estimate: %.3f", object$association$estimate))
    if ("ci_low" %in% names(object$association) && !is.na(object$association$ci_low)) {
      cat(sprintf(" [95%% CI: %.3f, %.3f]", object$association$ci_low, object$association$ci_high))
    }
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$association$magnitude))
  }

  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons ---\n"))
    cat(sprintf("P-value adjustment method: %s\n\n", object$parameters$p_adjust_method))
    print(object$posthoc_result, row.names = FALSE)
  }

  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT association detected (p < %.2f).\n", object$test_alpha))
    if (!is.null(object$association)) {
      cat(sprintf("The strength of association is %s (%s = %.3f).\n",
                  tolower(object$association$magnitude),
                  object$association$metric,
                  object$association$estimate))
    }
  } else {
    cat(sprintf("Result: NO significant association detected (p >= %.2f).\n", object$test_alpha))
  }

  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  }

  invisible(object)
}
