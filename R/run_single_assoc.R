#' Run Single-Variable Association Tests
#'
#' This function performs univariate statistical tests on categorical variables 
#' to assess the distribution of levels. It automatically applies \code{\link{binom.test}} 
#' for binomial variables (exactly 2 levels) and \code{\link{chisq.test}} for multinomial 
#' variables (> 2 levels).
#'
#' @param x A data frame containing the variables to be analyzed.
#' @param y A character vector of column names in \code{x} to analyze. 
#'   Defaults to all categorical (factor or character) variables.
#' @param subgroup A character vector of column names in \code{x} to use for 
#'   stratified analysis. If provided, tests are run within each level of each subgroup.
#' @param force_categorical A character vector of column names that should be 
#'   treated as categorical even if they are stored as other types (e.g., integers).
#' @param ... Additional parameters passed forward to both \code{\link{binom.test}} 
#'   and \code{\link{chisq.test}} (e.g., \code{p} for null hypothesis proportions 
#'   or \code{conf.level}).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{summary_table}: A data frame summarizing the test results.
#'     Columns include: \code{subgroup}, \code{group_level}, \code{feature}, 
#'     \code{method}, \code{estimate} (estimated proportion for binomial tests),
#'     \code{lower_95_ci}, \code{upper_95_ci} (confidence interval for binomial tests),
#'     \code{test_statistic} (e.g., Chi-squared value for chi-squared tests),
#'     \code{p_value}, \code{note}, and \code{interpretation}.
#'     For tests where a specific metric is not applicable, its column will contain \code{NA}.
#'   \item \code{raw_results}: A named list containing the full test objects 
#'     returned by the underlying R functions.
#' }
#' 
#' @importFrom stats binom.test
#' 
#' @author John Lennon L. Calorio
#' 
#' @export
#' @examples
#' # Example 1: Binomial test (exactly 2 levels)
#' # Testing if the distribution of 'am' (Transmission) in mtcars is 50/50
#' run_single_assoc(mtcars, y = "am", force_categorical = "am")
#'
#' # Example 2: Chi-squared test (more than 2 levels)
#' # Testing if the distribution of 'Species' in iris is uniform
#' run_single_assoc(iris, y = "Species")
#'
#' # Example 3: Subgroup analysis
#' # Testing 'am' distribution stratified by 'vs' (Engine shape)
#' run_single_assoc(mtcars, y = "am", subgroup = "vs", force_categorical = c("am", "vs"))
run_single_assoc <- function(x, 
                             y = NULL, 
                             subgroup = NULL,
                             force_categorical = NULL, 
                             ...) {
  
  # 1. Input Validation
  if (!is.data.frame(x)) {
    stop("Input 'x' must be a data frame.")
  }

  # Copy to local df to avoid side effects
  df_proc <- x

  # 2. Handle force_categorical
  if (!is.null(force_categorical)) {
    for (col in force_categorical) {
      if (col %in% names(df_proc)) {
        df_proc[[col]] <- as.factor(df_proc[[col]])
      } else {
        warning(paste("Variable", col, "not found in data frame. Skipping conversion."))
      }
    }
  }

  # 3. Identify target variables (y)
  if (is.null(y)) {
    # Detect factors, characters, and logicals
    y <- names(df_proc)[vapply(df_proc, function(v) is.factor(v) || is.character(v) || is.logical(v), logical(1))]
  } else {
    # Filter provided variables to keep only categorical ones (factors, characters, logicals)
    is_cat <- vapply(y, function(nm) {
      val <- df_proc[[nm]]
      is.factor(val) || is.character(val) || is.logical(val)
    }, logical(1))
    if (any(!is_cat)) {
      warning(paste("Skipping non-categorical variables:", paste(y[!is_cat], collapse = ", ")))
      y <- y[is_cat]
    }
  }

  if (length(y) == 0) {
    message("No categorical variables identified for analysis.")
    return(NULL)
  }

  # 4. Iterative Testing
  res_list <- list()
  table_rows <- list()

  # Helper to determine stratification
  strat_list <- if (is.null(subgroup)) {
    list(list(name = "Global", level = "All", data = df_proc))
  } else {
    unlist(lapply(subgroup, function(sub) {
      lvls <- if (is.factor(df_proc[[sub]])) {
        levels(droplevels(df_proc[[sub]]))
      } else {
        sort(unique(stats::na.omit(df_proc[[sub]])))
      }
      lapply(lvls, function(l) {
        list(name = sub, level = as.character(l), data = df_proc[df_proc[[sub]] == l, , drop = FALSE])
      })
    }), recursive = FALSE)
  }

  for (strat in strat_list) {
    current_df <- strat$data
    strat_name <- strat$name
    strat_lvl  <- strat$level

    for (var_name in y) {
      # Skip the variable if it is currently being used for stratification
      if (var_name == strat_name) next

      res_key <- if (strat_name == "Global") var_name else paste(strat_name, strat_lvl, var_name, sep = "_")
      vec <- current_df[[var_name]]
      counts <- table(vec, useNA = "no")
      n_levels <- length(counts)

      if (n_levels < 2) {
        table_rows[[res_key]] <- data.frame(
          subgroup = strat_name,
          group_level = strat_lvl,
          feature = var_name,
          method = "None",
          estimate = NA_real_,
          lower_95_ci = NA_real_,
          upper_95_ci = NA_real_,
          test_statistic = NA_real_,
          p_value = NA_real_,
          note = "Insufficient levels (<2)",
          interpretation = "Not applicable (insufficient levels)",
          stringsAsFactors = FALSE
        )
        next
      }

      if (n_levels == 2) {
        test_res <- binom.test(counts, ...)
        res_list[[res_key]] <- test_res

        lower_95_ci <- test_res$conf.int[1]
        upper_95_ci <- test_res$conf.int[2]

        interpretation_text <- ifelse(test_res$p.value < 0.05, 
                                      "Significantly different", 
                                      "Not significantly different")

        table_rows[[res_key]] <- data.frame(
          subgroup = strat_name,
          group_level = strat_lvl,
          feature = var_name,
          method = test_res$method,
          estimate = test_res$estimate[[1]],
          lower_95_ci = lower_95_ci,
          upper_95_ci = upper_95_ci,
          test_statistic = NA_real_,
          p_value = test_res$p.value,
          note = paste0("Prop: ", names(test_res$estimate)[1]),
          interpretation = interpretation_text,
          stringsAsFactors = FALSE
        )
      } else {
        test_res <- chisq.test(counts, ...)
        res_list[[res_key]] <- test_res

        interpretation_text <- ifelse(test_res$p.value < 0.05, 
                                      "Significantly different distribution", 
                                      "Not significantly different distribution")

        table_rows[[res_key]] <- data.frame(
          subgroup = strat_name,
          group_level = strat_lvl,
          feature = var_name,
          method = test_res$method,
          estimate = NA_real_,
          lower_95_ci = NA_real_,
          upper_95_ci = NA_real_,
          test_statistic = unname(test_res$statistic),
          p_value = test_res$p.value,
          note = "Goodness of fit",
          interpretation = interpretation_text,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary_table <- do.call(rbind, table_rows)
  rownames(summary_table) <- NULL

  return(list(summary_table = summary_table, raw_results = res_list))
}