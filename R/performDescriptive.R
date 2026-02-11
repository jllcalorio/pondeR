#' Performs Descriptive Statistics on a Dataframe
#'
#' @description
#' This function generates a comprehensive descriptive statistics table from a given dataframe,
#' optionally stratified by a categorical variable. It provides customizable summaries for both
#' continuous and categorical variables, including control over rounding, missing data reporting,
#' labeling, formatting, and table aesthetics. The output is a publication-ready `gtsummary` table
#' that can be easily integrated into reports, manuscripts, or exported (via copy-paste for now) for further use.
#'
#' @param df Dataframe. The input dataframe. The first row is assumed to contain column names.
#' @param summarize_what Character vector. A character vector of column names that will be included in the analysis, including what will be placed in `split_by` and `strata_by`.
#' @param split_by String. A column name in `df` (must be a categorical variable with at least 2 unique levels/factors) by which to stratify the descriptive statistics. If `NULL` (default), no stratification is performed.
#' @param split_by_header String. The header name to display for the `split_by` variable in the output table. Only applicable if `split_by` is not `NULL`, and takes on the value of `split_by`. Defaults to `NULL` initially.
#' @param strata_by String. A column name in `df` (must be a categorical variable with at least 2 unique levels/factors) by which to stratify the descriptive statistics 'before' `split_by`. If `NULL` (default), no stratification is performed.
#' @param rename_variables List. A list of the format `list("original_variable_name1" ~ "new_variable_name1", ...)`, used to rename variables in the output table.
#' @param continuous_statistics String. The type of statistics to report for continuous variables. Options are:
#'   \itemize{
#'     \item 'meanSD': mean ± standard deviation
#'     \item 'meanSD2': mean (standard deviation), which is another format of meanSD
#'     \item 'medianIQR': median (interquartile range) or equivalently median (p25, p75)
#'     \item 'mean': mean (average)
#'     \item 'sd': standard deviation
#'     \item 'p0': 0th percentile
#'     \item 'p25': 25th percentile (1st quartile)
#'     \item 'p50': 50th percentile (2nd quartile) (median)
#'     \item 'p75': 75th percentile (3rd quartile)
#'     \item 'p100': 100th percentile
#'     \item 'IQR': p25, p75 (interquartile range)
#'     }
#'     Defaults to 'meanSD'.
#' @param categorical_statistics String. The type of statistics to report for categorical variables.
#'   \itemize{
#'     \item 'n_percent': Displays count and percentage of valid values (e.g., '10 (5.0%)').
#'     \item 'n': Displays only the count of valid values.
#'     \item 'percent': Displays only the percentage of valid values.
#'     }
#'     Defaults to 'n_percent'.
#' Options are 'n_percent' (count (percentage)), 'n' (count only), or 'percent' (percentage only). Defaults to 'n_percent'.
#' @param force_continuous Vector. A character vector of column names that should be treated as continuous, even if their data type suggests otherwise (e.g., numeric variables with few unique values that might be misclassified as categorical).
#' @param force_categorical Vector. A character vector of column names that should be treated as categorical, even if their data type suggests otherwise (e.g., numeric variables that are truly categories, or binary variables that should be presented as full categories rather than dichotomous).
#' @param n_digits_continuous Numeric Vector or String.
#'   \itemize{
#'     \item c(2, 2): A vector of two integers specifying the number of digits to round for mean and standard deviation (or median and IQR) respectively, for continuous variables. For example, `c(2, 2)` means 2 decimal places for both (default).
#'     \item "dynamic": The other option is 'dynamic' which selects the number of decimal places dynamically.
#'   }
#' @param n_digits_categorical Numeric Vector or String.
#'   \itemize{
#'     \item c(0, 2): A vector of two integers specifying the number of digits to round for counts and percentages respectively, for categorical variables. For example, `c(0, 2)` means 0 decimal places for counts and 2 for percentages (default).
#'     \item "dynamic": Same as in `n_digits_continuous`.
#'   }
#' @param display_missing String. Controls how missing values are reported.
#'   \itemize{
#'     \item 'ifany': Report missing values only if they are present in the data.
#'     \item 'no': Do not report missing values.
#'     \item 'always': Always include a row for missing values, even if none are present.
#'   }
#'   Defaults to 'ifany'.
#' @param missing_text String. The label to use for missing values in the output table. Defaults to 'No/Missing data'.
#' @param missing_stat String. The format for reporting missing value statistics.
#'   \itemize{
#'     \item 'n_percent': Displays count and percentage of missing values (e.g., '10 (5.0%)').
#'     \item 'n': Displays only the count of missing values.
#'     \item 'percent': Displays only the percentage of missing values.
#'   }
#'   Defaults to 'n_percent'.
#' @param include_missing_in_splits Boolean. If `TRUE` (default), includes the missing values (missing text specified in `missing_text`) in the column after specifying `split_by` and `strata_by`.
#' @param sort_categorical_variables_by String. Determines the sorting order for levels within categorical variables.
#'   \itemize{
#'     \item 'alphanumeric': Sorts levels alphabetically or numerically. For custom sorting, use `factor(data, levels = c(level1, level2, ...))`.
#'     \item 'frequency': Sorts levels by their frequency (from highest to lowest).
#'   }
#'   Defaults to 'alphanumeric'.
#' @param calc_percent_by String. Specifies how percentages are calculated for categorical variables.
#'   \itemize{
#'     \item 'column': Percentages are calculated down columns (sum to 100% per column).
#'     \item 'row': Percentages are calculated across rows (sum to 100% per row, useful for `split_by`).
#'     \item 'cell': Percentages are calculated based on the total number of observations in the table.
#'   }
#'   Defaults to 'column' if `split_by` is not `NULL`, otherwise, automatically sets to 'row'.
#' @param calc_col_percent_using String. Specifies how the percentages are to be calculated.
#'   \itemize{
#'     \item 'n_in_column': Percentages are calculated based on the actual number of data points appearing in the column.
#'     \item 'n_valid_in_column': Percentages are calculated based on the valid (i.e., non missing) values.
#'   }
#'   Defaults to 'n_valid_in_column'.
#'   Note: This only applies when 'calc_percent_by = "column"'.
#' @param add_inferential_pvalues Boolean. If `TRUE`, adds the common comparative statistical tests such as t-tests, chi-square test, ANOVA, and their nonparametric counterparts.
#' @param force_statistical_test List. A list of the format `list("variable_name1" ~ "statistical_test1", ...)`, used to specify the test for the specific variable. Below are the list of tests. See type and enter ?gtsummary::tests for more details, and see `tbl_summary() %>% add_p()` section.
#'   \itemize{
#'     \item 't.test': Perform t-test.
#'     \item 'paired.t.test': Perform Paired t-test.
#'     \item 'wilcox.test': Perform Wilcoxon rank-sum test/Mann-Whitney U test.
#'     \item 'paired.wilcox.test': Perform Paired Wilcoxon rank-sum test.
#'     \item 'oneway.test': Perform One-way Analysis of Variance (ANOVA).
#'     \item 'kruskal.test': Perform Kruskal-Wallis test.
#'     \item 'chisq.test': Perform chi-square test of independence.
#'     \item 'chisq.test.no.correct': Perform chi-square test of independence. This specifies 'correct = FALSE'.
#'     \item 'fisher.test': Perform Fisher's exact test.
#'     \item 'mcnemar.test': Perform McNemar's test.
#'     \item 'mcnemar.test.wide': Perform McNemar's test.
#'     \item 'prop.test': Perform Test for equality of proportions.
#'     \item 'mood.test': Perform Mood two-sample test of scale.
#'     \item 'lme4': Perform random intercept logistic regression.
#'     \item 'ancova': Perform Analysis of Covariance (ANCOVA).
#'     \item 'emmeans': Perform Estimated Marginal Means or LS-means.
#'   }
#' @param n_digits_pvalues Numeric. Default is 3. The number of digits of p-values if `add_inferential_pvalues` is set to `TRUE`. This can be set to `NULL` which tells the code to automatically format the p-values. If the p-value < 0.10, there will be a maximum of three decimal places, otherwise, there will only be one.
#' @param bold_significant_pvalues Boolean. If `TRUE` (default), automatically detects if `add_inferential_pvalues = TRUE` then bolds the font of p-values under 0.05 (default) unless set in `bold_significant_pvalues_at`.
#' @param bold_significant_pvalues_at Numeric. The threshold for `bold_significant_pvalues` to bold font if significant.
#' @param header String. The header name for the 'Variable' column in the output table. Defaults to 'Variable'.
#' @param bold_labels Boolean. If `TRUE` (default), variable labels in the table will be bolded.
#' @param italicize_levels Boolean. If `TRUE` (default), categorical variable levels in the table will be italicized.
#' @param clean_table Boolean. If `TRUE` (default), cells with a count of 0 are removed (replaced with empty space) from the descriptive table, making it cleaner.
#' @param table_name String or NULL. The table name. Default is 'auto', where the table name depends if the table is descriptive or inferential in nature (has p-values), and automatically creates a table name. Set to `NULL` to have no table name. Can accept a string that will be used as the table name.
#' @param export_to_excel Boolean. If TRUE, exports the results to an Excel file. Defaults to FALSE. This is currently not recommended as the output is not formatted good.
#' @param export_to_html Boolean. If TRUE, exports the results to an HTML file. Defaults to FALSE.
#' @param export_to_png Boolean. If TRUE, exports the results to a PNG file. Defaults to FALSE.
#' @param export_to_pdf Boolean. If TRUE, exports the results to a PDF file. Defaults to FALSE.
#' @param export_to_word Boolean. If TRUE, exports the results to a Word file. Defaults to FALSE.
#' @param export_filename String. The file name of the exported file/s. Multiple exports can be done by setting each to TRUE. The default is "Summary Results" and will remain as it is if `export_filename = NULL` or `export_filename = ""`.
#'
#' @importFrom openxlsx saveWorkbook writeData addWorksheet createWorkbook
#' @importFrom gtsummary label_style_pvalue
#' @import webshot2
#'
#' @returns A `gtsummary` table object (inherits from `gt_tbl`), representing the simple descriptive statistics.
#' @export
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Load necessary library
#' library(gtsummary)
#' library(dplyr)
#' library(magrittr)
#'
#' # Create a sample dataframe
#' set.seed(123)
#' sample_df <- data.frame(
#'   Age = sample(c(18:65, NA), 100,
#'                replace = TRUE, prob = c(rep(1, 48), 5)),
#'   Gender = sample(c("Male", "Female", NA), 100,
#'                   replace = TRUE, prob = c(0.48, 0.48, 0.04)),
#'   Education = sample(c("Elementary", "High School", "Bachelors", "Masters", NA), 100,
#'                      replace = TRUE, prob = c(0.40, 0.30, 0.15, 0.10, 0.05)),
#'   Income = sample(c(seq(5000, 40000, by = 5000), NA), 100,
#'                   replace = TRUE, prob = c(rep(1, 8), 2)),
#'   Smoker = sample(c(1, 0, NA), 100,
#'                   replace = TRUE, prob = c(0.25, 0.65, 0.10)),
#'   Region = sample(c("North", "South", "East", "West", NA), 100,
#'                   replace = TRUE, prob = c(0.24, 0.24, 0.24, 0.24, .04))
#' )
#'
#' # Basic descriptive statistics
#' # CORRECT: Piped data goes to 'df' argument automatically
#' sample_df %>% performDescriptive()
#'
#' # INCORRECT: Causes an error because 'sample_df' is passed to 'summarize_what'
#' # sample_df %>% performDescriptive(sample_df)
#'
#' # Descriptive statistics split by Gender, with custom header and labels
#' sample_df %>%
#'   performDescriptive(
#'     split_by = "Gender",
#'     split_by_header = "Participant Gender",
#'     rename_variables = list(Age ~ "Age (Years)", Income ~ "Annual Income", Smoker ~ "Smoking Status")
#'   )
#'
#' # Descriptive statistics with median and IQR for continuous, and count only for categorical
#' sample_df %>%
#'   performDescriptive(
#'     continuous_statistics = "medianIQR",
#'     categorical_statistics = "n",
#'     n_digits_continuous = c(0, 0),
#'     n_digits_categorical = c(0, 0)
#'   )
#'
#' # Force 'Smoker' to be continuous (though conceptually it's not, for demonstration)
#' # And display missing always (with custom text), with bold labels and italicized levels
#' sample_df %>%
#'   performDescriptive(
#'     force_continuous = "Smoker",
#'     display_missing = "always",
#'     missing_text = "Data Unavailable",
#'     missing_stat = "n",
#'     bold_labels = TRUE,
#'     italicize_levels = TRUE
#'   )
#'
#' # Sort categorical variables by frequency and calculate percentages by row
#' # (not ideal, but for demonstration only)
#' sample_df %>%
#'   performDescriptive(
#'     sort_categorical_variables_by = "frequency",
#'     calc_percent_by = "row"
#'   )
#'
#' # Perform common comparative statistical tests
#' sample_df %>%
#'   performDescriptive(
#'     split_by = "Gender",
#'     add_inferential_pvalues = TRUE
#'   )
#'}
#'

performDescriptive <- function(
    df,
    summarize_what = everything(),
    split_by = NULL,
    split_by_header = NULL,
    strata_by = NULL,
    rename_variables = NULL,
    continuous_statistics = "meanSD",
    categorical_statistics = "n_percent",
    force_continuous = NULL,
    force_categorical = NULL,
    n_digits_continuous = c(2, 2),
    n_digits_categorical = c(0, 2),
    display_missing = "ifany",
    missing_text = "No data/missing",
    missing_stat = "n_percent",
    include_missing_in_splits = TRUE,
    sort_categorical_variables_by = "alphanumeric",
    calc_percent_by = "column",
    calc_col_percent_using = "n_valid_in_column",
    add_inferential_pvalues = FALSE,
    force_statistical_test = NULL,
    n_digits_pvalues = 3,
    bold_significant_pvalues = TRUE,
    bold_significant_pvalues_at = 0.05,
    header = "Variable",
    bold_labels = TRUE,
    italicize_levels = TRUE,
    clean_table = TRUE,
    table_name = "auto",
    export_to_excel = FALSE,
    export_to_html = FALSE,
    export_to_png = FALSE,
    export_to_pdf = FALSE,
    export_to_word = FALSE,
    export_filename = "Summary Results"
) {
  # Ensure df is a data frame
  df <- as.data.frame(df)

  # --- Parameter checks and improvements ---

  # NEW: Check for common piping error that causes Errors 1 & 2.
  # This happens when user calls `my_data %>% performDescriptive(my_data, ...)`
  if (is.data.frame(summarize_what)) {
    stop(paste(
      "Argument 'summarize_what' cannot be a data frame.",
      "This error often happens when piping data and also passing the data frame as an argument.",
      "Incorrect: my_data %>% performDescriptive(my_data, ...)",
      "Correct:   my_data %>% performDescriptive(...)",
      sep = "\n"
    ))
  }

  # Apply force_categorical to split_by and strata_by first if not NULL
  if (!is.null(force_categorical)) {
    for (col_name in force_categorical) {
      # Convert the column to factor
      df[[col_name]] <- as.factor(df[[col_name]])
    }
  }

  # 1. split_by
  if (!is.null(split_by)) {
    if (!is.character(split_by) || length(split_by) != 1) {
      stop("The 'split_by' parameter must be a single string representing a column name.")
    }
    if (!(split_by %in% names(df))) {
      stop(paste0("The 'split_by' parameter ('", split_by, "') is not a valid column name in the data frame."))
    }

    column_data <- df[[split_by]]

    # Check if the column is categorical (factor or character)
    if (!is.factor(column_data) && !is.character(column_data)) {
      stop(paste0("The column '", split_by, "' must be categorical (factor or character) for splitting."))
    }

    # Determine the number of unique levels, handling NA if present
    num_unique_levels <- length(unique(stats::na.omit(column_data)))

    if (num_unique_levels < 2) {
      stop(paste0(
        "The column '", split_by, "' must have 2 or more unique non-missing levels/factors for splitting. ",
        "It currently has only ", num_unique_levels, "."
      ))
    }

    # 1b. If split_by is not NULL and split_by_header is not provided, use split_by's value
    if (is.null(split_by_header)) {
      split_by_header <- split_by
    }
  }

  # 2. strata_by
  if (!is.null(strata_by)) {
    if (!is.character(strata_by) || length(strata_by) != 1) {
      stop("The 'strata_by' parameter must be a single string representing a column name.")
    }
    if (!(strata_by %in% names(df))) {
      stop(paste0("The 'strata_by' parameter ('", strata_by, "') is not a valid column name in the data frame."))
    }

    column_data <- df[[strata_by]]
    if (!is.factor(column_data) && !is.character(column_data)) {
      stop(paste0("The column '", strata_by, "' must be categorical (factor or character) for stratification."))
    }

    num_unique_levels <- length(unique(stats::na.omit(column_data)))
    if (num_unique_levels < 2) {
      stop(paste0(
        "The column '", strata_by, "' must have 2 or more unique non-missing levels/factors for stratification. ",
        "It currently has only ", num_unique_levels, "."
      ))
    }
  }

  # include_missing_in_splits
  if (!is.logical(include_missing_in_splits)) {
    stop("The 'include_missing_in_splits' parameter must be a single logical value (TRUE or FALSE).")
  }

  # calc_col_percent_using
  valid_method <- c("n_in_column", "n_valid_in_column")
  if (!(calc_col_percent_using %in% valid_method)) {
    stop(paste0("The 'calc_col_percent_using' parameter only accepts: '", paste(valid_method, collapse = "', '"), "'."))
  }

  # Handle missing values in split_by and strata_by columns based on 'include_missing_in_splits'
  if (include_missing_in_splits) {
    # Include split_by and strata_by columns when replacing missing values
    if (calc_percent_by == "column") {
      if (calc_col_percent_using == "n_in_column") {
        if (add_inferential_pvalues) {
          warning("`calc_col_percent_using = 'n_in_column'` is ignored when `add_inferential_pvalues = TRUE` because percentages must reflect valid data only. Missing values will be excluded.")
          calc_col_percent_using <- "n_valid_in_column"
        } else {
          df <- df %>%
            dplyr::mutate(across(where(is.character) | where(is.factor),
                                 ~forcats::fct_explicit_na(factor(.), na_level = missing_text)))
          #warning("Missing values replaced with `missing_text`, INCLUDING `split_by` and `strata_by` columns.")
        }
      } else if (calc_col_percent_using == "n_valid_in_column") {
        #warning("Data includes missing values, INCLUDING `split_by` and `strata_by` columns.")
      }
    }
  } else {
    # Exclude split_by and strata_by columns from missing value replacement
    exclude_cols <- unique(c(split_by, strata_by)) %>% .[!is.null(.)]

    if (calc_percent_by == "column") {
      if (calc_col_percent_using == "n_in_column") {
        if (add_inferential_pvalues) {
          warning("`calc_col_percent_using = 'n_in_column'` is ignored when `add_inferential_pvalues = TRUE`. Percentages will be based on valid data only.")
          calc_col_percent_using <- "n_valid_in_column"
        } else {
          df <- df %>%
            dplyr::mutate(across((where(is.character) | where(is.factor)) & !any_of(exclude_cols),
                                 ~forcats::fct_explicit_na(factor(.), na_level = missing_text)))
        }
        #warning("Missing values replaced with `missing_text`, EXCLUDING `split_by` and `strata_by` columns.")
      } else if (calc_col_percent_using == "n_valid_in_column") {
        #warning("Data includes missing values, EXCLUDING `split_by` and `strata_by` columns.")
      }
    }
  }


  # 3. rename_variables
  if (!is.null(rename_variables)) {
    if (!is.list(rename_variables) || any(!sapply(rename_variables, inherits, "formula"))) {
      stop("The 'rename_variables' parameter only accepts a list of formulas (e.g., list(original_variable_name ~ new_variable_name)).")
    }
    # Check if original variable names in rename_variables exist in df
    original_vars_in_label <- sapply(rename_variables, function(x) as.character(x)[[2]])
    if (!all(original_vars_in_label %in% names(df))) {
      missing_vars <- original_vars_in_label[!(original_vars_in_label %in% names(df))]
      stop(paste0("The following variables specified in 'rename_variables' do not exist in the dataframe: ", paste(missing_vars, collapse = ", "), "."))
    }
  }

  # 4. continuous_statistics
  valid_continuous_stats <- c("meanSD", "meanSD2", "medianIQR", "mean", "sd", "median", "p0", "p25", "p50", "p75", "p100", "IQR")

  if (is.character(continuous_statistics)) {
    if (!all(continuous_statistics %in% valid_continuous_stats)) {
      stop(paste0(
        "The 'continuous_statistics' parameter only accepts: ",
        paste(valid_continuous_stats, collapse = ", "), "."
      ))
    }

    # Construct gtsummary statistic string
    stat_exprs <- c(
      meanSD    = "{mean} ± {sd}",
      meanSD2   = "{mean} ({sd})",
      medianIQR = "{median} ({p25}, {p75})",
      mean      = "{mean}",
      sd        = "{sd}",
      median    = "{median}",
      p0        = "{p0}",
      p25       = "{p25}",
      p50       = "{p50}",
      p75       = "{p75}",
      p100      = "{p100}",
      IQR       = "{p25}, {p75}"
    )

    gts_continuous_stat <- paste(stat_exprs[continuous_statistics], collapse = ", ")
  } else {
    stop("The 'continuous_statistics' parameter must be a character vector of valid summary keywords.")
  }


  # 5. categorical_statistics
  valid_categorical_stats <- c("n_percent", "n", "percent")
  if (!(categorical_statistics %in% valid_categorical_stats)) {
    stop(paste0("The 'categorical_statistics' parameter only accepts: '", paste(valid_categorical_stats, collapse = "', '"), "'."))
  }
  # Map to gtsummary format
  gts_categorical_stat <- switch(categorical_statistics,
                                 "n_percent" = "{n} ({p}%)",
                                 "n"     = "{n}",
                                 "percent"   = "{p}%")

  # 6 and 7. force_continuous and force_categorical
  # Check if they are character vectors and if column names exist in df
  if (!is.null(force_continuous)) {
    if (!is.character(force_continuous) || !all(force_continuous %in% names(df))) {
      stop("The 'force_continuous' parameter must be a character vector of valid column names in the dataframe.")
    }
  }
  if (!is.null(force_categorical)) {
    if (!is.character(force_categorical) || !all(force_categorical %in% names(df))) {
      stop("The 'force_categorical' parameter must be a character vector of valid column names in the dataframe.")
    }
  }
  # Ensure no overlap between force_continuous and force_categorical
  if (!is.null(force_continuous) && !is.null(force_categorical)) {
    overlap_vars <- intersect(force_continuous, force_categorical)
    if (length(overlap_vars) > 0) {
      stop(paste0("The following variables are specified in both 'force_continuous' and 'force_categorical': ", paste(overlap_vars, collapse = ", "), ". A variable cannot be forced into both types simultaneously."))
    }
  }

  # Convert specified 'force_continuous' columns to numeric
  if (!is.null(force_continuous)) {
    for (col_name in force_continuous) {
      if (!is.numeric(df[[col_name]])) {
        # Attempt conversion
        converted_col <- as.numeric(df[[col_name]])
        # Check if conversion introduced NAs that weren't originally there
        if (any(is.na(converted_col) & !is.na(df[[col_name]]))) {
          warning(paste0("Column '", col_name, "' was specified as 'force_continuous' but contains non-numeric values that were converted to NA. Please inspect your data."))
        }
        df[[col_name]] <- converted_col
      }
    }
  }

  # 8 and 9. n_digits_continuous and n_digits_categorical
  if (
    !(
      (is.numeric(n_digits_continuous) && # Is it numeric?
       length(n_digits_continuous) == 2 && # Does it have exactly two elements?
       all(n_digits_continuous >= 0))     || # Are all elements non-negative? (use all() for vector)
      (is.character(n_digits_continuous) && # Is it a character string?
       length(n_digits_continuous) == 1 && # Is it a single string?
       n_digits_continuous == "dynamic")    # Is that string "dynamic"?
    )
  ) {
    stop("The 'n_digits_continuous' parameter must be a numeric vector of two non-negative integers (e.g., c(1, 1)), representing the number of decimal places of mean and SD, respectively, and median and IQR, respectively. The other option is 'dynamic' which selects the number of decimal places dynamically.")
  }
  if (
    !(
      (is.numeric(n_digits_categorical) && # Is it numeric?
       length(n_digits_categorical) == 2 && # Does it have exactly two elements?
       all(n_digits_categorical >= 0))     || # Are all elements non-negative? (use all() for vector)
      (is.character(n_digits_categorical) && # Is it a character string?
       length(n_digits_categorical) == 1 && # Is it a single string?
       n_digits_categorical == "dynamic")    # Is that string "dynamic"?
    )
  ) {
    stop("The 'n_digits_categorical' parameter must be a numeric vector of two non-negative integers (e.g., c(0, 1)), representing the number of decimal places of frequencies or counts and percentages, respectively.  The other option is 'dynamic' which selects the number of decimal places dynamically.")
  }

  # Set to NULL if "dynamic" is chosen
  if (isTRUE(identical(n_digits_continuous, "dynamic"))) {
    n_digits_continuous <- NULL
  }

  if (isTRUE(identical(n_digits_categorical, "dynamic"))) {
    n_digits_categorical <- NULL
  }

  # 10. display_missing
  valid_display_missing <- c("ifany", "no", "always")
  if (!(display_missing %in% valid_display_missing)) {
    stop(paste0("The 'display_missing' parameter only accepts: '", paste(valid_display_missing, collapse = "', '"), "'."))
  }

  # 11. missing_text is already a string, no specific check needed beyond default type.

  # 12. missing_stat (Updated for uniformity)
  valid_missing_stat <- c("n_percent", "n", "percent")
  if (!(missing_stat %in% valid_missing_stat)) {
    stop(paste0("The 'missing_stat' parameter only accepts: '", paste(valid_missing_stat, collapse = "', '"), "'."))
  }
  # Map to gtsummary format for missing_stat
  gts_missing_stat <- switch(missing_stat,
                             "n_percent" = "{N_miss} ({p_miss}%)",
                             "n"     = "{N_miss}",
                             "percent"   = "{p_miss}%")

  # 13. sort_categorical_variables_by
  valid_sort_cat_by <- c("alphanumeric", "frequency")
  if (!(sort_categorical_variables_by %in% valid_sort_cat_by)) {
    stop(paste0("The 'sort_categorical_variables_by' parameter only accepts: '", paste(valid_sort_cat_by, collapse = "', '"), "'."))
  }

  # 14. calc_percent_by
  valid_calc_percent_by <- c("column", "row", "cell")
  if (!(calc_percent_by %in% valid_calc_percent_by)) {
    stop(paste0("The 'calc_percent_by' parameter only accepts: '", paste(valid_calc_percent_by, collapse = "', '"), "'."))
  }

  # 15. header is already a string, no specific check needed.

  # 16. clean_table
  if (!is.logical(clean_table) || length(clean_table) != 1) {
    stop("The 'clean_table' parameter must be a single logical value (TRUE or FALSE).")
  }

  # 17 and 18. bold_labels and italicize_levels
  if (!is.logical(bold_labels) || length(bold_labels) != 1) {
    stop("The 'bold_labels' parameter must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(italicize_levels) || length(italicize_levels) != 1) {
    stop("The 'italicize_levels' parameter must be a single logical value (TRUE or FALSE).")
  }

  # 19, add_inferential_pvalues
  if (!is.logical(add_inferential_pvalues)) {
    stop("The 'add_inferential_pvalues' parameter must be a logical value (TRUE or FALSE).")
  }
  if (add_inferential_pvalues) {
    if (!is.null(strata_by)) {
      stop("The 'add_inferential_pvalues' cannot be TRUE when 'strata_by' is present. Set 'strata_by = NULL' if p-values are still desired.")
    }
  }

  # force_statistical_test
  if (!is.null(force_statistical_test)) {
    if (!is.list(force_statistical_test)) {
      stop("The parameter 'force_statistical_test' must be a list of the format 'list(variable_name1 ~ statistical_test1, ...)'.")
    }
  }

  # 20, n_digits_pvalues
  if (!is.null(n_digits_pvalues) && (!is.numeric(n_digits_pvalues) || n_digits_pvalues <= 0)) {
    stop("The 'n_digits_pvalues' parameter must be a positive numeric value or NULL.")
  }

  # 21, bold_significant_pvalues
  if (!is.logical(bold_significant_pvalues)) {
    stop("The 'bold_significant_pvalues' parameter must be a logical value (TRUE or FALSE).")
  }
  if (!is.null(strata_by)) {
    if (add_inferential_pvalues) {
      if (bold_significant_pvalues) {
        stop("The 'bold_significant_pvalues' parameter cannot be used when 'strata_by' is specified. Set 'strata_by = NULL' if p-values are still desired.")
      }
    } else {
      if (bold_significant_pvalues) {
        message("The 'bold_significant_pvalues = TRUE' (default or inputted) parameter is ignored since 'add_inferential_pvalues = FALSE'.")
      }
    }
  }

  # 22, bold_significant_pvalues_at
  if (!is.numeric(bold_significant_pvalues_at)) {
    stop("The 'bold_significant_pvalues_at' parameter must be a decimal/numeric value where bold p-values will be applied. Default is 0.05 (common). The other common threshold is 0.1, where p-values below 0.1 are bolded.")
  }

  # --- Analysis using gtsummary ---

  # Set 'calc_percent_by' automatically to 'row' if 'split_by' is not NULL
  if (is.null(calc_percent_by)) {
    if (!is.null(split_by)) {
      calc_percent_by = 'row'
    }
  }

  # Define 'type' argument for tbl_summary
  # OPTIMIZED: Replaced for-loop with lapply for vectorized list creation.
  type_list <- if (!is.null(force_continuous)) {
    lapply(force_continuous, function(col_name) {
      # Enclose column names in backticks if they contain spaces or special characters
      as.formula(paste0("`", col_name, "` ~ 'continuous'"))
    })
  } else {
    list()
  }

  # Commented since this is already done above after the transformation of df to data frame. Using this code will result to an error, where the force_categorical variables will be missing(?)
  if (!is.null(force_categorical)) {
    type_list_categorical <- lapply(force_categorical, function(col_name) {
      as.formula(paste0("`", col_name, "` ~ 'categorical'"))
    })
    type_list <- c(type_list, type_list_categorical) # Combine lists
  }

  # table_name
  if (!is.character(table_name) && !is.null(table_name)) {
    stop("The parameter 'table_name' must be a string or NULL.")

  }

  # Helper function to create summary table
  # NOTE: This internal function is required for the 'tbl_strata' .tbl_fun argument
  # and to avoid duplicating the large 'tbl_summary' call, per user request to avoid subfunctions.
  # This is a case where it is necessary for functionality.
  create_summary <- function(data) {
    gtsummary::tbl_summary(
      data = data,
      include = summarize_what,
      by = split_by,
      label = rename_variables,
      type = type_list,
      statistic = list(
        all_continuous() ~ gts_continuous_stat,
        all_categorical() ~ gts_categorical_stat
      ),
      digits = list(
        all_categorical() ~ n_digits_categorical,
        all_continuous() ~ n_digits_continuous
      ),
      missing = display_missing,
      missing_text = missing_text,
      missing_stat = gts_missing_stat,
      sort = list(all_categorical() ~ sort_categorical_variables_by),
      percent = calc_percent_by
    ) %>%
      gtsummary::modify_header(label = paste0("**", header, "**"))
  }

  # # Replace all missing values with NA if 'calc_col_percent_using' = 'n_in_column' because this will make the percentages equal to 100 with the missing values on.
  # if (calc_col_percent_using == "n_in_column") {
  #   warning("")
  #   df[df == missing_text] <- NA
  #
  # }

  # Main logic
  result_table <- if (!is.null(strata_by)) {
    df %>%
      dplyr::filter(!is.na(.data[[strata_by]])) %>%
      gtsummary::tbl_strata(
        strata = strata_by,
        include = summarize_what,
        .tbl_fun = ~ create_summary(.x)
      )
  } else {
    create_summary(df)
  }

  # Modify spanning header for split_by (only if strata_by is NULL)
  if (is.null(strata_by) && !is.null(split_by) && !is.null(split_by_header)) {
    result_table <- result_table %>%
      gtsummary::modify_spanning_header(
        gtsummary::all_stat_cols() ~ paste0("**", split_by_header, "**")
      )
  }

  # Apply bold_labels and italicize_levels
  if (bold_labels) {
    result_table <- result_table %>%
      gtsummary::bold_labels()
  }
  if (italicize_levels) {
    result_table <- result_table %>%
      gtsummary::italicize_levels()
  }

  # For the force_statistical_test
  # OPTIMIZED: Removed redundant line: 'if (!is.null(force_statistical_test)) { force_statistical_test = force_statistical_test }'

  # Add common statistical test p-values (and whether to bold it, and where)
  if (add_inferential_pvalues) {
    result_table <- result_table %>%
      gtsummary::add_p(
        test = force_statistical_test,
        pvalue_fun = label_style_pvalue(digits = n_digits_pvalues
        )) # Also adjust the no. of decimals of p-values
  }
  if (add_inferential_pvalues) {
    if (bold_significant_pvalues) { # Automatically bold p-values if less than 'bold_significant_pvalues_at' if either NULL or TRUE
      result_table <- result_table %>%
        gtsummary::bold_p(t = bold_significant_pvalues_at) # Also add where to add the bold
    }
  }

  # Clean table (remove "0 (0%)" or similar cells) using modify_table_body
  if (clean_table) {
    # Dynamically generate the "zero" string based on n_digits_categorical for percentages
    # n_digits_categorical[1] is for count, n_digits_categorical[2] is for percentage
    zero_percent_format <- paste0("0", ifelse(n_digits_categorical[2] > 0, paste0(".", paste(rep("0", n_digits_categorical[2]), collapse = "")), ""), "%")

    # Possible zero strings depending on categorical_statistics and n_digits_categorical
    # These are the exact strings gtsummary will produce when count or percent is zero.
    zero_strings_to_replace <- c()
    if (categorical_statistics == "n_percent") {
      zero_strings_to_replace <- c(zero_strings_to_replace, paste0("0 (", zero_percent_format, ")"))
    } else if (categorical_statistics == "n") {
      zero_strings_to_replace <- c(zero_strings_to_replace, "0")
    } else if (categorical_statistics == "percent") {
      zero_strings_to_replace <- c(zero_strings_to_replace, zero_percent_format)
    }

    # Also consider the missing_stat for consistency in cleaning "0" values if it's "n" or "percent"
    if (missing_stat == "n_percent") {
      # Note: gtsummary's missing_stat uses {p_miss} which defaults to 1 decimal place.
      # This can be hardcoded or made dynamic if gtsummary allowed setting digits for missing_stat.
      # For now, we'll anticipate "0 (0.0%)" and "0 (0%)" and "0.0%"
      zero_strings_to_replace <- c(zero_strings_to_replace, "0 (0.0%)", "0 (0%)")
    } else if (missing_stat == "n") {
      zero_strings_to_replace <- c(zero_strings_to_replace, "0")
    } else if (missing_stat == "percent") {
      zero_strings_to_replace <- c(zero_strings_to_replace, "0.0%", "0%") # Anticipate different formatting
    }

    # Remove duplicates
    zero_strings_to_replace <- unique(zero_strings_to_replace)

    result_table <- result_table %>%
      gtsummary::modify_table_body(
        ~ .x %>%
          dplyr::mutate(
            dplyr::across(
              dplyr::starts_with("stat_"), # Selects columns like stat_0, stat_1 etc.
              ~ dplyr::case_when(
                .x %in% zero_strings_to_replace ~ "", # Replace matched strings with empty
                TRUE ~ .x # Keep original for others
              )
            )
          )
      )
  }

  # Add table name header
  if (!is.null(table_name) && length(table_name) == 1 && table_name == "auto") {
    # Ensure 'gt' package is installed
    if (!requireNamespace("gt", quietly = TRUE)) {
      warning("Package 'gt' is required but not installed. Installing now...")
      install.packages("gt")
    }

    # Compose title
    title_prefix <- if (add_inferential_pvalues) {
      "Inferential Statistics Stratified by "
    } else {
      "Descriptive Statistics Stratified by "
    }

    if (!is.null(strata_by) && !is.null(split_by_header)) {
      title_text <- paste0(title_prefix, strata_by, " and ", split_by_header, " Columns")
    } else if (is.null(strata_by) && !is.null(split_by_header)) {
      title_text <- paste0(title_prefix, split_by_header, " Column")
    } else {
      title_text <- if (add_inferential_pvalues) {
        "Inferential Statistics"
      } else {
        "Descriptive Statistics"
      }
    }

    result_table <- result_table %>%
      gtsummary::as_gt() %>%
      gt::tab_header(title = gt::md(paste0("**", title_text, "**")))

  } else if (is.character(table_name) && length(table_name) == 1) {
    # Custom table name provided
    result_table <- result_table %>%
      gtsummary::as_gt() %>%
      gt::tab_header(title = gt::md(paste0("**", table_name, "**")))
  }

  # Export

  # Fixing filename
  if (is.null(export_filename) || export_filename == "") {
    export_filename = "Summary Results"
  }

  # To Excel
  if (export_to_excel) {
    # Ensure 'openxlsx' package is installed
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("Package 'openxlsx' is required but not installed. Installing now...")
      install.packages("openxlsx")
    }

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = export_filename)
    writeData(wb, export_filename, result_table)
    saveWorkbook(wb, paste0(export_filename, ".xlsx"), overwrite = TRUE)

    message(paste0("The results have been saved to: ", normalizePath(file.path(getwd(), paste0(export_filename, ".xlsx")))))
  }

  # To HTML
  if (export_to_html) {
    result_table %>%
      gt::gtsave(paste0(export_filename, ".html"))

    message(paste0("The results have been saved to: ", normalizePath(file.path(getwd(), paste0(export_filename, ".html")))))
  }

  # To PNG
  if (export_to_png) {
    # Ensure 'webshot2' package is installed
    if (!requireNamespace("webshot2", quietly = TRUE)) {
      warning("Package 'webshot2' is required but not installed. Installing now...")
      install.packages("webshot2")
    }

    result_table %>%
      gt::gtsave(paste0(export_filename, ".png"))

    message(paste0("The results have been saved to: ", normalizePath(file.path(getwd(), paste0(export_filename, ".png")))))
  }

  # To PDF
  if (export_to_pdf) {
    result_table %>%
      gt::gtsave(paste0(export_filename, ".pdf"))

    message(paste0("The results have been saved to: ", normalizePath(file.path(getwd(), paste0(export_filename, ".pdf")))))
  }

  # To Word
  if (export_to_word) {

    # Delete the file if it already exists
    # This is done since there's an error if the ".docx" file is already created
    if (file.exists(paste0(export_filename, ".docx"))) {
      file.remove(paste0(export_filename, ".docx"))
    }

    result_table %>%
      gt::gtsave(paste0(export_filename, ".docx"))

    message(paste0("The results have been saved to: ", normalizePath(file.path(getwd(), paste0(export_filename, ".docx")))))
  }

  # Final messages
  message("\nThe analysis has been completed.")

  if (export_to_excel | export_to_html | export_to_png | export_to_pdf | export_to_word) {
    message("\nExport/s has/have been generated. Check your current working directory.")
  }

  return(result_table) # Always return the gtsummary/gt table object
}
