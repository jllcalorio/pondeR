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
#' @param label List. A list of the format `list(original_variable_name1 ~ new_variable_name1, ...)`, used to rename variables in the output table.
#' @param continuous_statistics String. The type of statistics to report for continuous variables. Options are:
#'   \itemize{
#'     \item 'meanSD': mean ± standard deviation
#'     \item 'meanSD2': mean (standard deviation), which is another format of meanSD
#'     \item 'medianIQR': median (interquartile range) or median (p25, p75)
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
#'     \item 'n_pct': Displays count and percentage of valid values (e.g., '10 (5.0%)').
#'     \item 'n': Displays only the count of valid values.
#'     \item 'pct': Displays only the percentage of valid values.
#'     }
#'     Defaults to 'n_pct'.
#' Options are 'n_pct' (count (percentage)), 'n' (count only), or 'pct' (percentage only). Defaults to 'n_pct'.
#' @param force_continuous Vector. A character vector of column names that should be treated as continuous, even if their data type suggests otherwise (e.g., numeric variables with few unique values that might be misclassified as categorical).
#' @param force_categorical Vector. A character vector of column names that should be treated as categorical, even if their data type suggests otherwise (e.g., numeric variables that are truly categories, or binary variables that should be presented as full categories rather than dichotomous).
#' @param n_digits_continuous Numeric Vector. A vector of two integers specifying the number of digits to round for mean and standard deviation (or median and IQR) respectively, for continuous variables. For example, `c(2, 2)` means 2 decimal places for both (default).
#' @param n_digits_categorical Numeric Vector. A vector of two integers specifying the number of digits to round for counts and percentages respectively, for categorical variables. For example, `c(0, 2)` means 0 decimal places for counts and 2 for percentages (default).
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
#'     \item 'n_pct': Displays count and percentage of missing values (e.g., '10 (5.0%)').
#'     \item 'n': Displays only the count of missing values.
#'     \item 'pct': Displays only the percentage of missing values.
#'   }
#'   Defaults to 'n_pct'.
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
#' @param header String. The header name for the 'Variable' column in the output table. Defaults to 'Variable'.
#' @param clean_table Boolean. If `TRUE` (default), cells with a count of 0 are removed (replaced with empty space) from the descriptive table, making it cleaner.
#' @param bold_labels Boolean. If `TRUE` (default), variable labels in the table will be bolded.
#' @param italicize_levels Boolean. If `TRUE` (default), categorical variable levels in the table will be italicized.
#' @param add_inferential_pvalues Boolan. If `TRUE`, adds the common comparative statistical tests such as t-tests, chi-square test, ANOVA, and their nonparametric counterparts.
#' @param n_digits_pvalues Numeric. The number of digits of p-values if `add_inferential_pvalues` is set to `TRUE`.
#' @param bold_significant_pvalues Boolean/NULL. If `TRUE` or `NULL`, automatically detects if `add_inferential_pvalues = TRUE` then bolds the font of p-values under 0.05 (default) unless set in `bold_significant_pvalues_at`.
#' @param bold_significant_pvalues_at Numeric. The threshold for `bold_significant_pvalues` to bold font if significant.
#'
#' @returns A `gtsummary` table object (inherits from `gt_tbl`), representing the simple descriptive statistics.
#' @export
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
#'   Age = rnorm(100, 35, 10),
#'   Gender = sample(c("Male", "Female", NA), 100, replace = TRUE, prob = c(0.48, 0.48, 0.04)),
#'   Education = sample(c("High School", "Bachelors", "Masters", NA), 100, replace = TRUE, prob = c(0.25, 0.35, 0.35, 0.05)),
#'   Income = rgamma(100, shape = 5, rate = 0.1) * 1000,
#'   Smoker = sample(c(1, 0), 100, replace = TRUE, prob = c(0.3, 0.7)),
#'   Region = sample(c("North", "South", "East", "West"), 100, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
#' )
#'
#' # Basic descriptive statistics
#' performDescriptive(df = sample_df)
#'
#' # Descriptive statistics split by Gender, with custom header and labels
#' performDescriptive(
#'   df = sample_df,
#'   split_by = "Gender",
#'   split_by_header = "Participant Gender",
#'   label = list(Age ~ "Age (Years)", Income ~ "Annual Income", Smoker ~ "Smoking Status")
#' )
#'
#' # Descriptive statistics with median and IQR for continuous, and count only for categorical
#' performDescriptive(
#'   df = sample_df,
#'   continuous_statistics = "medianIQR",
#'   categorical_statistics = "n",
#'   n_digits_continuous = c(0, 0),
#'   n_digits_categorical = c(0, 0)
#' )
#'
#' # Force 'Smoker' to be continuous (though conceptually it's not, for demonstration)
#' # And display missing always (with custom text), with bold labels and italicized levels
#' performDescriptive(
#'   df = sample_df,
#'   force_continuous = "Smoker",
#'   display_missing = "always",
#'   missing_text = "Data Unavailable",
#'   missing_stat = "n",
#'   bold_labels = TRUE,
#'   italicize_levels = TRUE
#' )
#'
#' # Sort categorical variables by frequency and calculate percentages by row (not ideal, but for demonstration only)
#' performDescriptive(
#'   df = sample_df,
#'   sort_categorical_variables_by = "frequency",
#'   calc_percent_by = "row"
#' )
#' }
#'

performDescriptive <- function(
    df,
    summarize_what = everything(),
    split_by = NULL,
    split_by_header = NULL,
    strata_by = NULL,
    label = NULL,
    continuous_statistics = "meanSD",
    categorical_statistics = "n_pct",
    force_continuous = NULL,
    force_categorical = NULL,
    n_digits_continuous = c(2, 2),
    n_digits_categorical = c(0, 2),
    display_missing = "ifany",
    missing_text = "No data/missing",
    missing_stat = "n_pct",
    sort_categorical_variables_by = "alphanumeric",
    calc_percent_by = "column",
    header = "Variable",
    clean_table = TRUE,
    bold_labels = TRUE,
    italicize_levels = TRUE,
    add_inferential_pvalues = FALSE,
    n_digits_pvalues = 3,
    bold_significant_pvalues = NULL,
    bold_significant_pvalues_at = 0.05
) {
  # Load necessary packages
  if (!requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Package 'gtsummary' is required but not installed. Please install it with install.packages('gtsummary').")
  }
  # dplyr is needed for `across` and other tidyverse functions, so load it
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed. Please install it with install.packages('dplyr').")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' is required for the pipe operator (%>%) but not installed. Please install it with install.packages('magrittr').")
  }

  # Attach necessary packages
  library(magrittr)
  library(dplyr) # Explicitly load dplyr for `across`, `mutate`, `starts_with`

  # Ensure df is a data frame
  df <- as.data.frame(df)

  # --- Parameter checks and improvements ---

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

  # 3. label
  if (!is.null(label)) {
    if (!is.list(label) || any(!sapply(label, inherits, "formula"))) {
      stop("The 'label' parameter only accepts a list of formulas (e.g., list(original_variable_name ~ new_variable_name)).")
    }
    # Check if original variable names in label exist in df
    original_vars_in_label <- sapply(label, function(x) as.character(x)[[2]])
    if (!all(original_vars_in_label %in% names(df))) {
      missing_vars <- original_vars_in_label[!(original_vars_in_label %in% names(df))]
      stop(paste0("The following variables specified in 'label' do not exist in the dataframe: ", paste(missing_vars, collapse = ", "), "."))
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
  valid_categorical_stats <- c("n_pct", "n", "pct")
  if (!(categorical_statistics %in% valid_categorical_stats)) {
    stop(paste0("The 'categorical_statistics' parameter only accepts: '", paste(valid_categorical_stats, collapse = "', '"), "'."))
  }
  # Map to gtsummary format
  gts_categorical_stat <- switch(categorical_statistics,
                                 "n_pct" = "{n} ({p}%)",
                                 "n"     = "{n}",
                                 "pct"   = "{p}%")

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
  if (!is.numeric(n_digits_continuous) || length(n_digits_continuous) != 2 || any(n_digits_continuous < 0)) {
    stop("The 'n_digits_continuous' parameter must be a numeric vector of two non-negative integers (e.g., c(1, 1)).")
  }
  if (!is.numeric(n_digits_categorical) || length(n_digits_categorical) != 2 || any(n_digits_categorical < 0)) {
    stop("The 'n_digits_categorical' parameter must be a numeric vector of two non-negative integers (e.g., c(0, 1)).")
  }

  # 10. display_missing
  valid_display_missing <- c("ifany", "no", "always")
  if (!(display_missing %in% valid_display_missing)) {
    stop(paste0("The 'display_missing' parameter only accepts: '", paste(valid_display_missing, collapse = "', '"), "'."))
  }

  # 11. missing_text is already a string, no specific check needed beyond default type.

  # 12. missing_stat (Updated for uniformity)
  valid_missing_stat <- c("n_pct", "n", "pct")
  if (!(missing_stat %in% valid_missing_stat)) {
    stop(paste0("The 'missing_stat' parameter only accepts: '", paste(valid_missing_stat, collapse = "', '"), "'."))
  }
  # Map to gtsummary format for missing_stat
  gts_missing_stat <- switch(missing_stat,
                             "n_pct" = "{N_miss} ({p_miss}%)", # Consistent with categorical "n_pct"
                             "n"     = "{N_miss}",
                             "pct"   = "{p_miss}%") # Consistent with categorical "pct"

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

  # 20, n_digits_pvalues
  if (!is.numeric(n_digits_pvalues)) {
    stop("The 'n_digits_pvalues' parameter must be a numeric value.")
  }

  # 21, bold_significant_pvalues
  if (!is.null(bold_significant_pvalues) && !is.logical(bold_significant_pvalues)) {
    stop("The 'bold_significant_pvalues' parameter must be a logical value (TRUE or FALSE) or NULL.")
  }
  if (!is.null(bold_significant_pvalues)) { # If not NULL, do something
    if (bold_significant_pvalues == TRUE && !is.null(split_by)) {
      message("Warning: The 'bold_significant_pvalues' parameter cannot be used when 'strata_by' is specified. Set 'strata_by = NULL' if p-values are still desired.")
    }
  }

  # 22, bold_significant_pvalues_at
  if (!is.numeric(bold_significant_pvalues_at)) {
    stop("The 'bold_significant_pvalues_at' parameter must be a decimal/numeric value.")
  }
  if (is.null(split_by)) {
    if (is.null(bold_significant_pvalues) || bold_significant_pvalues == TRUE) {
      message("Warning: The 'bold_significant_pvalues_at' parameter cannot be used when 'strata_by' is specified. Set 'strata_by = NULL' if p-values are still desired.")
    }
  }


  # --- Analysis using gtsummary ---

  # Set 'calc_percent_by' automatically to 'row' if 'split_by' is not NULL
  if (is.null(calc_percent_by)) {
    if (!is.null(split_by)) {
      calc_percent_by = 'row'
    }
    print(calc_percent_by)
  }

  # Define 'type' argument for tbl_summary
  type_list <- list()

  if (!is.null(force_continuous)) {
    for (col_name in force_continuous) {
      # Enclose column names in backticks if they contain spaces or special characters
      type_list[[length(type_list) + 1]] <- as.formula(paste0("`", col_name, "` ~ 'continuous'"))
    }
  }

  if (!is.null(force_categorical)) {
    for (col_name in force_categorical) {
      # Enclose column names in backticks if they contain spaces or special characters
      type_list[[length(type_list) + 1]] <- as.formula(paste0("`", col_name, "` ~ 'categorical'"))
    }
  }

  # Helper function to create summary table
  create_summary <- function(data) {
    gtsummary::tbl_summary(
      data = data,
      include = summarize_what,
      by = split_by,
      label = label,
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

  # Add common statistical test p-values (and whether to bold it, and where)
  if (add_inferential_pvalues) {
    result_table <- result_table %>%
      gtsummary::add_p(
        pvalue_fun = label_style_pvalue(
          digits = n_digits_pvalues # Also adjust the no. of decimals of p-values
        ))
  }
  if (add_inferential_pvalues) {
    if (is.null(bold_significant_pvalues) || bold_significant_pvalues == TRUE) { # Automatically bold p-values if less than 'bold_significant_pvalues_at' if either NULL or TRUE
      result_table <- result_table %>%
        gtsummary::bold_p(t = bold_significant_pvalues_at) # Also add where to add the bold
    }
  }

  # Clean table (remove "0 (0%)" or similar cells) using modify_table_body
  if (clean_table) {
    # Dynamically generate the "zero" string based on n_digits_categorical for percentages
    # n_digits_categorical[1] is for count, n_digits_categorical[2] is for percentage
    zero_pct_format <- paste0("0", ifelse(n_digits_categorical[2] > 0, paste0(".", paste(rep("0", n_digits_categorical[2]), collapse = "")), ""), "%")

    # Possible zero strings depending on categorical_statistics and n_digits_categorical
    # These are the exact strings gtsummary will produce when count or percent is zero.
    zero_strings_to_replace <- c()
    if (categorical_statistics == "n_pct") {
      zero_strings_to_replace <- c(zero_strings_to_replace, paste0("0 (", zero_pct_format, ")"))
    } else if (categorical_statistics == "n") {
      zero_strings_to_replace <- c(zero_strings_to_replace, "0")
    } else if (categorical_statistics == "pct") {
      zero_strings_to_replace <- c(zero_strings_to_replace, zero_pct_format)
    }

    # Also consider the missing_stat for consistency in cleaning "0" values if it's "n" or "pct"
    if (missing_stat == "n_pct") {
      # Note: gtsummary's missing_stat uses {p_miss} which defaults to 1 decimal place.
      # This can be hardcoded or made dynamic if gtsummary allowed setting digits for missing_stat.
      # For now, we'll anticipate "0 (0.0%)" and "0 (0%)" and "0.0%"
      zero_strings_to_replace <- c(zero_strings_to_replace, "0 (0.0%)", "0 (0%)")
    } else if (missing_stat == "n") {
      zero_strings_to_replace <- c(zero_strings_to_replace, "0")
    } else if (missing_stat == "pct") {
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

  return(result_table) # Always return the gtsummary/gt table object
}
