# Performs Descriptive Statistics on a Dataframe

This function generates a comprehensive descriptive statistics table
from a given dataframe, optionally stratified by a categorical variable.
It provides customizable summaries for both continuous and categorical
variables, including control over rounding, missing data reporting,
labeling, formatting, and table aesthetics. The output is a
publication-ready `gtsummary` table that can be easily integrated into
reports, manuscripts, or exported (via copy-paste for now) for further
use.

## Usage

``` r
run_summarytable(
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
)
```

## Arguments

- df:

  Dataframe. The input dataframe. The first row is assumed to contain
  column names.

- summarize_what:

  Character vector. A character vector of column names that will be
  included in the analysis, including what will be placed in `split_by`
  and `strata_by`.

- split_by:

  String. A column name in `df` (must be a categorical variable with at
  least 2 unique levels/factors) by which to stratify the descriptive
  statistics. If `NULL` (default), no stratification is performed.

- split_by_header:

  String. The header name to display for the `split_by` variable in the
  output table. Only applicable if `split_by` is not `NULL`, and takes
  on the value of `split_by`. Defaults to `NULL` initially.

- strata_by:

  String. A column name in `df` (must be a categorical variable with at
  least 2 unique levels/factors) by which to stratify the descriptive
  statistics 'before' `split_by`. If `NULL` (default), no stratification
  is performed.

- rename_variables:

  List. A list of the format
  `list("original_variable_name1" ~ "new_variable_name1", ...)`, used to
  rename variables in the output table.

- continuous_statistics:

  String. The type of statistics to report for continuous variables.
  Options are:

  - 'meanSD': mean ± standard deviation

  - 'meanSD2': mean (standard deviation), which is another format of
    meanSD

  - 'medianIQR': median (interquartile range) or equivalently median
    (p25, p75)

  - 'mean': mean (average)

  - 'sd': standard deviation

  - 'p0': 0th percentile

  - 'p25': 25th percentile (1st quartile)

  - 'p50': 50th percentile (2nd quartile) (median)

  - 'p75': 75th percentile (3rd quartile)

  - 'p100': 100th percentile

  - 'IQR': p25, p75 (interquartile range)

  Defaults to 'meanSD'.

- categorical_statistics:

  String. The type of statistics to report for categorical variables.

  - 'n_percent': Displays count and percentage of valid values (e.g.,
    '10 (5.0%)').

  - 'n': Displays only the count of valid values.

  - 'percent': Displays only the percentage of valid values.

  Defaults to 'n_percent'. Options are 'n_percent' (count (percentage)),
  'n' (count only), or 'percent' (percentage only). Defaults to
  'n_percent'.

- force_continuous:

  Vector. A character vector of column names that should be treated as
  continuous, even if their data type suggests otherwise (e.g., numeric
  variables with few unique values that might be misclassified as
  categorical).

- force_categorical:

  Vector. A character vector of column names that should be treated as
  categorical, even if their data type suggests otherwise (e.g., numeric
  variables that are truly categories, or binary variables that should
  be presented as full categories rather than dichotomous).

- n_digits_continuous:

  Numeric Vector or String.

  - c(2, 2): A vector of two integers specifying the number of digits to
    round for mean and standard deviation (or median and IQR)
    respectively, for continuous variables. For example, `c(2, 2)` means
    2 decimal places for both (default).

  - "dynamic": The other option is 'dynamic' which selects the number of
    decimal places dynamically.

- n_digits_categorical:

  Numeric Vector or String.

  - c(0, 2): A vector of two integers specifying the number of digits to
    round for counts and percentages respectively, for categorical
    variables. For example, `c(0, 2)` means 0 decimal places for counts
    and 2 for percentages (default).

  - "dynamic": Same as in `n_digits_continuous`.

- display_missing:

  String. Controls how missing values are reported.

  - 'ifany': Report missing values only if they are present in the data.

  - 'no': Do not report missing values.

  - 'always': Always include a row for missing values, even if none are
    present.

  Defaults to 'ifany'.

- missing_text:

  String. The label to use for missing values in the output table.
  Defaults to 'No/Missing data'.

- missing_stat:

  String. The format for reporting missing value statistics.

  - 'n_percent': Displays count and percentage of missing values (e.g.,
    '10 (5.0%)').

  - 'n': Displays only the count of missing values.

  - 'percent': Displays only the percentage of missing values.

  Defaults to 'n_percent'.

- include_missing_in_splits:

  Boolean. If `TRUE` (default), includes the missing values (missing
  text specified in `missing_text`) in the column after specifying
  `split_by` and `strata_by`.

- sort_categorical_variables_by:

  String. Determines the sorting order for levels within categorical
  variables.

  - 'alphanumeric': Sorts levels alphabetically or numerically. For
    custom sorting, use `factor(data, levels = c(level1, level2, ...))`.

  - 'frequency': Sorts levels by their frequency (from highest to
    lowest).

  Defaults to 'alphanumeric'.

- calc_percent_by:

  String. Specifies how percentages are calculated for categorical
  variables.

  - 'column': Percentages are calculated down columns (sum to 100% per
    column).

  - 'row': Percentages are calculated across rows (sum to 100% per row,
    useful for `split_by`).

  - 'cell': Percentages are calculated based on the total number of
    observations in the table.

  Defaults to 'column' if `split_by` is not `NULL`, otherwise,
  automatically sets to 'row'.

- calc_col_percent_using:

  String. Specifies how the percentages are to be calculated.

  - 'n_in_column': Percentages are calculated based on the actual number
    of data points appearing in the column.

  - 'n_valid_in_column': Percentages are calculated based on the valid
    (i.e., non missing) values.

  Defaults to 'n_valid_in_column'. Note: This only applies when
  'calc_percent_by = "column"'.

- add_inferential_pvalues:

  Boolean. If `TRUE`, adds the common comparative statistical tests such
  as t-tests, chi-square test, ANOVA, and their nonparametric
  counterparts.

- force_statistical_test:

  List. A list of the format
  `list("variable_name1" ~ "statistical_test1", ...)`, used to specify
  the test for the specific variable. Below are the list of tests. See
  type and enter ?gtsummary::tests for more details, and see
  `tbl_summary() %>% add_p()` section.

  - 't.test': Perform t-test.

  - 'paired.t.test': Perform Paired t-test.

  - 'wilcox.test': Perform Wilcoxon rank-sum test/Mann-Whitney U test.

  - 'paired.wilcox.test': Perform Paired Wilcoxon rank-sum test.

  - 'oneway.test': Perform One-way Analysis of Variance (ANOVA).

  - 'kruskal.test': Perform Kruskal-Wallis test.

  - 'chisq.test': Perform chi-square test of independence.

  - 'chisq.test.no.correct': Perform chi-square test of independence.
    This specifies 'correct = FALSE'.

  - 'fisher.test': Perform Fisher's exact test.

  - 'mcnemar.test': Perform McNemar's test.

  - 'mcnemar.test.wide': Perform McNemar's test.

  - 'prop.test': Perform Test for equality of proportions.

  - 'mood.test': Perform Mood two-sample test of scale.

  - 'lme4': Perform random intercept logistic regression.

  - 'ancova': Perform Analysis of Covariance (ANCOVA).

  - 'emmeans': Perform Estimated Marginal Means or LS-means.

- n_digits_pvalues:

  Numeric. Default is 3. The number of digits of p-values if
  `add_inferential_pvalues` is set to `TRUE`. This can be set to `NULL`
  which tells the code to automatically format the p-values. If the
  p-value \< 0.10, there will be a maximum of three decimal places,
  otherwise, there will only be one.

- bold_significant_pvalues:

  Boolean. If `TRUE` (default), automatically detects if
  `add_inferential_pvalues = TRUE` then bolds the font of p-values under
  0.05 (default) unless set in `bold_significant_pvalues_at`.

- bold_significant_pvalues_at:

  Numeric. The threshold for `bold_significant_pvalues` to bold font if
  significant.

- header:

  String. The header name for the 'Variable' column in the output table.
  Defaults to 'Variable'.

- bold_labels:

  Boolean. If `TRUE` (default), variable labels in the table will be
  bolded.

- italicize_levels:

  Boolean. If `TRUE` (default), categorical variable levels in the table
  will be italicized.

- clean_table:

  Boolean. If `TRUE` (default), cells with a count of 0 are removed
  (replaced with empty space) from the descriptive table, making it
  cleaner.

- table_name:

  String or NULL. The table name. Default is 'auto', where the table
  name depends if the table is descriptive or inferential in nature (has
  p-values), and automatically creates a table name. Set to `NULL` to
  have no table name. Can accept a string that will be used as the table
  name.

- export_to_excel:

  Boolean. If TRUE, exports the results to an Excel file. Defaults to
  FALSE. This is currently not recommended as the output is not
  formatted good.

- export_to_html:

  Boolean. If TRUE, exports the results to an HTML file. Defaults to
  FALSE.

- export_to_png:

  Boolean. If TRUE, exports the results to a PNG file. Defaults to
  FALSE.

- export_to_pdf:

  Boolean. If TRUE, exports the results to a PDF file. Defaults to
  FALSE.

- export_to_word:

  Boolean. If TRUE, exports the results to a Word file. Defaults to
  FALSE.

- export_filename:

  String. The file name of the exported file/s. Multiple exports can be
  done by setting each to TRUE. The default is "Summary Results" and
  will remain as it is if `export_filename = NULL` or
  `export_filename = ""`.

## Value

A `gtsummary` table object (inherits from `gt_tbl`), representing the
simple descriptive statistics.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Load necessary library
library(gtsummary)
library(dplyr)
library(magrittr)

# Create a sample dataframe
set.seed(123)
sample_df <- data.frame(
  Age = sample(c(18:65, NA), 100,
               replace = TRUE, prob = c(rep(1, 48), 5)),
  Gender = sample(c("Male", "Female", NA), 100,
                  replace = TRUE, prob = c(0.48, 0.48, 0.04)),
  Education = sample(c("Elementary", "High School", "Bachelors", "Masters", NA), 100,
                     replace = TRUE, prob = c(0.40, 0.30, 0.15, 0.10, 0.05)),
  Income = sample(c(seq(5000, 40000, by = 5000), NA), 100,
                  replace = TRUE, prob = c(rep(1, 8), 2)),
  Smoker = sample(c(1, 0, NA), 100,
                  replace = TRUE, prob = c(0.25, 0.65, 0.10)),
  Region = sample(c("North", "South", "East", "West", NA), 100,
                  replace = TRUE, prob = c(0.24, 0.24, 0.24, 0.24, .04))
)

# Basic descriptive statistics
# CORRECT: Piped data goes to 'df' argument automatically
sample_df %>% run_summarytable()

# INCORRECT: Causes an error because 'sample_df' is passed to 'summarize_what'
# sample_df %>% run_summarytable(sample_df)

# Descriptive statistics split by Gender, with custom header and labels
sample_df %>%
  run_summarytable(
    split_by = "Gender",
    split_by_header = "Participant Gender",
    rename_variables = list(Age ~ "Age (Years)", 
                            Income ~ "Annual Income", 
                            Smoker ~ "Smoking Status")
  )

# Descriptive statistics with median and IQR for continuous, and count only for categorical
sample_df %>%
  run_summarytable(
    continuous_statistics = "medianIQR",
    categorical_statistics = "n",
    n_digits_continuous = c(0, 0),
    n_digits_categorical = c(0, 0)
  )

# Force 'Smoker' to be continuous (though conceptually it's not, for demonstration)
# And display missing always (with custom text), with bold labels and italicized levels
sample_df %>%
  run_summarytable(
    force_continuous = "Smoker",
    display_missing = "always",
    missing_text = "Data Unavailable",
    missing_stat = "n",
    bold_labels = TRUE,
    italicize_levels = TRUE
  )

# Sort categorical variables by frequency and calculate percentages by row
# (not ideal, but for demonstration only)
sample_df %>%
  run_summarytable(
    sort_categorical_variables_by = "frequency",
    calc_percent_by = "row"
  )

# Perform common comparative statistical tests
sample_df %>%
  run_summarytable(
    split_by = "Gender",
    add_inferential_pvalues = TRUE
  )
} # }
```
