# Automatic Statistical Comparison for Categorical Variables

Performs an automatic comparison between two categorical variables. It
intelligently selects the appropriate statistical test (Chi-square,
Fisher's Exact Test, or McNemar's Test), calculates measures of
association, and performs post-hoc tests when appropriate.

## Usage

``` r
run_assoc(
  x,
  var1,
  var2,
  weight = NULL,
  paired = FALSE,
  force_categorical = FALSE,
  test_type = c("auto", "chisq", "fisher"),
  expected_freq_threshold = 5,
  continuity_correction = TRUE,
  fisher_simulate = TRUE,
  test_alpha = 0.05,
  calculate_association = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  verbose = TRUE
)
```

## Arguments

- x:

  A data frame containing the variables for analysis.

- var1:

  A character vector of one or more column names in `x` specifying the
  categorical variable(s) to test. Each variable in `var1` will be
  independently paired with `var2` for a test of association. All listed
  columns must be categorical (factor or character), unless
  `force_categorical = TRUE`.

- var2:

  A character string specifying the name of a single categorical
  variable (typically the grouping or outcome column). Each variable in
  `var1` will be compared against this variable.

- weight:

  An optional character string specifying the name of a weight column.
  Use this if `x` is in a frequency-aggregated format (e.g., from
  [`as.data.frame.table()`](https://rdrr.io/r/base/table.html)).

- paired:

  A logical indicating whether the observations are paired (for
  McNemar's test). Default is `FALSE`.

- force_categorical:

  Logical. If `TRUE`, any column in `var1` or `var2` that is not already
  a factor or character will be coerced to a factor before analysis.
  This is useful when numeric codes represent categorical groups (e.g.,
  `0`/`1`, `1`/`2`/`3`). A warning is issued for each coerced column.
  Default is `FALSE`.

- test_type:

  A character string specifying the test strategy for independent
  samples. One of `"auto"` (default), `"chisq"`, or `"fisher"`. Ignored
  if `paired = TRUE`.

- expected_freq_threshold:

  The minimum expected cell frequency to justify using the Chi-square
  test when `test_type = "auto"`. Default is 5.

- continuity_correction:

  A logical indicating whether to apply Yates' continuity correction for
  2x2 Chi-square or McNemar's tests. Default is `TRUE`.

- fisher_simulate:

  A logical indicating whether to use Monte Carlo simulation for
  Fisher's Exact Test p-values. Recommended for tables larger than 2x2.
  Default is `TRUE`.

- test_alpha:

  Significance level for the final statistical test. Default is 0.05.

- calculate_association:

  Logical. Whether to calculate measures of association (Cramer's V,
  Phi, Odds Ratio). Requires `effectsize`. Default is `TRUE`.

- perform_posthoc:

  Logical. Whether to perform post-hoc pairwise Fisher's tests for
  significant tables with one binary dimension. Requires `rstatix`.
  Default is `TRUE`.

- p_adjust_method:

  P-value adjustment method for post-hoc comparisons. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default is
  `"BH"`.

- verbose:

  A logical indicating whether to print detailed messages. Default is
  `TRUE`.

## Value

A list with class `"run_assoc"` containing:

- test_used:

  Character string of the statistical test performed.

- test_result:

  An `htest` object from the selected test.

- contingency_table:

  Observed contingency table as a data frame.

- expected_table:

  Expected frequencies as a data frame (`NULL` for Fisher's/McNemar's).

- residuals:

  Pearson residuals as a data frame (`NULL` for Fisher's/McNemar's).

- stdres:

  Standardized residuals as a data frame (`NULL` for
  Fisher's/McNemar's). Values beyond ±1.96 or ±2.58 indicate cells
  contributing to significance.

- association:

  Data frame of association metrics (metric, estimate, CI, magnitude,
  interpretation), or `NULL`.

- posthoc_result:

  Data frame of post-hoc pairwise comparisons, or `NULL`.

- assumptions:

  Data frame summarising assumption check results.

- var1:

  Name of the row variable.

- var2:

  Name of the column variable.

- paired:

  The paired setting used.

- test_alpha:

  The significance level used.

- warnings:

  Character vector of any warnings generated.

- parameters:

  List of all parameters used in the analysis.

## Details

Automatic Statistical Comparison for Categorical Variables

**Test Family Selection:**

- If `paired = TRUE`, requires a 2x2 table and performs **McNemar's
  Test**.

- If `paired = FALSE`, assumes independent samples and builds a
  contingency table.

**Test Selection Logic (`test_type = "auto"` for independent samples):**
The choice between Chi-square and Fisher's test is based on Cochran's
rule.

- The function first calculates expected frequencies via a Chi-square
  check.

- If *any* cell has an expected frequency below
  `expected_freq_threshold` (default 5), Fisher's Exact Test is used.

- Otherwise, Pearson's Chi-square Test is used.

**Post-Hoc Tests:** For a significant result on a table where one
dimension is binary (2 x C or R x 2), pairwise Fisher's tests are
performed via
[`rstatix::pairwise_fisher_test`](https://rpkgs.datanovia.com/rstatix/reference/fisher_test.html).
P-values are adjusted using `p_adjust_method`.

**Measures of Association (Effect Sizes):**

- **Phi** and **Odds Ratio** for 2x2 independent tables.

- **Cramer's V** for RxC independent tables (larger than 2x2).

- **Odds Ratio (paired)** for 2x2 paired (McNemar's) tables.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Example 1: 2x2 Independent Test (Auto -> Chi-square) ---
data(mtcars)
res <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
                 force_categorical = TRUE)
print(res)
summary(res)

# --- Example 2: RxC Independent Test ---
data(HairEyeColor)
he <- as.data.frame.table(HairEyeColor)
he_male <- he[he$Sex == "Male", ]
res_hair <- run_assoc(x = he_male, var1 = "Hair", var2 = "Eye", weight = "Freq")
summary(res_hair)

# --- Example 3: Paired 2x2 (McNemar's) ---
set.seed(42)
N <- 100
before <- sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.7, 0.3))
after  <- before
after[before == "No"][1:14]  <- "Yes"
after[before == "Yes"][1:3]  <- "No"
paired_df <- data.frame(before = before, after = after)
res_paired <- run_assoc(x = paired_df, var1 = "before", var2 = "after",
                         paired = TRUE)
summary(res_paired)

# --- Example 4: Force Chi-square ---
res_chisq <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
                        test_type = "chisq", continuity_correction = TRUE,
                        force_categorical = TRUE)
print(res_chisq)
} # }
```
