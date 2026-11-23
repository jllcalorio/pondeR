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
  fisher_simulate = FALSE,
  test_alpha = 0.05,
  calculate_association = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  verbose = FALSE
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
  Default is `FALSE`.

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
  `FALSE`.

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
  Fisher's/McNemar's). Values beyond \u00b1 1.96 or \u00b1 2.58 indicate
  cells contributing to significance.

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

## References

Pearson K (1900). On the criterion that a given system of deviations
from the probable in the case of a correlated system of variables is
such that it can be reasonably supposed to have arisen from random
sampling. *Philosophical Magazine*, 50(302), 157-175.
[doi:10.1080/14786440009463897](https://doi.org/10.1080/14786440009463897)

Fisher RA (1922). On the interpretation of chi-squared from contingency
tables, and the calculation of P. *Journal of the Royal Statistical
Society*, 85(1), 87-94.
[doi:10.2307/2340521](https://doi.org/10.2307/2340521)

McNemar Q (1947). Note on the sampling error of the difference between
correlated proportions or percentages. *Psychometrika*, 12(2), 153-157.
[doi:10.1007/BF02295996](https://doi.org/10.1007/BF02295996)

Cochran WG (1954). Some methods for strengthening the common chi-squared
tests. *Biometrics*, 10(4), 417-451.
[doi:10.2307/3001616](https://doi.org/10.2307/3001616)

Cramer H (1946). *Mathematical Methods of Statistics*. Princeton
University Press.

Ben-Shachar MS, Ludecke D, Makowski D (2020). effectsize: Estimation of
effect size indices and standardized parameters. *Journal of Open Source
Software*, 5(56), 2815.
[doi:10.21105/joss.02815](https://doi.org/10.21105/joss.02815)

## See also

[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# --- Example 1: 2x2 Independent Test (Auto -> Chi-square) ---
data(mtcars)
res <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
                 force_categorical = TRUE)
print(res)
#> 
#> === Automatic Categorical Comparison ===
#> 
#> Test Used: Pearson's Chi-square test
#> Variables: am vs. vs
#> Paired: FALSE
#> 
#> Test Results:
#>   X-squared: 0.3475
#>   P-value: 0.5555
#>   Interpretation: No significant association detected.
#>   NA
#> 
#> Observed Frequencies:
#>    0 1
#> 0 12 7
#> 1  6 7
#> 
#> Warnings:
#>   - Column 'am' (numeric) coerced to factor via `force_categorical = TRUE`.
#>   - Column 'vs' (numeric) coerced to factor via `force_categorical = TRUE`.
summary(res)
#> 
#> === Detailed Categorical Comparison Summary ===
#> 
#> Row Variable:    am
#> Col Variable:    vs
#> Paired:          FALSE
#> Test Strategy:   auto (alpha = 0.05)
#> Test Selected:   Pearson's Chi-square test
#> 
#> --- Observed Frequencies ---
#>    0 1
#> 0 12 7
#> 1  6 7
#> 
#> --- Expected Frequencies ---
#>       0    1
#> 0 10.69 8.31
#> 1  7.31 5.69
#> 
#> --- Pearson Residuals ---
#>       0     1
#> 0  0.40 -0.46
#> 1 -0.49  0.55
#> 
#> --- Standardized Residuals ---
#>       0     1
#> 0  0.95 -0.95
#> 1 -0.95  0.95
#>   (Values beyond |1.96| or |2.58| indicate cells driving the association)
#> 
#> --- Assumption Checks ---
#>                 Check       Result Decision
#>  Expected freq. > 5.0 4/4 cells OK       OK
#> 
#> --- Statistical Test Results ---
#>   Method: Pearson's Chi-squared test with Yates' continuity correction
#>   X-squared: 0.3475
#>   df: 1
#>   p-value: 0.5555
#> 
#> --- Measures of Association ---
#> Metric:    Phi / Odds Ratio
#> Estimate:  0.000 [95% CI: 0.000, 1.000]
#> Magnitude: tiny
#> 
#> --- Interpretation ---
#> Result: NO significant association detected (p >= 0.05).
#> 
#> --- Warnings ---
#>   - Column 'am' (numeric) coerced to factor via `force_categorical = TRUE`.
#>   - Column 'vs' (numeric) coerced to factor via `force_categorical = TRUE`.

# --- Example 2: RxC Independent Test ---
data(HairEyeColor)
he <- as.data.frame.table(HairEyeColor)
he_male <- he[he$Sex == "Male", ]
res_hair <- run_assoc(x = he_male, var1 = "Hair", var2 = "Eye", weight = "Freq")
summary(res_hair)
#> 
#> === Detailed Categorical Comparison Summary ===
#> 
#> Row Variable:    Hair
#> Col Variable:    Eye
#> Paired:          FALSE
#> Test Strategy:   auto (alpha = 0.05)
#> Test Selected:   Fisher's Exact Test
#> 
#> --- Observed Frequencies ---
#>       Brown Blue Hazel Green
#> Black    32   11    10     3
#> Brown    53   50    25    15
#> Red      10   10     7     7
#> Blond     3   30     5     8
#> 
#> --- Assumption Checks ---
#>           Check                   Result Decision
#>  Test Selection Fisher's Exact Test used       OK
#> 
#> --- Statistical Test Results ---
#>   Method: Fisher's Exact Test for Count Data with simulated p-value
#>   (based on 2000 replicates)
#>   p-value: < 0.001
#> 
#> --- Measures of Association ---
#> Metric:    Cramer's V
#> Estimate:  NA [95% CI: 0.092, 1.000]
#> Magnitude: unknown
#> 
#> --- Interpretation ---
#> Result: SIGNIFICANT association detected (p < 0.05).
#> Strength: unknown (Cramer's V = NA).
#> 
#> --- Warnings ---
#>   - Exact Fisher's test failed (possibly out of memory). Switching to simulation.

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
#> 
#> === Detailed Categorical Comparison Summary ===
#> 
#> Row Variable:    before
#> Col Variable:    after
#> Paired:          TRUE
#> Test Strategy:   mcnemar (alpha = 0.05)
#> Test Selected:   McNemar's test
#> 
#> --- Observed Frequencies ---
#>     No Yes
#> No  52  14
#> Yes  3  31
#> 
#> --- Assumption Checks ---
#>           Check           Result Decision
#>  Data Structure 2x2 paired table       OK
#> 
#> --- Statistical Test Results ---
#>   Method: McNemar's Chi-squared test with continuity correction
#>   McNemar's chi-squared: 5.8824
#>   df: 1
#>   p-value: 0.0153
#> 
#> --- Measures of Association ---
#> Metric:    Odds Ratio (paired)
#> Estimate:  4.667 [95% CI: 1.341, 16.239]
#> Magnitude: positive change
#> 
#> --- Interpretation ---
#> Result: SIGNIFICANT association detected (p < 0.05).
#> Strength: positive change (Odds Ratio (paired) = 4.667).

# --- Example 4: Force Chi-square ---
res_chisq <- run_assoc(x = mtcars, var1 = "am", var2 = "vs",
                        test_type = "chisq", continuity_correction = TRUE,
                        force_categorical = TRUE)
print(res_chisq)
#> 
#> === Automatic Categorical Comparison ===
#> 
#> Test Used: Pearson's Chi-square test (forced)
#> Variables: am vs. vs
#> Paired: FALSE
#> 
#> Test Results:
#>   X-squared: 0.3475
#>   P-value: 0.5555
#>   Interpretation: No significant association detected.
#>   NA
#> 
#> Observed Frequencies:
#>    0 1
#> 0 12 7
#> 1  6 7
#> 
#> Warnings:
#>   - Column 'am' (numeric) coerced to factor via `force_categorical = TRUE`.
#>   - Column 'vs' (numeric) coerced to factor via `force_categorical = TRUE`.
```
