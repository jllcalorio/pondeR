# Automatic Statistical Comparison for Categorical Variables

This function performs an automatic comparison between two categorical
variables. It intelligently selects the appropriate statistical test
(Chi-square, Fisher's Exact Test, or McNemar's Test), calculates
measures of association, and performs post-hoc tests when appropriate.

## Usage

``` r
run_assoc(
  x,
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
  verbose = TRUE
)
```

## Arguments

- x:

  A data frame containing the variables for analysis.

- var1:

  A character string specifying the name of the first categorical
  variable (typically rows).

- var2:

  A character string specifying the name of the second categorical
  variable (typically columns).

- weight:

  An optional character string specifying the name of a weight column.
  Use this if your `x` is in a "long" format (e.g., from
  [`as.data.frame.table()`](https://rdrr.io/r/base/table.html)).

- paired:

  A logical indicating whether the observations are paired (for
  McNemar's test). Default is `FALSE`.

- test_type:

  A character string specifying the test strategy for independent
  samples. Must be one of `auto` (default), `chisq`, or `fisher`. This
  parameter is ignored if `paired = TRUE`.

- expected_freq_threshold:

  The minimum expected cell frequency to justify using the Chi-square
  test when `test_type = "auto"`. Default is 5.

- continuity_correction:

  A logical indicating whether to apply Yates' continuity correction for
  2x2 Chi-square tests. Default is `TRUE`.

- fisher_simulate:

  A logical indicating whether to use Monte Carlo simulation for
  Fisher's test p-value in tables larger than 2x2 (which can be
  computationally intensive). Default is `TRUE`.

- test_alpha:

  Significance level for the final statistical test. Default is 0.05.

- calculate_association:

  Logical indicating whether to calculate measures of association (e.g.,
  Cramer's V, Odds Ratio). Default is `TRUE`.

- perform_posthoc:

  Logical indicating whether to perform post-hoc tests for significant
  RxC tables. Default is `TRUE`.

- p_adjust_method:

  P-value adjustment method for post-hoc comparisons. See
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html).
  Default is "BH".

- verbose:

  A logical indicating whether to print detailed messages. Default is
  TRUE.

## Value

An object of class "run_assoc" which is a list containing:

- test_used:

  Character string of the statistical test performed.

- test_result:

  List containing results of the statistical test (a 'htest' object).

- contingency_table:

  The observed contingency table (as a data.frame).

- expected_table:

  The expected frequencies table (as a data.frame, NULL if
  Fisher's/McNemar's).

- residuals:

  Pearson residuals (as a data.frame, NULL if Fisher's/McNemar's).

- stdres:

  Standardized residuals (as a data.frame, NULL if Fisher's/McNemar's).

- association:

  Data frame containing association metrics, magnitude, and
  interpretation.

- posthoc_result:

  Data frame of post-hoc pairwise comparisons (NULL if not applicable).

- assumptions:

  Data frame summarizing assumption check results.

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

The function follows a decision-making process to determine the most
appropriate test:

**Test Family Selection:**

- If `paired = TRUE`, the function assumes the data is for matched pairs
  and requires a 2x2 table. It performs **McNemar's Test**.

- If `paired = FALSE`, the function assumes independent samples and
  builds a contingency table.

**Test Selection Logic (`test_type = "auto"` for independent samples):**
The choice between Chi-square and Fisher's test is based on Cochran's
rule.

- **Assumption Check:** The function first runs a Chi-square test to
  calculate the *expected frequencies* for each cell in the contingency
  table.

- **Decision:**

  - If *any* cell has an expected frequency less than
    `expected_freq_threshold` (default 5), the Chi-square approximation
    may be inaccurate. The function discards this result and performs a
    **Fisher's Exact Test**.

  - If all cells have an expected frequency \>=
    `expected_freq_threshold`, the function proceeds with the
    **Pearson's Chi-square Test**.

**Post-Hoc Tests:** For a significant result on an independent table
larger than 2x2, a post-hoc test is performed using
[`rstatix::pairwise_fisher_test`](https://rpkgs.datanovia.com/rstatix/reference/fisher_test.html)
to determine which specific pairs of categories are driving the
association. P-values are adjusted using the method specified in
`p_adjust_method`.

**Measures of Association (Effect Sizes):** The function automatically
calculates appropriate measures of association:

- **Phi (\\\phi\\)** and **Odds Ratio** for 2x2 independent tables.

- **Cramer's V** for RxC independent tables (larger than 2x2).

- **Odds Ratio (paired)** for 2x2 paired (McNemar's) tables.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Example 1: 2x2 Independent Test (Auto -> Chi-square) ---
# Using 'mtcars' to compare transmission (am) vs. engine shape (vs)
# This will use Chi-square as all expected counts are > 5.
data(mtcars)
res_mtcars <- run_assoc(x = mtcars, var1 = "am", var2 = "vs")
print(res_mtcars)
summary(res_mtcars)

# --- Example 2: RxC Independent Test (Auto -> Chi-square) ---
# Using 'HairEyeColor' data, which is a table.
# We must first convert it to a data frame and use the 'Freq' column as a weight.
data(HairEyeColor)
he_data <- as.data.frame.table(HairEyeColor)
# We'll just compare Hair vs. Eye
he_data_sub <- he_data[he_data$Sex == "Male", ] # Subset to avoid 3D

res_hair <- run_assoc(
  x = he_data_sub,
  var1 = "Hair",
  var2 = "Eye",
  weight = "Freq"
)
print(res_hair)
summary(res_hair) # Will include post-hoc results

# --- Example 3: Paired 2x2 Test (McNemar's) ---
# We must create a dataset for this.
# Imagine a campaign to get people to vote "Yes" (1).
set.seed(42)
N <- 100
before_vote <- sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.7, 0.3))
after_vote <- before_vote
# 20% of "No" voters switch to "Yes"
after_vote[before_vote == "No"][1:14] <- "Yes"
# 10% of "Yes" voters switch to "No"
after_vote[before_vote == "Yes"][1:3] <- "No"

paired_data <- data.frame(
  person_id = 1:N,
  before = before_vote,
  after = after_vote
)

res_paired <- run_assoc(
  x = paired_data,
  var1 = "before",
  var2 = "after",
  paired = TRUE
)
summary(res_paired)

# --- Example 4: Force a test_type ---
# Force Chi-square on the mtcars data, with continuity correction
res_mtcars_chisq <- run_assoc(
  x = mtcars,
  var1 = "am",
  var2 = "vs",
  test_type = "chisq",
  continuity_correction = TRUE
)
print(res_mtcars_chisq)
} # }
```
