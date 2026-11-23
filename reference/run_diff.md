# Automatic Statistical Comparison with Comprehensive Analysis

Performs an automatic comparison of a numeric outcome variable across
two or more groups. Intelligently selects the appropriate statistical
test (parametric or non-parametric) based on checks for normality (of
data or residuals) and homogeneity of variances, calculates effect
sizes, and performs post-hoc tests when appropriate.

When `outcome` is a character vector of length \> 1 and
`summary_table = TRUE`, an additional `$summary_table` element is
returned alongside the per-outcome list.

When `subgroup` is also specified and `summary_table = TRUE`, additional
named elements `$summary_table_<subgroup_var>` (one per subgroup
variable) are appended; each is a named list of per-level summary
tables, e.g. `$summary_table_Gender$Male`.

## Usage

``` r
run_diff(
  x,
  metadata = NULL,
  outcome,
  group = NULL,
  filter = NULL,
  within = NULL,
  subject_id = NULL,
  paired = FALSE,
  test_type = c("auto", "parametric", "nonparametric"),
  normality_method = c("auto", "shapiro", "lilliefors"),
  alpha_normality = 0.05,
  alpha_variance = 0.05,
  alpha_sphericity = 0.05,
  type_sosquares = 2,
  rm_effect_size = "ges",
  correction = "auto",
  test_alpha = 0.05,
  min_n_threshold = 3,
  calculate_effect_size = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  group_order = NULL,
  subgroup = NULL,
  summary_table = FALSE,
  simplify = TRUE,
  n_digits_pvalues = 3,
  n_digits = 2,
  zero_as_exp = TRUE,
  verbose = FALSE,
  num_cores = 1
)
```

## Arguments

- x:

  A data frame (or an object of class `"run_DIpreprocess"`) containing
  at least the numeric outcome column(s). Non-numeric columns are
  permitted and ignored; only the column(s) named in `outcome` must be
  numeric. When `metadata` is `NULL`, `x` is also used as the metadata
  source (single data-frame workflow).

- metadata:

  A data frame containing sample-level metadata (e.g., grouping,
  subgroups, subject IDs). When `NULL` (default), `x` is used for both
  feature values and metadata - useful when a single data frame holds
  all columns (see Examples).

- outcome:

  A character string or character vector naming the numeric outcome
  variable(s) found in `x`.

- group:

  A character string naming the grouping variable found in `metadata`.
  May be `NULL` for pure within-subjects designs.

- filter:

  Optional character vector of levels in `group` to exclude from the
  analysis. If `x` is a `run_DIpreprocess` object, defaults to the
  identified QC types.

- within:

  Character vector of within-subject factor column name(s) found in
  `metadata`. When supplied, `subject_id` must also be provided.
  Triggers repeated-measures or mixed-ANOVA logic.

- subject_id:

  Character string naming the subject identifier column found in
  `metadata`. Used in two contexts:

  1.  **Repeated-measures / mixed ANOVA** (`within` is non-`NULL`):
      required.

  2.  **Paired two-sample test** (`paired = TRUE`, two groups): optional
      but strongly recommended. When supplied, rows are sorted by
      subject within each group so that `group1[i]` and `group2[i]`
      refer to the same subject, and subject presence is verified across
      both groups. Without `subject_id`, row order is assumed to define
      the pairing (a warning is issued).

- paired:

  Logical; whether to treat observations as paired (two-sample designs
  only). Default `FALSE`.

  - Normality is assessed on the **within-subject differences** rather
    than per group.

  - Homogeneity of variance is not tested (irrelevant for paired
    designs).

  - Provide `subject_id` to guarantee correct pair alignment when rows
    are not already sorted by subject within each group.

- test_type:

  One of `"auto"` (default), `"parametric"`, or `"nonparametric"`.

- normality_method:

  One of `"auto"` (default), `"shapiro"`, or `"lilliefors"`.

- alpha_normality:

  Significance level for normality tests. Default `0.05`.

- alpha_variance:

  Significance level for the variance-homogeneity test. Default `0.05`.

- alpha_sphericity:

  Significance level for Mauchly's test. Default `0.05`.

- type_sosquares:

  Integer (1, 2, or 3) specifying the sums-of-squares type for RM /
  mixed ANOVA. Default `2`.

- rm_effect_size:

  Character vector; one or both of `"ges"` (generalized eta-squared,
  default) and `"pes"` (partial eta-squared). Generalized eta-squared is
  the default because it is directly comparable across designs that
  differ in the number of within- and between-subject factors. Use
  `"pes"` when matching output from software such as SPSS. Specify
  `c("ges", "pes")` to compute both; the first listed metric is treated
  as primary. Applies only to RM / mixed ANOVA designs.

- correction:

  Sphericity correction: `"auto"` (default; applies Greenhouse-Geisser
  when sphericity is violated), `"GG"`, `"HF"`, or `"none"`. Applies
  only to RM / mixed ANOVA.

- test_alpha:

  Significance level for the final test. Default `0.05`.

- min_n_threshold:

  Minimum per-group sample size. Default `3`.

- calculate_effect_size:

  Logical. Default `TRUE`.

- perform_posthoc:

  Logical. Default `TRUE`.

- p_adjust_method:

  P-value adjustment method for multiple testing correction. Applied to
  post-hoc comparison p-values (when available) and to main test
  p-values in summary tables when multiple outcomes are analyzed.
  Default `"BH"`.

- group_order:

  Optional character vector specifying group level order.

- subgroup:

  Optional character vector of column names in `metadata` for stratified
  analysis.

- summary_table:

  Logical; when `outcome` has length \> 1 and `TRUE`, appends
  `$summary_table` (and subgroup summary tables) to the returned list.
  Default `FALSE`.

- simplify:

  Logical; when `summary_table = TRUE`, controls whether
  `$summary_table` (and subgroup tables) are post-processed into a
  compact, presentation-ready form (class `"run_diff_summary"`, with an
  [`as_gt`](https://jllcalorio.github.io/pondeR/reference/as_gt.md)
  method): values are formatted, the effect size and its confidence
  interval are merged, and superscripts on the p-values identify the
  test used. Default `TRUE`.

- n_digits_pvalues:

  Numeric or `NULL`; decimal places for the p-value columns in
  `$summary_table` when `simplify = TRUE`. `NULL` auto-formats (3 d.p.
  when \\p \< 0.10\\, otherwise 1 d.p.), matching
  [`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)'s
  `n_digits_pvalues`. Default `3`.

- n_digits:

  Integer; decimal places for numeric `$summary_table` columns
  (statistic, group means, effect-size estimate/CI) when
  `simplify = TRUE`. Default `2`.

- zero_as_exp:

  Logical; when `TRUE` (default) and `simplify = TRUE`, numeric
  `$summary_table` values that round to `0` at `n_digits` but are
  genuinely non-zero are shown in scientific notation, matching
  [`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)'s
  `zero_as_exp`.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

- num_cores:

  Integer or `"max"`. Default `1`.

## Value

When `outcome` is a single string, a list of class `"run_diff"`
containing: `test_used`, `test_result`, `effect_size`, `posthoc_result`,
`assumptions`, `parametric`, `data_summary`, `outcome`, `group_var`,
`within_vars` (RM only), `subject_id_var` (RM only), `anova_table` (RM
only), `sphericity` (RM only), `sphericity_corrections` (RM only),
`paired`, `test_type`, `test_alpha`, `subgroup_analysis`, `warnings`,
`parameters`, `raw_data`.

When `outcome` is length \> 1, a named list of `"run_diff"` objects,
optionally with `$summary_table` and `$summary_table_<sg>` elements. The
summary table is ready for reports, with columns `Variable`,
`Statistic`, `df`, one `Mean (<group>)` per group, `ES measure`,
`ES (95% CI)`, `Magnitude`, `Significant pairs`, `Pairs tested`,
`p-value`, `Adjusted p-value`, `Significant` and `Interpretation`. With
`simplify = FALSE` it also has `Test`, and `ES`, `ES CI lower` and
`ES CI upper` replace `ES (95% CI)`.

## Details

Automatic Statistical Comparison with Comprehensive Analysis

**Test Family Selection:** If `paired = FALSE`, independent groups are
compared. With 2 groups an independent or paired two-sample test is
chosen; with \> 2 groups an ANOVA-family test is used.

**Assumption Checking (`test_type = "auto"`):**

- **Normality:** Shapiro-Wilk (n \< 50) or Lilliefors (n \>= 50) under
  `"auto"`. Per-group for 2-group comparisons; on model residuals for \>
  2 groups.

- **Homogeneity of Variance:** Levene's test (`car`); falls back to
  Bartlett's test if `car` is unavailable.

**Test Selection Logic (`test_type = "auto"`):**

- Independent two-sample: Mann-Whitney U (normality violated) -\>
  Welch's t (variance violated) -\> Student's t.

- Paired two-sample: Wilcoxon signed-rank (normality violated) -\>
  Paired t.

- Independent \> 2 groups: Kruskal-Wallis (normality violated) -\>
  Welch's ANOVA (variance violated) -\> One-way ANOVA.

**Post-Hoc Tests (\> 2 groups, significant omnibus):** Tukey HSD
(ANOVA), Games-Howell (Welch's ANOVA), Dunn's test (Kruskal-Wallis),
Pairwise t-test (RM / Mixed ANOVA).

**Effect Sizes:** Cohen's d / paired Cohen's d (t-tests), rank-biserial
correlation (Wilcoxon/Mann-Whitney), eta-squared (ANOVA),
epsilon-squared (Kruskal-Wallis), generalized or partial eta-squared (RM
/ Mixed ANOVA).

**Non-Parametric Interpretation (Stochastic Dominance):** For
non-parametric tests the `interpretation` column reads "X tends to have
larger/smaller values than Y".

## Note

**Paired-test workflow:** `paired = TRUE` applies only to two-group
comparisons. When the data are in long format (one row per observation,
a group column, and a subject ID column), pass the data frame as `x`
(with `metadata = NULL` or the same frame) and supply `subject_id`. The
function will sort both groups by subject, verify that every subject
appears in each group, and test normality on the within-subject
differences. Example:

      # sleep: 'extra' is outcome, 'group' is the two-level factor,
      # 'ID' identifies the subject.
      run_diff(sleep, outcome = "extra", group = "group",
               paired = TRUE, subject_id = "ID")

Without `subject_id`, pairs are assumed to be defined by row order
within each group (a warning is emitted). Use `within` + `subject_id`
instead of `paired` for repeated-measures ANOVA (\> 2 time points /
conditions).

## References

Shapiro SS, Wilk MB (1965). An analysis of variance test for normality
(complete samples). *Biometrika*, 52(3-4), 591-611.
[doi:10.1093/biomet/52.3-4.591](https://doi.org/10.1093/biomet/52.3-4.591)

Levene H (1960). Robust tests for equality of variances. In I. Olkin
(Ed.), *Contributions to Probability and Statistics* (pp. 278-292).
Stanford University Press.

Welch BL (1947). The generalization of 'Student's' problem when several
different population variances are involved. *Biometrika*, 34(1-2),
28-35.
[doi:10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)

Welch BL (1951). On the comparison of several mean values: an
alternative approach. *Biometrika*, 38(3-4), 330-336.
[doi:10.1093/biomet/38.3-4.330](https://doi.org/10.1093/biomet/38.3-4.330)

Mann HB, Whitney DR (1947). On a test of whether one of two random
variables is stochastically larger than the other. *The Annals of
Mathematical Statistics*, 18(1), 50-60.
[doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)

Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance
analysis. *Journal of the American Statistical Association*, 47(260),
583-621.
[doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)

Tukey JW (1949). Comparing individual means in the analysis of variance.
*Biometrics*, 5(2), 99-114.
[doi:10.2307/3001913](https://doi.org/10.2307/3001913)

Games PA, Howell JF (1976). Pairwise multiple comparison procedures with
unequal n's and/or variances: a Monte Carlo study. *Journal of
Educational Statistics*, 1(2), 113-125.
[doi:10.3102/10769986001002113](https://doi.org/10.3102/10769986001002113)

Dunn OJ (1964). Multiple comparisons using rank sums. *Technometrics*,
6(3), 241-252.
[doi:10.1080/00401706.1964.10490181](https://doi.org/10.1080/00401706.1964.10490181)

Cohen J (1988). *Statistical Power Analysis for the Behavioral Sciences*
(2nd ed.). Lawrence Erlbaum Associates.

Ben-Shachar MS, Ludecke D, Makowski D (2020). effectsize: Estimation of
effect size indices and standardized parameters. *Journal of Open Source
Software*, 5(56), 2815.
[doi:10.21105/joss.02815](https://doi.org/10.21105/joss.02815)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html),
[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md),
[`get_volcanodata`](https://jllcalorio.github.io/pondeR/reference/get_volcanodata.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# 3+ groups (auto test selection); x serves as its own metadata
res_iris <- run_diff(iris, outcome = "Sepal.Length", group = "Species")
res_iris
#> 
#> === Automatic Statistical Comparison ===
#> 
#> Test Used: Kruskal-Wallis test
#> Parametric: FALSE
#> 
#> Test Results:
#>   Statistic: 96.9374
#>   P-value: < 0.001 ***
#>   Effect size: large (0.651)
#> 
#> Data Summary:
#>       Group  n  Mean        SD Median
#>      setosa 50 5.006 0.3524897    5.0
#>  versicolor 50 5.936 0.5161711    5.9
#>   virginica 50 6.588 0.6358796    6.5
#> 
#> 
#> Post-Hoc Tests (Dunn's test):
#> Significant comparisons:
#>     Group 1    Group 2 Adjusted p-value Significance
#>      setosa versicolor     1.529257e-09          ***
#>      setosa  virginica     6.000296e-22          ***
#>  versicolor  virginica     2.774866e-04          ***
res_iris$posthoc_result
#>      Group 1    Group 2 n (group 1) n (group 2) Statistic      p-value
#> 1     setosa versicolor          50          50  6.106326 1.019504e-09
#> 2     setosa  virginica          50          50  9.741785 2.000099e-22
#> 3 versicolor  virginica          50          50  3.635459 2.774866e-04
#>   Adjusted p-value        Test Significance Significant Mean (group 1)
#> 1     1.529257e-09 Dunn's test          ***        TRUE          5.006
#> 2     6.000296e-22 Dunn's test          ***        TRUE          5.006
#> 3     2.774866e-04 Dunn's test          ***        TRUE          5.936
#>   Mean (group 2)                                         Interpretation
#> 1          5.936    setosa tends to have smaller values than versicolor
#> 2          6.588     setosa tends to have smaller values than virginica
#> 3          6.588 versicolor tends to have smaller values than virginica

# Paired (with subject alignment)
run_diff(sleep, outcome = "extra", group = "group",
         paired = TRUE, subject_id = "ID")
#> 
#> === Automatic Statistical Comparison ===
#> 
#> Test Used: Wilcoxon signed-rank test
#> Parametric: FALSE
#> 
#> Test Results:
#>   Statistic: 0.0000
#>   P-value: 0.0039 **
#>   Effect size: very large (-1.000)
#> 
#> Data Summary:
#>  Group  n Mean       SD Median
#>      1 10 0.75 1.789010   0.35
#>      2 10 2.33 2.002249   1.75

# \donttest{
# Several outcomes in one summary table, with separate x and metadata
outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
res_multi <- run_diff(iris[outcomes], iris["Species"], outcomes, "Species",
                      summary_table = TRUE)
res_multi$summary_table
#>       Variable Statistic     df Mean (setosa) Mean (versicolor)
#> 1 Sepal.Length     96.94      2          5.01              5.94
#> 2  Sepal.Width     49.16 2, 147          3.43              2.77
#> 3 Petal.Length    130.41      2          1.46              4.26
#> 4  Petal.Width    131.19      2          0.25              1.33
#>   Mean (virginica) ES measure    ES (95% CI) Magnitude Significant pairs
#> 1             6.59         ε² 0.65 (0.58, 1)     large                 3
#> 2             2.97         η² 0.40 (0.30, 1)     large                 3
#> 3             5.55         ε² 0.88 (0.86, 1)     large                 3
#> 4             2.03         ε² 0.88 (0.87, 1)     large                 3
#>   Pairs tested p-value Adjusted p-value Significant
#> 1            3  <.001ᵃ           <.001ᵃ        TRUE
#> 2            3  <.001ᵇ           <.001ᵇ        TRUE
#> 3            3  <.001ᵃ           <.001ᵃ        TRUE
#> 4            3  <.001ᵃ           <.001ᵃ        TRUE
#>                                                                                                                                                      Interpretation
#> 1 setosa tends to have smaller values than versicolor;\nsetosa tends to have smaller values than virginica;\nversicolor tends to have smaller values than virginica
#> 2                                         setosa has higher mean than versicolor;\nsetosa has higher mean than virginica;\nversicolor has lower mean than virginica
#> 3 setosa tends to have smaller values than versicolor;\nsetosa tends to have smaller values than virginica;\nversicolor tends to have smaller values than virginica
#> 4 setosa tends to have smaller values than versicolor;\nsetosa tends to have smaller values than virginica;\nversicolor tends to have smaller values than virginica

# Repeated-measures and mixed ANOVA
set.seed(1)
df_long <- data.frame(
  ID        = factor(rep(1:12, each = 3)),
  Time      = factor(rep(c("T1", "T2", "T3"), 12)),
  Group     = rep(c("A", "B"), each = 18),
  HeartRate = rnorm(36, 70, 5) + rep(c(0, 3, 6), 12)
)
run_diff(df_long, outcome = "HeartRate", within = "Time", subject_id = "ID")
#> 
#> === Automatic Statistical Comparison ===
#> 
#> Test Used: One-way Repeated Measures ANOVA
#> Parametric: TRUE
#> 
#> Test Results:
#>   Statistic: 3.8200
#>   P-value: 0.0380 *
#>   Effect size: large (ges = 0.191)
#> 
#> ANOVA Table (all effects):
#>       Effect DFn DFd        SSn     SSd        F        p p<.05   ges
#>  (Intercept)   1  11 193505.741 227.282 9365.294 1.79e-17     * 0.996
#>         Time   2  22    167.893 483.470    3.820 3.80e-02     * 0.191
#> 
#> 
#> Post-Hoc Tests (Pairwise t-test (paired)):
#> No significant pairwise differences found.
run_diff(df_long, outcome = "HeartRate", group = "Group",
         within = "Time", subject_id = "ID")
#> 
#> === Automatic Statistical Comparison ===
#> 
#> Test Used: Two-way Mixed ANOVA
#> Parametric: TRUE
#> 
#> Test Results:
#>   Statistic: 0.3930
#>   P-value: 0.6800
#>   Effect size: small (ges = 0.026)
#> 
#> ANOVA Table (all effects):
#>       Effect DFn DFd        SSn     SSd        F        p p<.05   ges
#>  (Intercept)   1  10 193505.741 222.875 8682.268 4.96e-16     * 0.996
#>        Group   1  10      4.407 222.875    0.198 6.66e-01       0.006
#>         Time   2  20    167.893 465.167    3.609 4.60e-02     * 0.196
#>   Group:Time   2  20     18.303 465.167    0.393 6.80e-01       0.026
#> 
#> 
#> Post-Hoc Tests (Pairwise t-test (paired)):
#> No significant pairwise differences found.
# }
```
