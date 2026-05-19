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
  verbose = TRUE,
  num_cores = 1
)
```

## Arguments

- x:

  A data frame, or an object of class `"run_DIpreprocess"`. If a data
  frame, it should contain only the numeric features/metabolites to be
  analysed.

- metadata:

  A data frame containing sample-level metadata (e.g., grouping,
  subgroups, subject IDs). Required if `x` is a data frame.

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
  `metadata`. Required when `within` is non-`NULL`.

- paired:

  Logical; whether observations are paired. Default `FALSE`.

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

  P-value adjustment method for post-hoc comparisons. Default `"BH"`.

- group_order:

  Optional character vector specifying group level order.

- subgroup:

  Optional character vector of column names in `metadata` for stratified
  analysis.

- summary_table:

  Logical; when `outcome` has length \> 1 and `TRUE`, appends
  `$summary_table` (and subgroup summary tables) to the returned list.
  Default `FALSE`.

- verbose:

  Logical. Default `TRUE`.

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
optionally with `$summary_table` and `$summary_table_<sg>` elements.

## Details

Automatic Statistical Comparison with Comprehensive Analysis

**Test Family Selection:** If `paired = FALSE`, independent groups are
compared. With 2 groups an independent or paired two-sample test is
chosen; with \> 2 groups an ANOVA-family test is used.

**Assumption Checking (`test_type = "auto"`):**

- **Normality:** Shapiro-Wilk (n \< 50) or Lilliefors (n ≥ 50) under
  `"auto"`. Per-group for 2-group comparisons; on model residuals for \>
  2 groups.

- **Homogeneity of Variance:** Levene's test (`car`); falls back to
  Bartlett's test if `car` is unavailable.

**Test Selection Logic (`test_type = "auto"`):**

- Independent two-sample: Mann-Whitney U (normality violated) → Welch's
  t (variance violated) → Student's t.

- Paired two-sample: Wilcoxon signed-rank (normality violated) → Paired
  t.

- Independent \> 2 groups: Kruskal-Wallis (normality violated) → Welch's
  ANOVA (variance violated) → One-way ANOVA.

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

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# --- 3+ Groups (auto) ---
res_iris <- run_diff(iris, "Sepal.Length", "Species")
print(res_iris); summary(res_iris)

# --- Two-sample independent ---
mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
run_diff(mtcars, "mpg", "am")

# --- Paired ---
run_diff(sleep, "extra", "group", paired = TRUE)

# --- Multi-outcome with summary table ---
outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
res_multi <- run_diff(iris, outcomes, "Species",
                      summary_table = TRUE, verbose = FALSE)
res_multi$summary_table

# --- One-way RM ANOVA ---
run_diff(df_long, "HeartRate", within = "Time", subject_id = "ID")

# --- Mixed ANOVA ---
run_diff(df_long, "HeartRate", group = "Group",
         within = "Time", subject_id = "ID")
} # }
```
