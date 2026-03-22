# Automatic Statistical Comparison with Comprehensive Analysis

This function performs an automatic comparison of a numeric outcome
variable across two or more groups. It intelligently selects the
appropriate statistical test (parametric or non-parametric) based on
checks for normality (of data or residuals) and homogeneity of
variances, calculates effect sizes, and performs post-hoc tests when
appropriate.

## Usage

``` r
run_diff(
  x,
  outcome,
  group,
  paired = FALSE,
  test_type = c("auto", "parametric", "nonparametric"),
  normality_method = c("auto", "shapiro", "lilliefors"),
  alpha_normality = 0.05,
  alpha_variance = 0.05,
  test_alpha = 0.05,
  min_n_threshold = 3,
  calculate_effect_size = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  group_order = NULL,
  subgroup = NULL,
  verbose = TRUE,
  num_cores = 1
)
```

## Arguments

- x:

  A data frame containing the variables for analysis.

- outcome:

  A character string or character vector specifying the name(s) of the
  numeric outcome variable(s). When a vector is supplied, the function
  loops over each outcome and returns a named list of `"run_diff"`
  objects instead of a single object.

- group:

  A character string specifying the name of the grouping factor
  variable.

- paired:

  A logical indicating whether the observations are paired. Default is
  FALSE.

- test_type:

  A character string specifying the test strategy. Must be one of
  `"auto"` (default), `"parametric"`, or `"nonparametric"`.

- normality_method:

  Method for assessing normality: "shapiro" (Shapiro-Wilk test),
  "lilliefors" (Lilliefors (Kolmogorov-Smirnov) test), or "auto" (uses
  Shapiro-Wilk for n \< 50, Lilliefors for n \>= 50). Default is "auto".

- alpha_normality:

  Significance level for normality tests. Default is 0.05.

- alpha_variance:

  Significance level for variance homogeneity tests. Default is 0.05.

- test_alpha:

  Significance level for the final statistical test. Default is 0.05.

- min_n_threshold:

  Minimum sample size required per group. Default is 3.

- calculate_effect_size:

  Logical indicating whether to calculate effect sizes. Default is TRUE.

- perform_posthoc:

  Logical indicating whether to perform post-hoc tests for \>2 groups.
  Default is TRUE.

- p_adjust_method:

  The p-value adjustment method for post-hoc comparisons: "holm",
  "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default
  is "BH".

- group_order:

  Character vector specifying the order of groups. If NULL, uses factor
  levels in alphabetical order.

- subgroup:

  Character vector or NULL. Column name(s) in `x` to use for subgroup
  analysis. Each must be a categorical column (factor or character) in
  `x`. When specified, the full analysis is repeated independently for
  each level of each subgroup column. Results are returned under
  `$subgroup_analysis` in the output, nested as
  `subgroup_var -> level -> run_diff result`. Default is NULL.

- verbose:

  A logical indicating whether to print detailed messages. Default is
  TRUE.

- num_cores:

  Integer specifying the number of cores to use for parallel processing,
  or the character string "max". If "max", the function uses the number
  of available cores minus 2 (to reserve system resources), with a
  minimum of 1. Default is 1. Supports both Windows and POSIX
  (Mac/Linux) systems.

## Value

When `outcome` is a single string, returns an object of class
`"run_diff"` (a list) containing:

- test_used:

  Character string of the statistical test performed.

- test_result:

  List containing results of the statistical test.

- effect_size:

  List containing effect size estimate, confidence interval, magnitude,
  and interpretation.

- posthoc_result:

  Data frame of post-hoc pairwise comparisons (NULL if not applicable).

- assumptions:

  Data frame summarizing assumption check results.

- data_summary:

  Data frame with comprehensive descriptive statistics per group.

- outcome:

  Name of the outcome variable.

- group_var:

  Name of the group variable.

- paired:

  The paired setting used.

- test_type:

  The test_type strategy used.

- test_alpha:

  The significance level used.

- subgroup_analysis:

  Named list of subgroup results, nested as
  `subgroup_var -> level -> run_diff result`. NULL if no subgroup
  specified.

- warnings:

  Character vector of any warnings generated.

- parameters:

  List of all parameters used in the analysis.

- raw_data:

  Data frame containing the variables used.

When `outcome` is a character vector of length \> 1, returns a named
list where each element is a `"run_diff"` object for the corresponding
outcome variable.

## Details

The function follows a decision-making process to determine the most
appropriate test:

**Test Family Selection:** If `paired = FALSE`, it compares independent
groups. If there are 2 groups, it selects a two-sample test family. If
there are more than 2 groups, it selects an ANOVA-family test.

**Assumption Checking (only for `test_type = "auto"`):**

- **Normality:** Assessed using the Shapiro-Wilk test for small samples
  (n \< 50). For larger samples (n \>= 50), uses the Lilliefors
  (Kolmogorov-Smirnov) test from the `nortest` package. For two groups,
  each group is tested. For more than two groups, the normality of model
  residuals is tested.

- **Homogeneity of Variance:** For independent tests, assessed using
  Levene's test from the `car` package. If `car` is not installed, falls
  back to Bartlett's test.

**Test Selection Logic (`test_type = "auto"`):**

- **Independent Two-Sample:**

  - If normality is violated -\> Mann-Whitney U test (Wilcoxon
    rank-sum).

  - If normality is OK but variance is violated -\> Welch's t-test.

  - If both are OK -\> Student's t-test.

- **Paired Two-Sample:**

  - If normality is violated -\> Wilcoxon signed-rank test.

  - If normality is OK -\> Paired t-test.

- **Independent ANOVA (\>2 groups):**

  - If normality is violated -\> Kruskal-Wallis test.

  - If normality is OK but variance violated -\> Welch's ANOVA.

  - If both OK -\> One-way ANOVA.

**Post-Hoc Tests:** For significant results with \>2 groups, appropriate
post-hoc tests are automatically performed:

- Tukey HSD for one-way ANOVA

- Games-Howell for Welch's ANOVA

- Dunn's test for Kruskal-Wallis

**Effect Sizes:** The function automatically calculates appropriate
effect sizes:

- Cohen's d for t-tests

- Rank-biserial correlation for Mann-Whitney U

- Eta-squared for ANOVA

- Epsilon-squared for Kruskal-Wallis

## References

Aaron Schlege, *Games-Howell Post-Hoc Test*.
<https://rpubs.com/aaronsc32/games-howell-test>.

Albers, C., & Lakens, D. (2018). *When power analyses based on pilot
data are biased: Inaccurate effect size estimators and follow-up bias*.
Journal of Experimental Social Psychology. URL
<https://doi.org/10.1016/j.jesp.2017.09.004>

Algina, J., Keselman, H. J., & Penfield, R. D. (2006). *Confidence
intervals for an effect size when variances are not equal*. Journal of
Modern Applied Statistical Methods, 5(1), 2. URL
<https://jmasm.com/index.php/jmasm/article/view/222> FULL PDF
<https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1256&context=jmasm>

Allen, R. (2017). *Statistics and Experimental Design for Psychologists:
A Model Comparison Approach*. World Scientific Publishing Company.

Bache S, Wickham H (2022). **magrittr: A Forward-Pipe Operator for R**.
doi:10.32614/CRAN.package.magrittr
<https://doi.org/10.32614/CRAN.package.magrittr>, R package version
2.0.3, <https://CRAN.R-project.org/package=magrittr>.

Bartlett, M. S. (1937). *Properties of sufficiency and statistical
tests*. Proceedings of the Royal Society of London Series A 160,
268–282. doi:10.1098/rspa.1937.0109.

Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) *The New S
Language*. Wadsworth & Brooks/Cole.

Benjamini, Y., and Hochberg, Y. (1995). *Controlling the false discovery
rate: a practical and powerful approach to multiple testing*. Journal of
the Royal Statistical Society Series B, 57, 289–300.
doi:10.1111/j.2517-6161.1995.tb02031.x.

Benjamini, Y., and Yekutieli, D. (2001). *The control of the false
discovery rate in multiple testing under dependency*. Annals of
Statistics, 29, 1165–1188. doi:10.1214/aos/1013699998.

Ben-Shachar M, Lüdecke D, Makowski D (2020). *effectsize: Estimation of
Effect Size Indices and Standardized Parameters*. Journal of Open Source
Software, 5(56), 2815. doi: 10.21105/joss.02815

Chambers, J. M., Freeny, A and Heiberger, R. M. (1992) *Analysis of
variance; designed experiments*. Chapter 5 of Statistical Models in S
eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

Chambers, J. M. and Hastie, T. J. (1992) *Statistical Models in S*.
Wadsworth & Brooks/Cole.

Cliff, N. (1993). *Dominance statistics: Ordinal analyses to answer
ordinal questions*. Psychological bulletin, 114(3), 494.

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences (2nd Ed.)*. New York: Routledge.

Cureton, E. E. (1956). *Rank-biserial correlation*. Psychometrika,
21(3), 287-290.

Dallal, G.E. and Wilkinson, L. (1986): *An analytic approximation to the
distribution of Lilliefors' test for normality*. The American
Statistician, 40, 294–296.

David F. Bauer (1972). *Constructing confidence sets using rank
statistics*. Journal of the American Statistical Association 67,
687–690. doi:10.1080/01621459.1972.10481279.

Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021, May 7).
*Why Hedges’ g\*s based on the non-pooled standard deviation should be
reported with Welch's t-test*. doi:10.31234/osf.io/tu6mp

Dunn, O. J. (1964) *Multiple comparisons using rank sums Technometrics*,
6(3):241-252.

Eric Langford (2006) *Quartiles in Elementary Statistics*, Journal of
Statistics Education 14 3; doi:10.1080/10691898.2006.11910589

Fox, J. (2016) *Applied Regression Analysis and Generalized Linear
Models*, Third Edition. Sage.

Fox J, Weisberg S (2019). **An R Companion to Applied Regression**,
Third edition. Sage, Thousand Oaks CA.
<https://www.john-fox.ca/Companion/>.

Glass, G. V. (1965). *A ranking variable analogue of biserial
correlation: Implications for short-cut item analysis*. Journal of
Educational Measurement, 2(1), 91-95.

Gross J, Ligges U (2015). **nortest: Tests for Normality**.
doi:10.32614/CRAN.package.nortest
<https://doi.org/10.32614/CRAN.package.nortest>, R package version
1.0-4, <https://CRAN.R-project.org/package=nortest>.

Hedges, L. V. & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Orlando, FL: Academic Press.

Hochberg, Y. (1988). *A sharper Bonferroni procedure for multiple tests
of significance*. Biometrika, 75, 800–803. doi:10.2307/2336325

Holm, S. (1979). *A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics*, 6, 65–70.
<https://www.jstor.org/stable/4615733>.

Hommel, G. (1988). *A stagewise rejective multiple test procedure based
on a modified Bonferroni test*. Biometrika, 75, 383–386.
doi:10.2307/2336190.

Hunter, J. E., & Schmidt, F. L. (2004). *Methods of meta-analysis:
Correcting error and bias in research findings*. Sage.

Hyndman, R. J. and Fan, Y. (1996) *Sample quantiles in statistical
packages*, American Statistician 50, 361–365. doi:10.2307/2684934.

Kassambara A (2023). **rstatix: Pipe-Friendly Framework for Basic
Statistical Tests**. doi:10.32614/CRAN.package.rstatix
<https://doi.org/10.32614/CRAN.package.rstatix>, R package version
0.7.2, <https://CRAN.R-project.org/package=rstatix>.

Kelley, T. (1935) *An unbiased correlation ratio measure*. Proceedings
of the National Academy of Sciences. 21(9). 554-559.

Kerby, D. S. (2014). *The simple difference formula: An approach to
teaching nonparametric correlation*. Comprehensive Psychology, 3, 11-IT.

King, B. M., & Minium, E. W. (2008). *Statistical reasoning in the
behavioral sciences*. John Wiley & Sons Inc.

Komsta L, Novomestky F (2022). **moments: Moments, Cumulants, Skewness,
Kurtosis and Related Tests**. doi:10.32614/CRAN.package.moments
<https://doi.org/10.32614/CRAN.package.moments>, R package version
0.14.1, <https://CRAN.R-project.org/package=moments>.

Makkonen, L. and Pajari, M. (2014) *Defining Sample Quantiles by the
True Rank Probability*, Journal of Probability and Statistics; Hindawi
Publ.Corp. doi:10.1155/2014/326579

Myles Hollander and Douglas A. Wolfe (11973). *Nonparametric Statistical
Methods*. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75
(two-sample). Or second edition (1999).

Myles Hollander and Douglas A. Wolfe (1973), *Nonparametric Statistical
Methods*. New York: John Wiley & Sons. Pages 115–120.

Müller K, Wickham H (2025). **tibble: Simple Data Frames**.
doi:10.32614/CRAN.package.tibble
<https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0,
<https://CRAN.R-project.org/package=tibble>.

Olejnik, S., & Algina, J. (2003). *Generalized eta and omega squared
statistics: measures of effect size for some common research designs*.
Psychological methods, 8(4), 434.

Patrick Royston (1982). *Algorithm AS 181: The W test for Normality*.
Applied Statistics, 31, 176–180. doi:10.2307/2347986.

Patrick Royston (1982). *An extension of Shapiro and Wilk's W test for
normality to large samples*. Applied Statistics, 31, 115–124.
doi:10.2307/2347973.

Patrick Royston (1995). *Remark AS R94: A remark on Algorithm AS 181:
The W test for normality*. Applied Statistics, 44, 547–551.
doi:10.2307/2986146.

Ruxton, G.D., and Beauchamp, G. (2008) *‘Time for some a priori thinking
about post hoc testing’*, Behavioral Ecology, 19(3), pp. 690-693. doi:
10.1093/beheco/arn020. In-text citations: (Ruxton and Beauchamp, 2008)

Sangseok Lee, Dong Kyu Lee. *What is the proper way to apply the
multiple comparison test?*. Korean J Anesthesiol. 2018;71(5):353-360.

Sarkar, S. (1998). *Some probability inequalities for ordered MTP2
random variables: a proof of Simes conjecture*. Annals of Statistics,
26, 494–504. doi:10.1214/aos/1028144846.

Sarkar, S., and Chang, C. K. (1997). *The Simes method for multiple
hypothesis testing with positively dependent test statistics*. Journal
of the American Statistical Association, 92, 1601–1608.
doi:10.2307/2965431.

Shaffer, J. P. (1995). *Multiple hypothesis testing. Annual Review of
Psychology*, 46, 561–584. doi:10.1146/annurev.ps.46.020195.003021. (An
excellent review of the area.)

Steiger, J. H. (2004). *Beyond the F test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis*. Psychological Methods, 9, 164-182.

Stephens, M.A. (1974): *EDF statistics for goodness of fit and some
comparisons*. Journal of the American Statistical Association, 69,
730–737.

Thode Jr., H.C. (2002): *Testing for Normality*. Marcel Dekker, New
York.

Tomczak, M., & Tomczak, E. (2014). *The need to report effect size
estimates revisited*. An overview of some recommended measures of effect
size.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E,
Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K,
Vaughan D, Wilke C, Woo K, Yutani H (2019). *“Welcome to the
tidyverse.”* *Journal of Open Source Software*, *4*(43), 1686.
doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

Wicklin, R. (2017) *Sample quantiles: A comparison of 9 definitions*;
SAS Blog.
https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.htm

Wright, S. P. (1992). *Adjusted P-values for simultaneous inference*.
Biometrics, 48, 1005–1013. doi:10.2307/2532694. (Explains the adjusted
P-value approach.)

R Core Team (2025). **R: A Language and Environment for Statistical
Computing**. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Basic Usage with 3+ Groups ---
# Using the classic 'iris' dataset. Sepal.Length is numeric, Species is a factor.
# This will likely trigger the One-way ANOVA -> Tukey HSD pathway.
data(iris)
res_iris <- run_diff(x = iris, outcome = "Sepal.Length", group = "Species")
print(res_iris)
summary(res_iris)

# --- Two-Sample Independent Comparison ---
# Using 'mtcars' to compare miles per gallon (mpg) between automatic (0) and manual (1) cars.
# The variances are unequal, so this will trigger the Welch's t-test.
data(mtcars)
mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
res_mtcars <- run_diff(x = mtcars, outcome = "mpg", group = "am")
print(res_mtcars)

# --- Paired-Sample Comparison ---
# Using the 'sleep' dataset, which is designed for a paired t-test.
# It compares the extra hours of sleep for the same subjects under two drugs.
data(sleep)
# We specify the subject ID is implicit in the row order for a paired test.
res_sleep <- run_diff(x = sleep, outcome = "extra", group = "group", paired = TRUE)
summary(res_sleep)

# --- Forcing a Non-Parametric Test ---
# Use 'iris' data again, but force the non-parametric Kruskal-Wallis test.
data(iris)
res_iris_np <- run_diff(
  x = iris,
  outcome = "Sepal.Length",
  group = "Species",
  test_type = "nonparametric"
)
print(res_iris_np)

# --- Forcing a Parametric Test & Custom P-Value Adjustment ---
# Use the 'ToothGrowth' data with 'dose' as the grouping factor.
# We force an ANOVA (even if normality is slightly off) and use Bonferroni correction for post-hocs.
data(ToothGrowth)
# Convert dose to a factor to be treated as a group
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
res_tooth <- run_diff(
  x = ToothGrowth,
  outcome = "len",
  group = "dose",
  test_type = "parametric",
  p_adjust_method = "bonferroni"
)
# Look at the post-hoc results and the adjusted p-values
print(res_tooth$posthoc_result)

# --- Demonstrating Custom Group Ordering ---
# Use 'PlantGrowth' and specify a non-alphabetical order for the groups.
data(PlantGrowth)
custom_order <- c("trt1", "ctrl", "trt2")
res_plant <- run_diff(
  x = PlantGrowth,
  outcome = "weight",
  group = "group",
  group_order = custom_order,
  verbose = FALSE # Suppress messages for cleaner output
)
# The summary and plots will respect the custom order
print(res_plant$data_summary)

# --- Disabling Features for Quicker Analysis ---
# If you only need the main test result, you can disable effect sizes and post-hoc tests.
# Using 'mtcars' again with the 'cyl' variable as a group.
data(mtcars)
mtcars$cyl <- as.factor(mtcars$cyl)
res_mtcars_fast <- run_diff(
  x = mtcars,
  outcome = "mpg",
  group = "cyl",
  calculate_effect_size = FALSE,
  perform_posthoc = FALSE
)
print(res_mtcars_fast)
} # }
```
