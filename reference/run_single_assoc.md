# Run Single-Variable Association Tests

This function performs univariate statistical tests on categorical
variables to assess the distribution of levels. It automatically applies
[`binom.test`](https://rdrr.io/r/stats/binom.test.html) for binomial
variables (exactly 2 levels) and
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) for multinomial
variables (\> 2 levels).

## Usage

``` r
run_single_assoc(x, y = NULL, subgroup = NULL, force_categorical = NULL, ...)
```

## Arguments

- x:

  A data frame containing the variables to be analyzed.

- y:

  A character vector of column names in `x` to analyze. Defaults to all
  categorical (factor or character) variables.

- subgroup:

  A character vector of column names in `x` to use for stratified
  analysis. If provided, tests are run within each level of each
  subgroup.

- force_categorical:

  A character vector of column names that should be treated as
  categorical even if they are stored as other types (e.g., integers).

- ...:

  Additional parameters passed forward to both
  [`binom.test`](https://rdrr.io/r/stats/binom.test.html) and
  [`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) (e.g., `p` for
  null hypothesis proportions or `conf.level`).

## Value

A list with the following components:

- `summary_table`: A data frame summarizing the test results. Columns
  include: `subgroup`, `group_level`, `feature`, `method`, `estimate`
  (estimated proportion for binomial tests), `lower_95_ci`,
  `upper_95_ci` (confidence interval for binomial tests),
  `test_statistic` (e.g., Chi-squared value for chi-squared tests),
  `p_value`, `note`, and `interpretation`. For tests where a specific
  metric is not applicable, its column will contain `NA`.

- `raw_results`: A named list containing the full test objects returned
  by the underlying R functions.

## Examples

``` r
# Example 1: Binomial test (exactly 2 levels)
# Testing if the distribution of 'am' (Transmission) in mtcars is 50/50
run_single_assoc(mtcars, y = "am", force_categorical = "am")
#> $summary_table
#>   subgroup group_level feature              method estimate lower_95_ci
#> 1   Global         All      am Exact binomial test  0.59375   0.4064492
#>   upper_95_ci test_statistic   p_value                         note
#> 1   0.7630159             NA 0.3770856 Prop: probability of success
#>                interpretation
#> 1 Not significantly different
#> 
#> $raw_results
#> $raw_results$am
#> 
#>  Exact binomial test
#> 
#> data:  counts
#> number of successes = 19, number of trials = 32, p-value = 0.3771
#> alternative hypothesis: true probability of success is not equal to 0.5
#> 95 percent confidence interval:
#>  0.4064492 0.7630159
#> sample estimates:
#> probability of success 
#>                0.59375 
#> 
#> 
#> 

# Example 2: Chi-squared test (more than 2 levels)
# Testing if the distribution of 'Species' in iris is uniform
run_single_assoc(iris, y = "Species")
#> $summary_table
#>   subgroup group_level feature                                   method
#> 1   Global         All Species Chi-squared test for given probabilities
#>   estimate lower_95_ci upper_95_ci test_statistic p_value            note
#> 1       NA          NA          NA              0       1 Goodness of fit
#>                             interpretation
#> 1 Not significantly different distribution
#> 
#> $raw_results
#> $raw_results$Species
#> 
#>  Chi-squared test for given probabilities
#> 
#> data:  counts
#> X-squared = 0, df = 2, p-value = 1
#> 
#> 
#> 

# Example 3: Subgroup analysis
# Testing 'am' distribution stratified by 'vs' (Engine shape)
run_single_assoc(mtcars, y = "am", subgroup = "vs", force_categorical = c("am", "vs"))
#> $summary_table
#>   subgroup group_level feature              method  estimate lower_95_ci
#> 1       vs           0      am Exact binomial test 0.6666667   0.4099252
#> 2       vs           1      am Exact binomial test 0.5000000   0.2303605
#>   upper_95_ci test_statistic   p_value                         note
#> 1   0.8665726             NA 0.2378845 Prop: probability of success
#> 2   0.7696395             NA 1.0000000 Prop: probability of success
#>                interpretation
#> 1 Not significantly different
#> 2 Not significantly different
#> 
#> $raw_results
#> $raw_results$vs_0_am
#> 
#>  Exact binomial test
#> 
#> data:  counts
#> number of successes = 12, number of trials = 18, p-value = 0.2379
#> alternative hypothesis: true probability of success is not equal to 0.5
#> 95 percent confidence interval:
#>  0.4099252 0.8665726
#> sample estimates:
#> probability of success 
#>              0.6666667 
#> 
#> 
#> $raw_results$vs_1_am
#> 
#>  Exact binomial test
#> 
#> data:  counts
#> number of successes = 7, number of trials = 14, p-value = 1
#> alternative hypothesis: true probability of success is not equal to 0.5
#> 95 percent confidence interval:
#>  0.2303605 0.7696395
#> sample estimates:
#> probability of success 
#>                    0.5 
#> 
#> 
#> 
```
