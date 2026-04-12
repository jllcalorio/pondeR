# Get Features with Non-Zero Coefficients from a Regularized Regression Result

Extracts the names of predictor features that received non-zero
coefficients from a fitted regularized regression model produced by
[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
By default the best-performing model (as determined internally by
`run_regreg`) is used. A specific alpha value may be requested instead.

Intercept terms are excluded from the returned vector unless
`include_intercept = TRUE`. For multinomial models, where the same
feature may appear in multiple group comparisons, the union of all
non-zero features across comparisons is returned by default; set
`by_comparison = TRUE` to obtain a named list broken down by comparison.

## Usage

``` r
get_sigfeatures(
  res,
  alpha = NULL,
  include_intercept = FALSE,
  by_comparison = FALSE,
  sort_by = "name",
  x_names = NULL
)
```

## Arguments

- res:

  A list of class `"run_regreg"` produced by
  [`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).

- alpha:

  A single numeric value in `[0, 1]` specifying which fitted model to
  query. Must match one of the alpha values that was passed to
  `run_regreg`. When `NULL` (default), the best-performing model
  selected by `run_regreg` is used.

- include_intercept:

  Logical. If `TRUE`, intercept term(s) (`"(Intercept)"`) are retained
  in the output. Default: `FALSE`.

- by_comparison:

  Logical. Only relevant for multinomial classification results, which
  contain a `Comparison` column in their coefficient table. When `TRUE`,
  a named list is returned where each element contains the non-zero
  feature names for one group comparison (e.g., `"B vs A"`). When
  `FALSE` (default), the union of non-zero features across all
  comparisons is returned as a single character vector. Ignored for
  binary classification and regression results.

- sort_by:

  A single string controlling how the returned feature names are
  ordered. One of:

  `"name"`

  :   Alphabetical order (default).

  `"coef"`

  :   Descending absolute coefficient magnitude. Not available when
      `by_comparison = TRUE` and the result is multinomial (the union
      step discards magnitudes); a warning is issued and the function
      falls back to `"name"` in that case.

  `"none"`

  :   Preserve the original row order from the coefficient table.

- x_names:

  An optional character vector of column names from the original `x`
  that was passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  When supplied, only features whose names appear in `x_names` are
  returned, silently dropping any unpenalized metadata columns that were
  appended via `not_penalized`. If `NULL` (default), the function first
  checks whether `res` stores the original `x` column names internally
  (via `res$x_names`); if found, those are used automatically. If
  neither is available, no name-based filtering is applied and a note is
  issued. The simplest usage is `x_names = colnames(x)`.

## Value

- Default (`by_comparison = FALSE`):

  A character vector of feature names with non-zero coefficients, in the
  order specified by `sort_by`. Returns `character(0)` with a message
  when no non-zero features are found (e.g., the null model).

- `by_comparison = TRUE` (multinomial only):

  A named list of character vectors, one per group comparison. Each
  element is named by its comparison label (e.g., `"B vs A"`).

## Details

Extract Significant Features from a pondeR Model Result

## See also

[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## ---- Setup ---------------------------------------------------------------
set.seed(42)
n  <- 80
p  <- 30

feat <- as.data.frame(
  matrix(rnorm(n * p), nrow = n,
         dimnames = list(NULL, paste0("feat_", seq_len(p))))
)
meta <- data.frame(
  Group = sample(c("Control", "Case"), n, replace = TRUE)
)

## ---- Binary classification -----------------------------------------------
res_bin <- run_regreg(
  x    = feat, metadata = meta,
  pred = "Group", ref = "Control",
  alpha = c(0, 0.5, 1), seed = 123
)

# Best model features (default)
get_sigfeatures(res_bin)

# Features from LASSO specifically
get_sigfeatures(res_bin, alpha = 1)

# Sorted by absolute coefficient magnitude
get_sigfeatures(res_bin, sort_by = "coef")

# Include the intercept term
get_sigfeatures(res_bin, include_intercept = TRUE)

## ---- Multinomial classification ------------------------------------------
meta3 <- data.frame(
  Group = sample(c("A", "B", "C"), n, replace = TRUE)
)
res_multi <- run_regreg(
  x    = feat, metadata = meta3,
  pred = "Group", ref = "A",
  alpha = 0.5, seed = 1
)

# Union of non-zero features across all comparisons
get_sigfeatures(res_multi)

# Per-comparison breakdown
get_sigfeatures(res_multi, by_comparison = TRUE)

## ---- Gaussian regression -------------------------------------------------
meta_num <- data.frame(Score = rnorm(n))
res_num  <- run_regreg(
  x    = feat, metadata = meta_num,
  pred = "Score", alpha = 0.5, seed = 7
)
get_sigfeatures(res_num)
} # }
```
