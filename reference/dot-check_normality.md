# Check normality of a numeric vector

Check normality of a numeric vector

## Usage

``` r
.check_normality(
  data_vector,
  group_label = NULL,
  method = "auto",
  alpha_normality = 0.05,
  warnings_list = character(0)
)
```

## Arguments

- data_vector:

  Numeric vector to test.

- group_label:

  Optional label for the group (used in the returned label).

- method:

  One of `"auto"`, `"shapiro"`, `"lilliefors"`.

- alpha_normality:

  Significance level.

- warnings_list:

  Character vector of accumulated warnings; updated in-place via
  explicit return (no superassignment).

## Value

A named list with elements `violated`, `p_value`, `method`, `label`, and
`warnings_list`.
