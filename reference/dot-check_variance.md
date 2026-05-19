# Run variance-homogeneity test (Levene or Bartlett fallback)

Run variance-homogeneity test (Levene or Bartlett fallback)

## Usage

``` r
.check_variance(y, grp, alpha_variance, verbose, warnings_list)
```

## Arguments

- y:

  Numeric outcome vector.

- grp:

  Factor grouping vector.

- alpha_variance:

  Significance level.

- verbose:

  Logical.

- warnings_list:

  Character vector of accumulated warnings.

## Value

A list with `p_value`, `test_name`, `violated`, and `warnings_list`.
