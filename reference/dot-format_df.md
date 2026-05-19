# Format a degrees-of-freedom string consistently across test types

For t-tests / Wilcoxon the `parameter` slot is a single named number;
for ANOVA-like tests it is a named vector of length 2 or more.

## Usage

``` r
.format_df(parameter)
```

## Arguments

- parameter:

  Named numeric vector from an `htest` or equivalent.

## Value

A single character string, e.g. `"df = 18"` or `"DFn = 2, DFd = 27"`.
