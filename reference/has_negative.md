# Check if an Object Contains Negative Values

Tests whether an object contains at least one element with a value
strictly less than `0`. Accepts tabular objects, plain lists, or atomic
vectors. Returns a single logical value (`TRUE` or `FALSE`).

## Usage

``` r
has_negative(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), a warning is emitted when `x` is
  `NULL` or of an unsupported type.

## Value

`TRUE` if at least one strictly negative value is found, `FALSE`
otherwise. Never `NA`.

## Details

Only numeric (or coercible-to-numeric) elements are tested. `NA` values
are excluded from the comparison. `-Inf` is treated as negative.

## Author

John Lennon L. Calorio

## Examples

``` r
has_negative(c(1, -2, 3))       # TRUE
#> [1] TRUE
has_negative(c(1, 2, 3))        # FALSE
#> [1] FALSE
has_negative(c(0, -0.001, NA))  # TRUE
#> [1] TRUE
has_negative(c(-Inf, 1, 2))     # TRUE
#> [1] TRUE

df <- data.frame(x = c(1, -1), y = c(3, 4))
has_negative(df)                 # TRUE
#> [1] TRUE
```
