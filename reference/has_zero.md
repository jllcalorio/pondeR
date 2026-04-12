# Check if an Object Contains Zero Values

Tests whether an object contains at least one element equal to `0`.
Accepts tabular objects, plain lists, or atomic vectors. Returns a
single logical value (`TRUE` or `FALSE`).

## Usage

``` r
has_zero(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), a warning is emitted when `x` is
  `NULL` or of an unsupported type.

## Value

`TRUE` if at least one `0` is found among numeric elements, `FALSE`
otherwise. Never `NA`.

## Details

Only numeric (or coercible-to-numeric) elements are tested. Non-numeric
elements (e.g., character, logical) are silently skipped. `NA` values
are not counted as zero.

## Author

John Lennon L. Calorio

## Examples

``` r
has_zero(c(1, 0, 3))            # TRUE
#> [1] TRUE
has_zero(c(1, 2, 3))            # FALSE
#> [1] FALSE
has_zero(c(NA, 0, 1))           # TRUE - NA is not zero, but 0 is present
#> [1] TRUE
has_zero(c(NA, NA))             # FALSE
#> [1] FALSE

df <- data.frame(x = c(1, 0), y = c(3, 4))
has_zero(df)                     # TRUE
#> [1] TRUE
```
