# Check if an Object Contains Infinite Values

Tests whether an object contains at least one `Inf` or `-Inf` value.
Accepts tabular objects, plain lists, or atomic vectors. Returns a
single logical value (`TRUE` or `FALSE`).

## Usage

``` r
has_inf(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), a warning is emitted when `x` is
  `NULL` or of an unsupported type.

## Value

`TRUE` if at least one `Inf` or `-Inf` is found, `FALSE` otherwise.
Never `NA`.

## Details

Both positive (`Inf`) and negative (`-Inf`) infinite values trigger a
`TRUE` result. `NaN` and `NA` are not infinite (consistent with
[`is.infinite()`](https://rdrr.io/r/base/is.finite.html) in base R). For
tabular inputs the entire table is scanned; for vectors or list inputs
all elements are scanned after unlisting.

## Author

John Lennon L. Calorio

## Examples

``` r
has_inf(c(1, Inf, 3))           # TRUE
#> [1] TRUE
has_inf(c(1, -Inf, 3))          # TRUE
#> [1] TRUE
has_inf(c(1, NaN, NA))          # FALSE - NaN/NA are not Inf
#> [1] FALSE
has_inf(c(1, 2, 3))             # FALSE
#> [1] FALSE

df <- data.frame(x = c(1, Inf), y = c(3, 4))
has_inf(df)                      # TRUE
#> [1] TRUE

mat <- matrix(c(1, 2, -Inf, 4), nrow = 2)
has_inf(mat)                     # TRUE
#> [1] TRUE
```
