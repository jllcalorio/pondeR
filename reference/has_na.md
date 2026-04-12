# Check if an Object Contains Missing Values

Tests whether an object contains at least one `NA` or `NaN` value.
Accepts tabular objects (`data.frame`, tibble, `data.table`, `matrix`),
plain lists, or atomic vectors (including individual table columns).
Returns a single logical value (`TRUE` or `FALSE`).

## Usage

``` r
has_na(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), a warning is emitted when `x` is
  `NULL` or of an unsupported type.

## Value

`TRUE` if at least one `NA` or `NaN` is found, `FALSE` otherwise. Never
`NA`.

## Details

`NaN` is treated as a missing value (consistent with
`is.na(NaN) == TRUE` in base R). Zero (`0`) is *not* considered missing.
For tabular inputs the entire table is scanned; for vectors or list
inputs all elements are scanned after unlisting.

## Author

John Lennon L. Calorio

## Examples

``` r
has_na(c(1, 2, NA, 4))          # TRUE
#> [1] TRUE
has_na(c(1, 2, NaN, 4))         # TRUE
#> [1] TRUE
has_na(c(1, 2, 0, 4))           # FALSE - 0 is not NA
#> [1] FALSE
has_na(c(1, 2, 3))              # FALSE
#> [1] FALSE

df <- data.frame(x = c(1, NA), y = c(3, 4))
has_na(df)                       # TRUE
#> [1] TRUE

has_na(mtcars)                   # FALSE
#> [1] FALSE
```
