# Check if an Object has Names

Tests whether an object has names. For tabular objects (`data.frame`,
tibble, `data.table`, `matrix`), column names are checked. For plain
lists, element names are checked. For atomic vectors, element names are
checked. Returns a single logical value (`TRUE` or `FALSE`).

## Usage

``` r
has_names(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test. Supported types: `data.frame`, tibble,
  `data.table`, `matrix`, plain `list`, or atomic vector.

- warn:

  Logical scalar. If `TRUE` (default), informative warnings are emitted
  describing why `x` fails the check.

## Value

A single logical value: `TRUE` if all names are present and non-empty,
`FALSE` otherwise. Never `NA`.

## Details

An object is considered named if
[`names()`](https://rdrr.io/r/base/names.html) (or
[`colnames()`](https://rdrr.io/r/base/colnames.html) for tabular
objects) returns a non-`NULL` character vector where every element is a
non-empty, non-`NA` string. A partially named object (some names are
`""` or `NA`) returns `FALSE`, and a warning is emitted identifying the
problematic positions.

`NULL` always returns `FALSE`. Unsupported object types (e.g.,
environments, S4 objects) also return `FALSE` with a warning.

## Author

John Lennon L. Calorio

## Examples

``` r
# Tabular - TRUE
has_names(mtcars)
#> [1] TRUE

# Matrix - TRUE
has_names(matrix(1:4, 2, dimnames = list(NULL, c("A", "B"))))
#> [1] TRUE

# Matrix without column names - FALSE
has_names(matrix(1:4, 2))
#> Warning: 'x' has no column names attribute.
#> [1] FALSE

# Named list - TRUE
has_names(list(a = 1, b = 2))
#> [1] TRUE

# Unnamed list - FALSE
has_names(list(1, 2, 3))
#> Warning: 'x' has no names attribute.
#> [1] FALSE

# Named vector - TRUE
has_names(c(x = 1, y = 2))
#> [1] TRUE

# Partially named list - FALSE
has_names(list(a = 1, 2, b = 3))
#> Warning: 'x' has missing or empty names at position(s): 2.
#> [1] FALSE

# NULL - FALSE
has_names(NULL)
#> Warning: 'x' is NULL.
#> [1] FALSE
```
