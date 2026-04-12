# Check if an Object is a Tabular Data Structure

Tests whether an object belongs to a recognized tabular data class -
`data.frame`, `tbl_df` (tibble), `data.table`, or `matrix` - and
confirms that it has exactly two dimensions. Returns a single logical
value (`TRUE` or `FALSE`).

## Usage

``` r
is_tabular(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), informative warnings are emitted
  when `x` fails the check, describing the specific reason for failure.
  Set to `FALSE` to suppress all warnings and obtain a silent `FALSE`.

## Value

A single logical value: `TRUE` if `x` is a recognized tabular structure
with two dimensions, `FALSE` otherwise. The return value is never `NA`.

## Details

An object is considered tabular if it satisfies **both** of the
following conditions:

1.  Its class is one of `"data.frame"`, `"tbl_df"`, `"tbl"` (generic
    tibble superclass), `"data.table"`, or `"matrix"`.

2.  It has exactly two dimensions (i.e., `length(dim(x)) == 2`).

The class check uses [`inherits()`](https://rdrr.io/r/base/class.html),
so subclasses of `data.frame` (e.g., tibbles, which inherit from
`"data.frame"`) are naturally detected. A `matrix` always has two
dimensions by construction, but the dimension check is applied uniformly
for consistency and to guard against unusual edge cases such as
zero-dimensional matrices produced by certain external packages.

`NULL` and objects with no class attribute are handled gracefully: they
always return `FALSE`.

Note that `is_tabular()` does **not** check whether the object is
non-empty. A zero-row or zero-column tabular object will still return
`TRUE`, as it retains the correct class and dimensions. Callers that
require non-empty data should perform that check separately.

## Author

John Lennon L. Calorio

## Examples

``` r
# --- TRUE cases ---

# Base data.frame
is_tabular(mtcars)
#> [1] TRUE

# Matrix
is_tabular(matrix(1:12, nrow = 3))
#> [1] TRUE

# tibble (if tibble is installed)
if (requireNamespace("tibble", quietly = TRUE)) {
  is_tabular(tibble::tibble(x = 1:3, y = letters[1:3]))
}
#> [1] TRUE

# data.table (if data.table is installed)
if (requireNamespace("data.table", quietly = TRUE)) {
  is_tabular(data.table::data.table(a = 1:5, b = rnorm(5)))
}
#> [1] TRUE

# Zero-row data.frame - still tabular
is_tabular(mtcars[0, ])
#> [1] TRUE

# --- FALSE cases ---

# Atomic vector
is_tabular(1:10)
#> Warning: 'x' has class <integer>, which is not a recognized tabular type. Expected one of: data.frame, tibble (tbl_df), data.table, or matrix.
#> [1] FALSE

# List
is_tabular(list(a = 1, b = 2))
#> Warning: 'x' has class <list>, which is not a recognized tabular type. Expected one of: data.frame, tibble (tbl_df), data.table, or matrix.
#> [1] FALSE

# 3-D array (more than two dimensions)
is_tabular(array(1:24, dim = c(2, 3, 4)))
#> Warning: 'x' has class <array>, which is not a recognized tabular type. Expected one of: data.frame, tibble (tbl_df), data.table, or matrix.
#> [1] FALSE

# NULL
is_tabular(NULL)
#> Warning: 'x' is NULL. A tabular object (data.frame, tibble, data.table, or matrix) was expected.
#> [1] FALSE

# Suppress warnings for programmatic use
is_tabular("not a table", warn = FALSE)
#> [1] FALSE
```
