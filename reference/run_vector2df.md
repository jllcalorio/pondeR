# Convert a Named Vector to a Two-Column Data Frame

Transforms a named vector into a tidy two-column data frame, where one
column holds the names and the other holds the corresponding values.
Optionally filters rows by a value threshold or range. Useful for
quickly tabulating named outputs such as VIP scores, loadings, or any
named numeric/character vector.

## Usage

``` r
run_vector2df(
  x,
  name = "Variable",
  value = "Value",
  sort_by = c("none", "name", "value_asc", "value_desc"),
  min = NULL,
  max = NULL,
  range_type = c("inner", "outer"),
  stringsAsFactors = FALSE
)
```

## Arguments

- x:

  A named vector (numeric, integer, character, or logical). Must have
  names; an error is raised if `names(x)` is `NULL`.

- name:

  A single character string giving the name of the column that will
  contain the vector names. Defaults to `"Variable"`.

- value:

  A single character string giving the name of the column that will
  contain the vector values. Defaults to `"Value"`.

- sort_by:

  A single character string controlling row ordering. One of `"none"`
  (original order), `"name"` (alphabetical by name column),
  `"value_asc"` (ascending by value), or `"value_desc"` (descending by
  value). Defaults to `"none"`.

- min:

  A single numeric value used as a lower threshold for filtering the
  value column. Behaviour depends on `range_type`:

  - `"inner"`: rows with `value >= min` are kept (when `max` is `NULL`),
    or rows satisfying `value >= min & value <= max` (when both are
    supplied).

  - `"outer"`: rows with `value <= min` are kept (when `max` is `NULL`),
    or rows satisfying `value <= min | value >= max` (when both are
    supplied).

  Defaults to `NULL` (no lower filtering).

- max:

  A single numeric value used as an upper threshold for filtering the
  value column. Behaviour depends on `range_type`:

  - `"inner"`: rows with `value <= max` are kept (when `min` is `NULL`),
    or rows satisfying `value >= min & value <= max` (when both are
    supplied).

  - `"outer"`: rows with `value >= max` are kept (when `min` is `NULL`),
    or rows satisfying `value <= min | value >= max` (when both are
    supplied).

  Defaults to `NULL` (no upper filtering).

- range_type:

  A single character string, either `"inner"` or `"outer"`, that
  controls how `min` and `max` are applied when both are supplied.
  `"inner"` keeps values *between* the thresholds (inclusive); `"outer"`
  keeps values *outside* the thresholds (inclusive). Ignored when only
  one of `min`/`max` is supplied. Defaults to `"inner"`.

- stringsAsFactors:

  Logical. If `TRUE`, the name column is returned as a `factor`.
  Defaults to `FALSE`.

## Value

A `data.frame` with two columns named according to `name` and `value`,
and one row per (retained) element of `x`. Row names are reset to
sequential integers.

## Author

John Lennon L. Calorio

## Examples

``` r
vip <- c(`623.0149@9.699` = 0.3616770, `769.4391@9.699` = 1.0151602,
         `1147.7305@9.699` = 1.2068827, `397.3033@9.7028` = 0.5754963)
run_vector2df(vip, name = "Feature", value = "VIP",
              min = 1.0, sort_by = "value_desc")
#>           Feature      VIP
#> 1 1147.7305@9.699 1.206883
#> 2  769.4391@9.699 1.015160
```
