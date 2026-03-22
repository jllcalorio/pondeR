# Reduce Multiple Data Frames by Common or Unique Column/Row Names

Accepts two or more data frames (or tibbles or matrices with dimnames)
and returns a single data frame whose columns (or rows) are either the
*intersection* (names present in **all** supplied objects) or the
*union* (names present in **any** supplied object) of the column or row
names across all inputs.

## Usage

``` r
run_reduce(..., check_cols = TRUE, return_common = TRUE)
```

## Arguments

- ...:

  Two or more objects coercible to a data frame: `data.frame`, `tbl_df`
  (tibble), or named `matrix`. If exactly one object is supplied a
  warning is issued and that object is returned as-is (after coercion to
  `data.frame`).

- check_cols:

  Logical scalar (default `TRUE`). When `TRUE` the reduction is
  performed over **column** names; when `FALSE` it is performed over
  **row** names. All supplied objects must carry the relevant dimnames
  when `check_cols = FALSE`.

- return_common:

  Logical scalar (default `TRUE`). When `TRUE` only names present in
  **all** data frames are retained (intersection). When `FALSE` names
  present in **any** data frame are retained (union); missing entries
  are filled with `NA`.

## Value

A named list of class `"run_reduce"` with the following elements:

- `data`:

  A `data.frame` containing the reduced result.

- `summary`:

  A named list with one element per input data frame (named `df1`,
  `df2`, ...). Each element is a character vector of the column (or row)
  names that were **removed** from that input during the reduction, or
  `character(0)` if nothing was removed.

- `call`:

  The matched call, for reproducibility.

- `check_cols`:

  Logical; mirrors the `check_cols` argument.

- `return_common`:

  Logical; mirrors the `return_common` argument.

- `n_inputs`:

  Integer; the number of data frames supplied.

- `target_names`:

  Character vector of the names retained in the final result.

## Details

### Single input

Passing a single data frame raises a warning and returns that frame
unchanged (coerced to `data.frame`). No reduction is possible.

### Column-wise reduction (`check_cols = TRUE`)

The set of column names to keep is computed as the intersection or union
of [`colnames()`](https://rdrr.io/r/base/colnames.html) across all
inputs. Each frame is then subset (or expanded with `NA` columns for the
union case) to the target set, and the results are row-bound with
`do.call(rbind, ...)`.

### Row-wise reduction (`check_cols = FALSE`)

The set of row names to keep is computed as the intersection or union of
[`rownames()`](https://rdrr.io/r/base/colnames.html) across all inputs.
Each frame is then subset (or expanded with `NA` rows for the union
case) to the target set, and the results are column-bound with
`do.call(cbind, ...)`. All inputs must have non-`NULL`, non-empty, and
non-duplicated row names.

### Duplicate names

Duplicated column names within a single data frame, or duplicated row
names when `check_cols = FALSE`, will trigger an error because the
subsetting behaviour would be ambiguous.

### Performance note

Name intersection/union is computed once via `Reduce(intersect, ...)` or
`Reduce(union, ...)` — an idiomatic base-R approach that avoids repeated
pairwise comparisons and scales well with the number of inputs.

## Author

John Lennon L. Calorio

## Examples

``` r
## --- Example data ---------------------------------------------------
df1 <- data.frame(a = 1:3, b = 4:6,  c = 7:9)
df2 <- data.frame(a = 10:12, b = 13:15, d = 16:18)
df3 <- data.frame(a = 19:21, b = 22:24, e = 25:27)

## 1. Return columns common to ALL data frames (default)
result_common <- run_reduce(df1, df2, df3)
result_common$data          # columns: a, b
#>    a  b
#> 1  1  4
#> 2  2  5
#> 3  3  6
#> 4 10 13
#> 5 11 14
#> 6 12 15
#> 7 19 22
#> 8 20 23
#> 9 21 24
result_common$summary       # columns dropped per frame
#> $df1
#> [1] "c"
#> 
#> $df2
#> [1] "d"
#> 
#> $df3
#> [1] "e"
#> 

## 2. Return columns from ANY data frame (union, fill with NA)
result_union <- run_reduce(df1, df2, df3, return_common = FALSE)
result_union$data           # columns: a, b, c, d, e  (NAs where absent)
#>    a  b  c  d  e
#> 1  1  4  7 NA NA
#> 2  2  5  8 NA NA
#> 3  3  6  9 NA NA
#> 4 10 13 NA 16 NA
#> 5 11 14 NA 17 NA
#> 6 12 15 NA 18 NA
#> 7 19 22 NA NA 25
#> 8 20 23 NA NA 26
#> 9 21 24 NA NA 27

## 3. Row-wise reduction -------------------------------------------
m1 <- data.frame(x = c(1, 2, 3), row.names = c("r1", "r2", "r3"))
m2 <- data.frame(y = c(4, 5, 6), row.names = c("r1", "r3", "r4"))

result_rows <- run_reduce(m1, m2, check_cols = FALSE)
result_rows$data            # rows: r1, r3 (common to both)
#>    x y
#> r1 1 4
#> r3 3 5

result_rows_union <- run_reduce(m1, m2, check_cols = FALSE,
                                return_common = FALSE)
result_rows_union$data      # rows: r1, r2, r3, r4 (NAs where absent)
#>     x  y
#> r1  1  4
#> r2  2 NA
#> r3  3  5
#> r4 NA  6

## 4. Single data frame (warning, returned as-is) ------------------
result_single <- run_reduce(df1)
#> Warning: Only one data frame was supplied. Returning it unchanged (no reduction performed).
```
