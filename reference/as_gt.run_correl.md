# Convert a run_correl to a gt table

Passes the data frame to
[`gt::gt()`](https://gt.rstudio.com/reference/gt.html) and appends the
method footnote (if present) as a source note via
[`gt::tab_source_note()`](https://gt.rstudio.com/reference/tab_source_note.html).

## Usage

``` r
# S3 method for class 'run_correl'
as_gt(x, ...)
```

## Arguments

- x:

  A `run_correl` object.

- ...:

  Additional arguments forwarded to
  [`gt::gt()`](https://gt.rstudio.com/reference/gt.html).

## Value

A `gt_tbl` object.
