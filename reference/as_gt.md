# Convert a run_correl to a gt table

Passes the data frame to
[`gt::gt()`](https://gt.rstudio.com/reference/gt.html) and appends the
method footnote (if present) as a source note via
[`gt::tab_source_note()`](https://gt.rstudio.com/reference/tab_source_note.html).

## Usage

``` r
as_gt(x, ...)

# S3 method for class 'run_correl'
as_gt(x, ...)
```

## Arguments

- x:

  An object to convert.

- ...:

  Additional arguments.

## Value

A `gt_tbl` object. Convert an object to a gt table

A gt table object.
