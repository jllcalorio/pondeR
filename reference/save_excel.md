# Save Data Frames, Tibbles, or Matrices to an Excel File

Exports one or more data frames, tibbles, or matrices to a single
`.xlsx` file (or `.csv` for single-object exports). When `x` is a named
list, each tabular element is written to its own worksheet. pondeR
result objects (classes beginning with `run_`) are automatically
unpacked to extract their tabular components. If `x` contains no direct
tables but holds a list of tables one level deep, that nested list is
transparently unwrapped.

## Usage

``` r
save_excel(
  x,
  filename = "output",
  folder = NULL,
  overwrite = FALSE,
  freeze_first_row = TRUE,
  freeze_first_col = FALSE,
  auto_col_width = TRUE,
  fd_date_stamp = FALSE,
  fd_time_stamp = FALSE,
  fn_date_stamp = FALSE,
  fn_time_stamp = FALSE
)
```

## Arguments

- x:

  A data frame, tibble, matrix, a pondeR result object (e.g. output of
  [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md),
  [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)),
  or a named list whose elements are any combination of those types. One
  level of nesting is automatically unpacked when no direct tabular
  elements are found at the top level.

- filename:

  A single character string for the output file name. Used as-is when
  `x` is a single object. Ignored when `x` is a list with multiple
  elements — each element's name is used as its sheet name and the list
  object's name becomes the file name. Illegal characters are removed
  silently. Defaults to `"output"` if not supplied.

- folder:

  A single character string giving the target directory. Supports nested
  folders (e.g., `"Main folder/sub folder"`). Created recursively if
  absent. Defaults to `NULL` (current working directory).

- overwrite:

  Logical. When `TRUE` an existing file is silently overwritten. When
  `FALSE` (default) a Windows-style counter suffix is appended, e.g.
  `"results (1).xlsx"`.

- freeze_first_row:

  Logical. When `TRUE` (default) the first row of every worksheet is
  frozen.

- freeze_first_col:

  Logical. When `TRUE` the first column of every worksheet is frozen.
  Defaults to `FALSE`.

- auto_col_width:

  Logical. When `TRUE` (default) column widths are auto-fitted.

- fd_date_stamp:

  Logical. Appends `"MmmDDYYYY"` date to `folder`. Default `FALSE`.

- fd_time_stamp:

  Logical. Appends `"HHMM"` time to `folder`. Default `FALSE`.

- fn_date_stamp:

  Logical. Appends date to the file name. Default `FALSE`.

- fn_time_stamp:

  Logical. Appends time to the file name. Default `FALSE`.

## Value

Called for its side effect. Returns the full normalised file path
invisibly as a character string.

## Details

**Sheet naming.** Sheet names are taken from `names(x)`. Illegal Excel
characters are stripped, names are truncated to 25 characters (vowels
removed from the end first, then consonants if needed), the reserved
name `"History"` becomes `"History_"`, and duplicates get a `" (N)"`
suffix (up to `" (999)"`, hence the 25-character base limit).

**Single vs. multi-object behaviour.** When `x` is a single tabular
object, the supplied `filename` is used exactly (after sanitisation).
When `x` is a list with multiple elements, the `filename` argument is
ignored and each element's name is used as its sheet name; the workbook
file is named after the list variable.

**pondeR result objects.** Objects whose class starts with `run_` (e.g.
`run_auc`, `run_diff`) are unpacked automatically and their tabular
slots are written as separate sheets.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Single data frame — filename is respected
save_excel(mtcars, filename = "mtcars_export")

# Named list — filename is ignored; sheet names come from list names
sheets <- list(Cars = mtcars, Iris = iris)
save_excel(sheets, filename = "ignored")

# pondeR result object
res <- run_auc(data = mydata, response = "group", predictors = c("A","B"))
save_excel(res, filename = "auc_results")

# With date + time stamp on the file name
save_excel(iris, filename = "iris", fn_date_stamp = TRUE, fn_time_stamp = TRUE)
} # }
```
