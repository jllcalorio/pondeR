# Balanced Subsetting of Rows or Columns

Subsets a data frame, tibble, or matrix either row-wise or column-wise,
with optional balancing across categorical groupings. When balancing is
requested via `match`, the function ensures that every category within
each specified column is represented by the same number of observations.
An optional `n` argument acts as a row ceiling: the per-category count
is derived so that the chosen set is as close to `n` as possible while
remaining perfectly balanced within each `match` column.

When `match` contains more than one column, the `independent` argument
controls how they are combined:

- **`independent = TRUE` (default):** Each `match` column is balanced
  separately. A row is selected if it satisfies the balance requirement
  of **any** `match` column (union). This maximises the number of chosen
  rows.

- **`independent = FALSE`:** All `match` columns are crossed into a
  single Cartesian product of categories (e.g. `"1-2 years × Male"`),
  and every unique combination is represented equally. This guarantees
  perfect joint balance.

## Usage

``` r
run_subset(
  x,
  n = NULL,
  match = NULL,
  independent = TRUE,
  balance = "within_match",
  group = NULL,
  by = "row",
  random = TRUE,
  top_to_bottom = TRUE,
  left_to_right = TRUE,
  seed = 123
)
```

## Arguments

- x:

  A data frame, tibble, or matrix to be subsetted.

- n:

  An integer specifying the desired number of rows (when `by = "row"`)
  or columns (when `by = "col"`) in `$chosen`.

- match:

  `NULL` (default) or a character vector of one or more column names in
  `x` containing categorical variables to balance on.

- independent:

  Logical; `TRUE` (default). When `match` has multiple columns: If
  `TRUE`, each column is balanced independently and their results are
  unioned. If `FALSE`, columns are crossed into a single interaction and
  balanced jointly.

- balance:

  A string; either `"within_match"` (default) or `"within_group"`.
  Controls the direction of balancing when both `match` and `group` are
  provided.

  - **`"within_match"`**: For every category in `match`, the levels of
    `group` will have an equal number of observations. For example, if
    balancing by Age Group, each Age Group will have an equal number of
    Diseased and Healthy patients.

  - **`"within_group"`**: For every category in `group`, the levels of
    `match` will have an equal number of observations. For example,
    within the Diseased group, there will be an equal number of patients
    from each Age Group.

- group:

  `NULL` (default) or a character vector of column names that define
  additional grouping strata.

- by:

  A string; either `"row"` (default) or `"col"`.

- random:

  Logical; `TRUE` (default) to select rows or columns randomly.

- top_to_bottom:

  Logical; `TRUE` (default). For deterministic row pulls.

- left_to_right:

  Logical; `TRUE` (default). For deterministic col pulls.

- seed:

  An integer; `123` (default). The random seed. `NULL` disables seeding.

## Value

A named list of class `"run_subset"`.
