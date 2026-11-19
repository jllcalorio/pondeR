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

  A string; `"within_match"` (default), `"within_group"`, or `"all"`.
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

  - **`"all"`**: Every `match` x `group` combination has exactly the
    same number of observations (the size of the smallest cell). Use
    with `independent = FALSE` so multiple `match` columns are crossed.
    If a combination is absent from the data, nothing can be selected.

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

## Author

John Lennon L. Calorio

## Examples

``` r
set.seed(1)
d <- data.frame(
  id = 1:300,
  `Age Group` = sample(c("2-5", "6-10"), 300, TRUE, prob = c(.6, .4)),
  Sex = sample(c("Male", "Female"), 300, TRUE),
  Status = sample(c("Stunted", "Non-stunted"), 300, TRUE, prob = c(.3, .7)),
  check.names = FALSE
)
tab <- function(r) with(r$chosen, table(`Age Group`, Sex, Status))

# 1. Simple random subset of 50 rows
nrow(run_subset(d, n = 50)$chosen)
#> [1] 50

# 2. Joint Age Group x Sex cells; Status is equal within every cell
#    (Stunted == Non-stunted in each Age Group x Sex combination), but
#    different cells may have different sizes.
tab(run_subset(d, match = c("Age Group", "Sex"), group = "Status",
               independent = FALSE))
#> , , Status = Non-stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      25   34
#>      6-10     13   24
#> 
#> , , Status = Stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      25   34
#>      6-10     13   24
#> 

# 3. Default independent = TRUE balances each match column separately and
#    unions the rows, so the joint cells are NOT guaranteed to be balanced.
tab(run_subset(d, match = c("Age Group", "Sex"), group = "Status"))
#> , , Status = Non-stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      36   58
#>      6-10     25   26
#> 
#> , , Status = Stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      25   34
#>      6-10     13   24
#> 

# 4. Flip the direction: within each Status, Age Group x Sex cells are equal
tab(run_subset(d, match = c("Age Group", "Sex"), group = "Status",
               independent = FALSE, balance = "within_group"))
#> , , Status = Non-stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      33   33
#>      6-10     33   33
#> 
#> , , Status = Stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      13   13
#>      6-10     13   13
#> 

# 4b. Fully equal: all Age Group x Sex x Status cells have the same count
tab(run_subset(d, match = c("Age Group", "Sex"), group = "Status",
               independent = FALSE, balance = "all"))
#> , , Status = Non-stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      13   13
#>      6-10     13   13
#> 
#> , , Status = Stunted
#> 
#>          Sex
#> Age Group Female Male
#>      2-5      13   13
#>      6-10     13   13
#> 

# 5. Cap the size; a warning is raised if n cannot be reached
nrow(run_subset(d, n = 40, match = "Sex", group = "Status")$chosen)
#> [1] 40

# 6. Deterministic: first rows of each stratum
run_subset(d, n = 10, match = "Sex", random = FALSE)$chosen
#>    id Age Group    Sex      Status
#> 1   1       2-5   Male     Stunted
#> 2   2       2-5   Male     Stunted
#> 3   3       2-5   Male Non-stunted
#> 4   4      6-10   Male     Stunted
#> 5   5       2-5   Male     Stunted
#> 7   7      6-10 Female Non-stunted
#> 8   8      6-10 Female Non-stunted
#> 9   9      6-10 Female     Stunted
#> 10 10       2-5 Female Non-stunted
#> 20 20      6-10 Female Non-stunted

# 7. Column subsetting
names(run_subset(d, n = 2, by = "col", random = FALSE)$chosen)
#> [1] "id"        "Age Group"
```
