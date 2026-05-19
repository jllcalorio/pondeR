# Run normality checks in parallel (or sequentially)

Encapsulates the Windows / POSIX branching so it is not duplicated
between `run_diff` and `.run_diff_single`.

## Usage

``` r
.parallel_normality(
  groups,
  y,
  grp,
  normality_method,
  alpha_normality,
  num_cores
)
```

## Arguments

- groups:

  Character vector of group labels.

- y:

  Numeric outcome vector.

- grp:

  Factor grouping vector (same length as `y`).

- normality_method:

  Passed to `.check_normality`.

- alpha_normality:

  Passed to `.check_normality`.

- num_cores:

  Number of cores to use.

## Value

A list of normality-check result lists, one per group. Each element
carries its own `warnings_list` slot; callers must harvest these.
