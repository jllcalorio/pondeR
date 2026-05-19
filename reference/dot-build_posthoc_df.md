# Build a standardised post-hoc result data frame

Accepts the raw output from any of
[`rstatix::tukey_hsd`](https://rpkgs.datanovia.com/rstatix/reference/tukey_hsd.html),
[`rstatix::games_howell_test`](https://rpkgs.datanovia.com/rstatix/reference/games_howell_test.html),
or
[`rstatix::dunn_test`](https://rpkgs.datanovia.com/rstatix/reference/dunn_test.html)
and adds `n1`, `n2`, `mean1`, `mean2`, formatted p columns,
`significant`, `p.adj.signif`, and `interpretation`.

## Usage

``` r
.build_posthoc_df(
  ph,
  posthoc_label,
  y,
  grp,
  group_ns,
  groups,
  parametric,
  test_alpha,
  p_adjust_method
)
```

## Arguments

- ph:

  Raw data frame from rstatix.

- posthoc_label:

  Character string naming the post-hoc test used.

- y:

  Numeric outcome vector.

- grp:

  Factor grouping vector.

- group_ns:

  Named integer table of per-group n.

- groups:

  Character vector of group names (levels).

- parametric:

  Logical — use mean-based interpretation when `TRUE`.

- test_alpha:

  Significance threshold.

- p_adjust_method:

  p-adjustment method (used only for Dunn's; others already carry p.adj
  from rstatix).

## Value

A data frame with a standardised column set.
