# Plot Relative Abundance as Stacked Bar Chart

Generates a stacked bar chart of relative abundances from a feature
matrix, optionally restricting the display to the top or bottom `limit`
features by *total* abundance and collapsing the remainder into an
`"Others"` category. Samples may be shown individually or aggregated by
a grouping variable using the mean or median. Both horizontal and
vertical orientations are supported, along with full control over
colours, fonts, legend layout, and axis labels. Column names containing
spaces or special characters are fully supported.

## Usage

``` r
plot_relabund(
  x,
  metadata,
  limit = 10L,
  sort = NULL,
  others_position = NULL,
  show_x_as = NULL,
  group_aggregate = NULL,
  normalize_to_relabund = TRUE,
  ignore_normalization = FALSE,
  order = NULL,
  angle = NULL,
  abundance_in_y = TRUE,
  legend_position = NULL,
  legend_nrow = NULL,
  legend_ncol = NULL,
  italicize_legend_items = FALSE,
  exclude_others = FALSE,
  unknown_before_others = TRUE,
  theme = "nature",
  colors = NULL,
  plot_title = "",
  plot_subtitle = "",
  legend_title = "",
  xlab = NULL,
  ylab = "Relative abundance",
  global_font_size = 20,
  font_family = "sans",
  title_size = NULL,
  subtitle_size = NULL,
  legend_size = NULL,
  xlab_size = NULL,
  ylab_size = NULL,
  ...
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or `matrix` with named columns representing
  features (e.g., taxa, metabolites). Rows are samples; columns are
  features. Column names may contain spaces and special characters.
  **When `group_aggregate = NULL` and `ignore_normalization = FALSE`**,
  `x` is expected to already contain row-wise relative abundances (each
  row summing to 1). A warning is issued for any column whose sum is not
  0 or 1, indicating that column-wise relative abundance normalisation
  may not have been applied.

- metadata:

  A `data.frame` containing sample-level metadata. Must have the same
  number of rows as `x`; row order is assumed to correspond to `x`.

- limit:

  Integer or `NULL`. Number of features to display, selected by *total*
  abundance across all samples. When `sort = "high"` (or
  `NULL`/`"alphanumeric"` with a limit), the top `limit` highest-total
  features are shown and the rest collapsed to `"Others"`. When
  `sort = "low"`, the bottom `limit` lowest-total features are shown
  instead. If `NULL`, all features are plotted. Default is `10`.

- sort:

  Character or `NULL`. Controls both which features are selected (when
  `limit` is set) and the order features appear in the legend/stack.
  `NULL` (default) preserves the original column order of `x`. Options:
  `"alphanumeric"` (lexicographic), `"high"` (highest total abundance
  first — selects top `limit`), `"low"` (lowest total abundance first —
  selects bottom `limit`). **Note:** `"high"` and `"low"` rank by column
  totals and are meaningful only when `x` contains raw counts or
  intensities (row sums not already 0 or 1). If normalised data is
  detected, a warning is issued and `sort` is ignored. When using ranked
  sorting, supply raw count data and set `normalize_to_relabund = TRUE`.

- others_position:

  Character or `NULL`. Where to place the `"Others"` category in the
  legend/stack order. `NULL` (default) treats `"Others"` according to
  `sort`. `"first"` forces it to the top; `"last"` forces it to the
  bottom.

- show_x_as:

  Character or `NULL`. Column name in `metadata` to use as the x-axis
  variable (sample IDs or group labels). If `NULL`, the function
  searches for `"UID"`, `"ID"`, or `"SampleID"` (case-insensitive) in
  `metadata`. An error is raised if none are found.

- group_aggregate:

  Character or `NULL`. If `"mean"` or `"median"`, `show_x_as` is treated
  as a grouping variable and raw feature values in `x` are aggregated
  per group before normalising to relative abundances. **Assumes `x`
  contains raw (unnormalised) data.** See `normalize_to_relabund`. If
  `NULL` (default), each row in `x` is an individual sample and `x` must
  already contain row-wise relative abundances (unless
  `ignore_normalization = TRUE`).

- normalize_to_relabund:

  Logical. Relevant when `group_aggregate` is non-`NULL` or `sort` is
  `"high"`/`"low"`. If `TRUE` (default),
  `run_normalize(method = "col_rel_abundance")` is applied to `x` before
  aggregation or ranking. Set to `FALSE` if `x` is already appropriately
  normalised. Ignored when `ignore_normalization = TRUE`.

- ignore_normalization:

  Logical. If `TRUE`, all normalisation checks and pre-normalisation
  steps are bypassed and `x` is plotted as supplied. Default is `FALSE`.
  **Note:** Setting this to `TRUE` means the plotted values are not
  relative abundances; use this when plotting already-transformed data
  such as CLR-transformed values where the data originated from
  normalised data and should be visualised as-is.

- order:

  Character vector, `"alphanumeric"`, or `NULL`. Left-to-right order of
  samples/groups on the x-axis. A character vector gives an explicit
  sequence. `"alphanumeric"` sorts lexicographically. `NULL` (default)
  preserves natural order of appearance in `metadata`.

- angle:

  Numeric or `NULL`. Rotation angle (0–360 degrees) for sample/group
  tick labels. `NULL` (default): `0` when `group_aggregate` is
  non-`NULL`, `45` otherwise.

- abundance_in_y:

  Logical. If `TRUE` (default), relative abundance is on the y-axis and
  samples/groups are on the x-axis. If `FALSE`, axes are flipped:
  relative abundance on the x-axis, samples/groups on the y-axis. The
  legend order is also reversed when `FALSE` so the top-of-stack feature
  appears first (top) in the legend.

- legend_position:

  Character or `NULL`. One of `"top"`, `"bottom"`, `"left"`, `"right"`.
  Auto-detected when `NULL`: `"top"` for `abundance_in_y = TRUE`;
  `"right"` for `abundance_in_y = FALSE`.

- legend_nrow:

  Integer or `NULL`. Rows in legend key grid.

- legend_ncol:

  Integer or `NULL`. Columns in legend key grid.

- italicize_legend_items:

  Logical. Render legend labels in italic. Requires `ggtext`; falls back
  to `element_text(face = "italic")` with a warning. Default `FALSE`.

- exclude_others:

  Logical. If `TRUE`, the `"Others"` category is omitted from the plot
  entirely. Default `FALSE`.

- unknown_before_others:

  Logical. If `TRUE` (default), any feature whose name is exactly
  `"unknown"` (case-insensitive, e.g., `"Unknown"`, `"UNKNOWN"`,
  `"unknown"`) is placed immediately before `"Others"` in the
  legend/stack, or last when `exclude_others = TRUE`. Partial matches
  are not considered; only column names that are *entirely* the word
  "unknown" are affected.

- theme:

  Character. ggplot2 theme. One of `"nature"`, `"minimal"`, `"classic"`,
  `"bw"`, `"light"`, `"dark"`. Default `"nature"`.

- colors:

  Character vector or `NULL`. Custom fill colours. Length must equal the
  number of plotted fill levels. Named vectors matched by name. `NULL`
  uses the Okabe-Ito palette (extended with HCL colours).

- plot_title:

  Character. Main title. Supports `ggtext` markdown: `**bold**`,
  `*italic*`. Default `""`.

- plot_subtitle:

  Character. Subtitle. Supports same markdown. Default `""`.

- legend_title:

  Character. Legend title. Default `""`.

- xlab:

  Character or `NULL`. X-axis label. Defaults to `show_x_as`.

- ylab:

  Character. Y-axis label. Default `"Relative abundance"`.

- global_font_size:

  Numeric or `NULL`. Scales all text elements proportionally relative to
  an internal 11 pt reference. Default `20`.

- font_family:

  Character. Font family for all text. Default `"sans"`.

- title_size:

  Numeric or `NULL`. Title font size (pt).

- subtitle_size:

  Numeric or `NULL`. Subtitle font size (pt).

- legend_size:

  Numeric or `NULL`. Legend title font size (pt).

- xlab_size:

  Numeric or `NULL`. X-axis label font size (pt).

- ylab_size:

  Numeric or `NULL`. Y-axis label font size (pt).

- ...:

  Additional arguments passed to
  [`ggplot2::geom_bar()`](https://ggplot2.tidyverse.org/reference/geom_bar.html).

## Value

A `ggplot` object, also drawn to the active graphics device.

## Details

Plot Relative Abundance as Stacked Bar Chart

**Special characters in column names.** Column names containing spaces,
slashes, parentheses, or other special characters are fully supported.
The function avoids
[`stats::reshape()`](https://rdrr.io/r/stats/reshape.html) (which can
mangle such names) and instead builds the long-format data frame
manually.

**Feature selection with `limit`.** Features are ranked by their *column
total* across all samples (or aggregated groups). When `sort = "high"`
or `NULL` or `"alphanumeric"`, the top `limit` features by total are
retained; all others are collapsed into `"Others"`. When `sort = "low"`,
the bottom `limit` features by total are retained instead.

**Normalisation flow.** When `ignore_normalization = FALSE` and
`normalize_to_relabund = TRUE`, row-relative abundance is applied via
`run_normalize(method = "row_rel_abundance")` before any aggregation or
ranking. Row-wise proportions are then computed so each bar sums to
100\\ `ignore_normalization = TRUE`, the data are plotted as supplied
and the y-axis label should be interpreted accordingly.

**Legend reversal when flipped.** When `abundance_in_y = FALSE`,
`coord_flip()` causes the stack to render left-to-right. The legend fill
levels are reversed so the first (top-of-stack) feature also appears
first (top) in the legend, maintaining visual correspondence.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
library(pondeR)
set.seed(42)
n    <- 20
taxa <- c("Firmicutes", "Bacteroidetes", "Proteobacteria (class A)",
          "Actino/bacteroides", "Verrucomicrobia", "Fusobacteria sp.",
          "Tenericutes", "Spirochaetes", "Unknown", "Cyanobacteria",
          "Chloroflexi", "Acidobacteria")

counts <- as.data.frame(
  matrix(sample(0L:5000L, n * length(taxa), replace = TRUE),
         nrow = n, dimnames = list(NULL, taxa))
)
meta <- data.frame(
  SampleID = paste0("S", seq_len(n)),
  Group    = rep(c("Control", "Treatment"), each = n / 2L)
)

# --- Pipeline: col_rel_abundance → CLR → plot as-is ---
norm_out <- run_normalize(counts, meta, method = "col_rel_abundance",
                          verbose = FALSE)
clr_out  <- run_transform(norm_out$data, method = "clr")

plot_relabund(
  x                    = clr_out$data,
  metadata             = meta,
  limit                = 10L,
  ignore_normalization = TRUE,
  plot_title           = "Gut Microbiome (CLR-transformed)"
)

# --- Group-aggregated from raw counts, bold title ---
plot_relabund(
  x               = counts,
  metadata        = meta,
  show_x_as       = "Group",
  group_aggregate = "mean",
  order           = c("Treatment", "Control"),
  limit           = 8L,
  sort            = "high",
  plot_title      = "**Relative Abundance** by Group",
  legend_title    = "Phylum"
)

# --- Horizontal bars ---
plot_relabund(
  x                     = counts,
  metadata              = meta,
  show_x_as             = "Group",
  group_aggregate       = "mean",
  limit                 = 6L,
  exclude_others        = FALSE,
  unknown_before_others = TRUE,
  abundance_in_y        = FALSE
)
} # }
```
