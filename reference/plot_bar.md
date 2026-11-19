# Create Publication-Ready Bar Plots

Generate highly customizable bar plots using categorical variables from
a data frame, tibble, or matrix. The function supports stacked, dodged,
and filled bar plots, optional grouping, automatic color-blind friendly
palettes, percentage or count displays, optional coordinate flipping,
and automatic handling of missing numeric axis levels.

## Usage

``` r
plot_bar(
  x,
  by = NULL,
  plot_type = "stacked",
  colors = "auto",
  x_axis,
  show_y_as = "n",
  show_values = FALSE,
  values_position = "auto",
  fill_missing_levels = TRUE,
  arrange_x = NULL,
  x_axis_label = NULL,
  y_axis_label = NULL,
  title = "",
  subtitle = "",
  caption = "",
  legend = NULL,
  legend_position = "bottom",
  flip = FALSE,
  global_font_size = 25,
  title_size = global_font_size,
  subtitle_size = global_font_size,
  caption_size = global_font_size,
  legend_size = global_font_size,
  axis_title_size = global_font_size,
  axis_text_size = global_font_size,
  theme = "nature",
  ...
)
```

## Arguments

- x:

  A data frame, tibble, or matrix containing the variables to plot.

- by:

  A single character string specifying a categorical grouping variable
  used to split the bars. Default is `NULL`.

- plot_type:

  Character string specifying the bar plot type. Supported options are:

  - `"stacked"`: stacked bar plot

  - `"dodged"`: side-by-side grouped bars

  - `"filled"`: proportional stacked bar plot

  Default is `"stacked"`.

- colors:

  Either `"auto"` for automatic Okabe-Ito color palette generation or a
  character vector of valid color names or HEX color codes. If `by` is
  not `NULL`, the length of `colors` must match the number of unique
  categories in `by`.

- x_axis:

  A single character string specifying the column name to use on the
  x-axis.

- show_y_as:

  Character string specifying y-axis display:

  - `"n"`: counts

  - `"p"`: percentages

  Default is `"n"`.

- show_values:

  Logical. If `TRUE`, count or percentage labels are displayed on the
  bars. Default is `FALSE`.

- values_position:

  Character string specifying label placement. Supported options are:

  - `"auto"`: automatic placement

  - `"inside"`: inside bars

  - `"outside"`: outside bars

  Default is `"auto"`.

- fill_missing_levels:

  Logical. Default is `TRUE`. If `x_axis` is numeric, all integer values
  between the minimum and maximum observed values will be displayed on
  the axis, including values with zero counts. This is not a
  missing-value imputation procedure.

- arrange_x:

  Optional vector specifying the desired ordering of x-axis categories.
  Categories not listed will generate an error. Missing values can be
  included using `NA`.

- x_axis_label:

  Optional character string for x-axis label. Defaults to `x_axis`.

- y_axis_label:

  Optional character string for y-axis label. Defaults to `"Count"` or
  `"Percent"` depending on `show_y_as`.

- title:

  Optional character string specifying plot title.

- subtitle:

  Optional character string specifying plot subtitle.

- caption:

  Optional character string specifying plot caption.

- legend:

  Optional character string specifying legend title. Defaults to `by`.

- legend_position:

  Character string specifying legend position. Supported options are
  `"bottom"`, `"top"`, `"left"`, and `"right"`. Default is `"bottom"`.

- flip:

  Logical. If `TRUE`, flips the plot coordinates using
  [`ggplot2::coord_flip()`](https://ggplot2.tidyverse.org/reference/coord_flip.html).
  Default is `FALSE`.

- global_font_size:

  A single positive numeric value specifying the base font size for the
  plot. Default is `25`.

- title_size:

  A single positive numeric value specifying title font size. Defaults
  to `global_font_size`.

- subtitle_size:

  A single positive numeric value specifying subtitle font size.
  Defaults to `global_font_size`.

- caption_size:

  A single positive numeric value specifying caption font size. Defaults
  to `global_font_size`.

- legend_size:

  A single positive numeric value specifying legend font size. Defaults
  to `global_font_size`.

- axis_title_size:

  A single positive numeric value specifying axis title font size.
  Defaults to `global_font_size`.

- axis_text_size:

  A single positive numeric value specifying axis text font size.
  Defaults to `global_font_size`.

- theme:

  Character string specifying ggplot2 theme style. Supported options
  are:

  - `"nature"`

  - `"minimal"`

  - `"classic"`

  - `"bw"`

  - `"light"`

  - `"dark"`

  Default is `"nature"`.

- ...:

  Additional arguments passed to
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `ggplot2` object representing the generated bar plot.

## Details

This function is designed for publication-ready categorical bar plots
with minimal user setup. Numeric x-axis variables are treated as
categorical while preserving numeric ordering.

The default `"auto"` color mode uses the Okabe-Ito color-blind friendly
palette:

    "#0072B2", "#D55E00", "#009E73", "#56B4E9",
    "#E69F00", "#CC79A7", "#F0E442", "#000000"

Missing values are explicitly displayed as `"NA"`.

Gridlines are removed by default for a cleaner publication-ready
appearance.

Requires the packages:

- ggplot2

- scales

- ggrepel

## References

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

Okabe M, Ito K (2008). *Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people*.
<https://jfly.uni-koeln.de/color/>

## Author

John Lennon L. Calorio

## Examples

``` r
# Basic count plot
plot_bar(
  x = mtcars,
  x_axis = "cyl"
)


# Dodged grouped plot
plot_bar(
  x = as.data.frame(Titanic),
  x_axis = "Class",
  by = "Sex",
  plot_type = "dodged"
)


# Numeric x-axis with missing levels shown
set.seed(1)
df <- data.frame(
  Age = c(18:38, 40),
  Group = sample(c("A", "B"), 22, TRUE)
)

plot_bar(
  x = df,
  x_axis = "Age",
  by = "Group"
)

```
