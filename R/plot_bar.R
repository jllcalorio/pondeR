#' @title Create Publication-Ready Bar Plots
#'
#' @description
#' Generate highly customizable bar plots using categorical variables from a
#' data frame, tibble, or matrix. The function supports stacked, dodged,
#' and filled bar plots, optional grouping, automatic color-blind friendly
#' palettes, percentage or count displays, optional coordinate flipping,
#' and automatic handling of missing numeric axis levels.
#'
#' @param x A data frame, tibble, or matrix containing the variables to plot.
#'
#' @param by A single character string specifying a categorical grouping
#' variable used to split the bars. Default is `NULL`.
#'
#' @param plot_type Character string specifying the bar plot type.
#' Supported options are:
#' \itemize{
#'   \item `"stacked"`: stacked bar plot
#'   \item `"dodged"`: side-by-side grouped bars
#'   \item `"filled"`: proportional stacked bar plot
#' }
#' Default is `"stacked"`.
#'
#' @param colors Either `"auto"` for automatic Okabe-Ito color palette generation
#' or a character vector of valid color names or HEX color codes.
#' If `by` is not `NULL`, the length of `colors` must match the number
#' of unique categories in `by`.
#'
#' @param x_axis A single character string specifying the column name
#' to use on the x-axis.
#'
#' @param show_y_as Character string specifying y-axis display:
#' \itemize{
#'   \item `"n"`: counts
#'   \item `"p"`: percentages
#' }
#' Default is `"n"`.
#'
#' @param show_values Logical. If `TRUE`, count or percentage labels are
#' displayed on the bars. Default is `FALSE`.
#'
#' @param fill_missing_levels Logical. Default is `TRUE`.
#' If `x_axis` is numeric, all integer values between the minimum and maximum
#' observed values will be displayed on the axis, including values with zero
#' counts. This is not a missing-value imputation procedure.
#'
#' @param arrange_x Optional vector specifying the desired ordering of x-axis
#' categories. Categories not listed will generate an error.
#' Missing values can be included using `NA`.
#'
#' @param x_axis_label Optional character string for x-axis label.
#' Defaults to `x_axis`.
#'
#' @param y_axis_label Optional character string for y-axis label.
#' Defaults to `"Count"` or `"Percent"` depending on `show_y_as`.
#'
#' @param title Optional character string specifying plot title.
#'
#' @param subtitle Optional character string specifying plot subtitle.
#'
#' @param caption Optional character string specifying plot caption.
#'
#' @param legend Optional character string specifying legend title.
#' Defaults to `by`.
#'
#' @param legend_position Character string specifying legend position.
#' Supported options are `"bottom"`, `"top"`, `"left"`, and `"right"`.
#' Default is `"bottom"`.
#'
#' @param flip Logical. If `TRUE`, flips the plot coordinates using
#' `ggplot2::coord_flip()`. Default is `FALSE`.
#'
#' @param base_size Base font size in points; every text element scales from
#' it. Default is `11`, which suits a default 7 x 5 inch Quarto figure.
#'
#' @param title_size,subtitle_size,caption_size,legend_size,xlab_size,ylab_size,axis_text_size,label_size
#' Font sizes in points for the title, subtitle, caption, legend, x-axis label,
#' y-axis label, tick labels and value labels. Default `NULL` derives each from
#' `base_size`.
#'
#' @param theme Character string specifying ggplot2 theme style.
#' Supported options are:
#' \itemize{
#'   \item `"nature"`
#'   \item `"minimal"`
#'   \item `"classic"`
#'   \item `"bw"`
#'   \item `"light"`
#'   \item `"dark"`
#' }
#' Default is `"nature"`.
#'
#' @details
#' This function is designed for publication-ready categorical bar plots with
#' minimal user setup. Numeric x-axis variables are treated as categorical
#' while preserving numeric ordering.
#'
#' The default `"auto"` color mode uses the Okabe-Ito color-blind friendly
#' palette:
#'
#' \preformatted{
#' "#0072B2", "#D55E00", "#009E73", "#56B4E9",
#' "#E69F00", "#CC79A7", "#F0E442", "#000000"
#' }
#'
#' Missing values are explicitly displayed as `"NA"`.
#'
#' Gridlines are removed by default for a cleaner publication-ready appearance.
#'
#' @return
#' A `ggplot2` object representing the generated bar plot.
#'
#' @references
#' Wickham H (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}.
#'   Springer-Verlag New York. \doi{10.1007/978-3-319-24277-4}
#'
#' Okabe M, Ito K (2008). \emph{Color Universal Design (CUD): How to make
#'   figures and presentations that are friendly to colorblind people}.
#'   \url{https://jfly.uni-koeln.de/color/}
#'
#' @examples
#' # Basic count plot
#' plot_bar(
#'   x = mtcars,
#'   x_axis = "cyl"
#' )
#'
#' # Dodged grouped plot
#' plot_bar(
#'   x = as.data.frame(Titanic),
#'   x_axis = "Class",
#'   by = "Sex",
#'   plot_type = "dodged"
#' )
#'
#' # Numeric x-axis with missing levels shown
#' set.seed(1)
#' df <- data.frame(
#'   Age = c(18:38, 40),
#'   Group = sample(c("A", "B"), 22, TRUE)
#' )
#'
#' plot_bar(
#'   x = df,
#'   x_axis = "Age",
#'   by = "Group"
#' )
#'
#' @author
#' John Lennon L. Calorio
#'
#' @export
plot_bar <- function(
    x,
    by = NULL,
    plot_type = "stacked",
    colors = "auto",
    x_axis,
    show_y_as = "n",
    show_values = FALSE,
    fill_missing_levels = TRUE,
    arrange_x = NULL,
    x_axis_label = NULL,
    y_axis_label = NULL,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    legend = NULL,
    legend_position = "bottom",
    flip = FALSE,
    base_size = 11,
    title_size = NULL,
    subtitle_size = NULL,
    caption_size = NULL,
    legend_size = NULL,
    xlab_size = NULL,
    ylab_size = NULL,
    axis_text_size = NULL,
    label_size = NULL,
    theme = "nature"
) {
  sizes <- .text_sizes(base_size, title_size, subtitle_size, xlab_size, ylab_size,
                       axis_text_size, legend_size, caption_size, label = label_size)

  # ===========================================================================
  # Input validation
  # ===========================================================================

  if (!inherits(x, c("data.frame", "matrix"))) {
    stop(
      "`x` must be a data.frame, tibble, or matrix.",
      call. = FALSE
    )
  }

  x <- as.data.frame(x, check.names = FALSE)

  if (nrow(x) == 0L) {
    stop(
      "`x` contains zero rows.",
      call. = FALSE
    )
  }

  if (!is.character(x_axis) ||
      length(x_axis) != 1L ||
      is.na(x_axis)) {
    stop(
      "`x_axis` must be a single non-missing character string.",
      call. = FALSE
    )
  }

  if (!x_axis %in% colnames(x)) {
    stop(
      sprintf("Column '%s' was not found in `x`.", x_axis),
      call. = FALSE
    )
  }

  if (!is.null(by)) {

    if (!is.character(by) ||
        length(by) != 1L ||
        is.na(by)) {
      stop(
        "`by` must be a single non-missing character string.",
        call. = FALSE
      )
    }

    if (!by %in% colnames(x)) {
      stop(
        sprintf("Column '%s' was not found in `x`.", by),
        call. = FALSE
      )
    }
  }

  valid_plot_types <- c(
    "stacked",
    "dodged",
    "filled"
  )

  if (!plot_type %in% valid_plot_types) {
    stop(
      sprintf(
        "`plot_type` must be one of: %s",
        paste(valid_plot_types, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  valid_show_y <- c("n", "p")

  if (!show_y_as %in% valid_show_y) {
    stop(
      "`show_y_as` must be either 'n' or 'p'.",
      call. = FALSE
    )
  }

  valid_legend_positions <- c(
    "bottom",
    "top",
    "left",
    "right"
  )

  if (!legend_position %in% valid_legend_positions) {
    stop(
      sprintf(
        "`legend_position` must be one of: %s",
        paste(valid_legend_positions, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.logical(show_values) ||
      length(show_values) != 1L ||
      is.na(show_values)) {
    stop(
      "`show_values` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(fill_missing_levels) ||
      length(fill_missing_levels) != 1L ||
      is.na(fill_missing_levels)) {
    stop(
      "`fill_missing_levels` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(flip) ||
      length(flip) != 1L ||
      is.na(flip)) {
    stop(
      "`flip` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  # ===========================================================================
  # Extract variables
  # ===========================================================================

  x_var <- x[[x_axis]]

  if (all(is.na(x_var))) {
    stop(
      sprintf(
        "Column '%s' contains only missing values.",
        x_axis
      ),
      call. = FALSE
    )
  }

  # ===========================================================================
  # Determine ordering
  # ===========================================================================

  if (is.numeric(x_var)) {

    finite_values <- x_var[is.finite(x_var)]

    if (length(finite_values) == 0L) {
      stop(
        sprintf(
          "Column '%s' contains no finite numeric values.",
          x_axis
        ),
        call. = FALSE
      )
    }

    if (isTRUE(fill_missing_levels)) {

      ordered_levels <- seq(
        from = floor(min(finite_values)),
        to = ceiling(max(finite_values)),
        by = 1
      )

    } else {

      ordered_levels <- sort(
        unique(finite_values)
      )
    }

  } else {

    ordered_levels <- sort(
      unique(as.character(x_var)),
      na.last = TRUE
    )
  }

  ordered_levels <- as.character(ordered_levels)

  # ===========================================================================
  # Explicit NA handling
  # ===========================================================================

  x_var_chr <- as.character(x_var)
  x_var_chr[is.na(x_var_chr)] <- "NA"

  if (!is.null(arrange_x)) {

    arrange_x_chr <- as.character(arrange_x)
    arrange_x_chr[is.na(arrange_x_chr)] <- "NA"

    missing_levels <- setdiff(
      unique(x_var_chr),
      arrange_x_chr
    )

    if (length(missing_levels) > 0L) {
      stop(
        sprintf(
          "The following categories are missing from `arrange_x`: %s",
          paste(missing_levels, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    ordered_levels <- arrange_x_chr

  } else {

    if ("NA" %in% x_var_chr) {
      ordered_levels <- c(
        ordered_levels,
        "NA"
      )
    }
  }

  x_factor <- factor(
    x_var_chr,
    levels = ordered_levels
  )

  # ===========================================================================
  # Grouping variable
  # ===========================================================================

  if (!is.null(by)) {

    by_var <- as.character(x[[by]])
    by_var[is.na(by_var)] <- "NA"

    group_levels <- unique(by_var)

  } else {

    by_var <- NULL
    group_levels <- NULL
  }

  # ===========================================================================
  # Color handling
  # ===========================================================================

  okabe_ito <- c(
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#56B4E9",
    "#E69F00",
    "#CC79A7",
    "#F0E442",
    "#000000"
  )

  if (!is.null(by)) {

    n_groups <- length(unique(group_levels))

    if (identical(colors, "auto")) {

      colors <- rep(
        okabe_ito,
        length.out = n_groups
      )

    } else {

      if (!is.character(colors)) {
        stop(
          "`colors` must be a character vector.",
          call. = FALSE
        )
      }

      if (length(colors) != n_groups) {
        stop(
          sprintf(
            "Length of `colors` (%d) must equal the number of `%s` groups (%d).",
            length(colors),
            by,
            n_groups
          ),
          call. = FALSE
        )
      }
    }
  }

  # ===========================================================================
  # Build summary table
  # ===========================================================================

  if (is.null(by)) {

    summary_df <- as.data.frame(
      table(
        x = x_factor,
        useNA = "no"
      )
    )

    colnames(summary_df) <- c(
      "x",
      "count"
    )

    summary_df$prop <- summary_df$count /
      sum(summary_df$count)

  } else {

    summary_df <- as.data.frame(
      table(
        x = x_factor,
        group = by_var,
        useNA = "no"
      )
    )

    colnames(summary_df) <- c(
      "x",
      "group",
      "count"
    )

    total_by_x <- stats::ave(
      summary_df$count,
      summary_df$x,
      FUN = sum
    )

    summary_df$prop <- summary_df$count /
      total_by_x
  }

  # ===========================================================================
  # Labels
  # ===========================================================================

  summary_df$label <- if (show_y_as == "n") {

    summary_df$count

  } else {

    paste0(
      round(summary_df$prop * 100, 1),
      "%"
    )
  }

  # ===========================================================================
  # Axis labels
  # ===========================================================================

  if (is.null(x_axis_label)) {
    x_axis_label <- x_axis
  }

  # Filled bars always show proportions, whatever `show_y_as` is
  y_is_pct <- show_y_as == "p" || plot_type == "filled"

  if (is.null(y_axis_label)) {
    y_axis_label <- if (y_is_pct) "Percent" else "Count"
  }

  if (is.null(legend)) {
    legend <- by
  }

  # ===========================================================================
  # Position setting
  # ===========================================================================

  position_setting <- switch(
    plot_type,
    stacked = "stack",
    dodged = ggplot2::position_dodge(width = 0.9),
    filled = "fill"
  )

  # ---------------------------------------------------------------------------
  # Plot construction
  # ---------------------------------------------------------------------------

  p <- if (is.null(by)) {
    ggplot2::ggplot(summary_df,
                    ggplot2::aes(x = x, y = if (show_y_as == "n") count else prop)) +
      ggplot2::geom_col(fill = okabe_ito[1])
  } else {
    ggplot2::ggplot(summary_df,
                    ggplot2::aes(x = x, y = if (show_y_as == "n") count else prop,
                                 fill = group)) +
      ggplot2::geom_col(position = position_setting) +
      ggplot2::scale_fill_manual(values = colors, name = legend)
  }

  # ---------------------------------------------------------------------------
  # Value labels: centred in stacked segments, above dodged or single bars
  # ---------------------------------------------------------------------------

  if (isTRUE(show_values)) {
    lab_aes <- if (is.null(by)) ggplot2::aes(label = label) else
      ggplot2::aes(label = label, group = group)
    centred <- !is.null(by) && plot_type != "dodged"
    lab_pos <- if (is.null(by)) "identity" else switch(
      plot_type,
      stacked = ggplot2::position_stack(vjust = 0.5),
      dodged  = ggplot2::position_dodge(width = 0.9),
      filled  = ggplot2::position_fill(vjust = 0.5)
    )
    p <- p + ggplot2::geom_text(
      data = summary_df[summary_df$count > 0, , drop = FALSE],
      mapping = lab_aes,
      position = lab_pos,
      vjust = if (centred || flip) 0.5 else -0.3,
      hjust = if (!centred && flip) -0.1 else 0.5,
      size = .pt_to_mm(sizes$label),
      show.legend = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Axes, labels and theme
  # ---------------------------------------------------------------------------

  p <- p +
    ggplot2::scale_y_continuous(
      labels = if (y_is_pct) function(v) paste0(round(v * 100), "%") else ggplot2::waiver(),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    )

  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }

  p +
    ggplot2::labs(
      x = x_axis_label,
      y = y_axis_label,
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    .ponder_theme(theme, sizes) +
    ggplot2::theme(
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}
