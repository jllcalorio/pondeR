#' Plot Odds Ratios
#'
#' @description
#' Creates a publication-ready odds ratio (forest) plot from a data frame
#' containing odds ratios and optional confidence interval bounds. Points are
#' colored according to whether the odds ratio is below, equal to, or above 1
#' (the null value). A vertical reference line is drawn at OR = 1. When
#' confidence interval columns are supplied, horizontal error bars are added.
#'
#' The function accepts data frames, matrices, and tibbles, and safely handles
#' column names containing special characters (spaces, hyphens, parentheses,
#' etc.).
#'
#' @param x A data frame, matrix, or tibble containing at minimum a column of
#'   odds ratios. Column names with special characters are supported.
#' @param or A single character string giving the name of the numeric column in
#'   \code{x} that contains the odds ratios. **Required.** Values must be
#'   strictly positive (> 0).
#' @param low Optional. A single character string giving the name of the numeric
#'   column in \code{x} that contains the lower confidence interval bound (e.g.,
#'   lower 95% CI). Must be supplied together with \code{high}; supplying only
#'   one of \code{low} / \code{high} raises an error.
#' @param high Optional. A single character string giving the name of the
#'   numeric column in \code{x} that contains the upper confidence interval
#'   bound (e.g., upper 95% CI). Must be supplied together with \code{low};
#'   supplying only one of \code{low} / \code{high} raises an error.
#' @param ci Numeric scalar (default \code{95}). The confidence level expressed
#'   as a percentage (e.g., \code{95}, \code{99}). Used only in axis and legend
#'   labels — it does not compute any CI values.
#' @param feature_col Optional. A single character string giving the name of the
#'   column in \code{x} to use as feature/row labels on the y-axis. If
#'   \code{NULL} (default), row names of \code{x} are used; if \code{x} has no
#'   row names, sequential integers are used instead.
#' @param colors A character vector of length 3 specifying the colors for odds
#'   ratios that are (1) less than 1, (2) equal to 1, and (3) greater than 1,
#'   respectively. Defaults to \code{c("blue", "grey", "red")}.
#' @param dot_size Numeric scalar. Size of the point (dot) for each odds ratio.
#'   Defaults to \code{3}.
#' @param line_size Numeric scalar. Line width of the confidence interval error
#'   bars. Only used when both \code{low} and \code{high} are supplied. Defaults
#'   to \code{0.8}.
#' @param theme Character string specifying the ggplot2 theme. One of
#'   \code{"nature"} (default), \code{"bw"}, \code{"classic"}, \code{"minimal"},
#'   \code{"gray"}, or \code{"dark"}.
#' @param global_font_size Numeric scalar. Base font size (in pt) applied to all
#'   text elements in the plot (axis text, axis titles, plot title, subtitle,
#'   caption, legend). Defaults to \code{15}.
#' @param plot_title Character string. Main title of the plot. Defaults to
#'   \code{NULL} (no title).
#' @param plot_subtitle Character string. Subtitle of the plot. Defaults to
#'   \code{NULL} (no subtitle).
#' @param x_label Character string. Label for the x-axis. When \code{NULL}
#'   (default), the label is auto-generated as \code{"Odds Ratio"} (when no CI
#'   columns are given) or \code{"Odds Ratio (95\% CI)"} (or the corresponding
#'   percentage from \code{ci}).
#' @param y_label Character string. Label for the y-axis. Defaults to
#'   \code{NULL} (no label).
#' @param font_family Character string. Font family for all text elements.
#'   Defaults to \code{"Helvetica"}.
#'
#' @return A \code{ggplot2} object that can be further customized, saved with
#'   \code{\link[ggplot2]{ggsave}}, or passed to \code{save_plots()}.
#'
#' @details
#' **Color logic**
#'
#' Each odds ratio point is colored by comparing its value to 1:
#' \itemize{
#'   \item OR < 1  → first color in \code{colors} (default \code{"blue"})
#'   \item OR = 1  → second color in \code{colors} (default \code{"grey"})
#'   \item OR > 1  → third color in \code{colors} (default \code{"red"})
#' }
#' A small numerical tolerance (\code{.Machine$double.eps^0.5}) is applied when
#' testing for exact equality to 1.
#'
#' **Axis centering**
#'
#' The x-axis is expanded symmetrically around 1 on the log scale to ensure
#' that 1 always sits at the visual center of the plot. This makes it
#' immediately clear whether each OR falls in the "decreased odds" or
#' "increased odds" half of the panel.
#'
#' **Missing CI values**
#'
#' Individual rows where \code{low} and/or \code{high} are \code{NA} will have
#' their error bars silently omitted; the OR dot is still plotted.
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{run_logreg}}, \code{\link{plot_volcano}},
#'   \code{\link{plot_diff}}
#'
#' @examples
#' # ── Basic usage: OR only ─────────────────────────────────────────────────────
#' df_basic <- data.frame(
#'   feature  = c("Age", "BMI", "Smoking", "Hypertension", "Diabetes"),
#'   OR       = c(1.45, 0.82, 2.31, 0.55, 1.02),
#'   stringsAsFactors = FALSE
#' )
#' rownames(df_basic) <- df_basic$feature
#'
#' plot_oddsratio(x = df_basic, or = "OR")
#'
#' # ── With confidence intervals ────────────────────────────────────────────────
#' df_ci <- data.frame(
#'   feature  = c("Age", "BMI", "Smoking", "Hypertension", "Diabetes"),
#'   OR       = c(1.45, 0.82, 2.31, 0.55, 1.02),
#'   Lower95  = c(1.10, 0.61, 1.70, 0.31, 0.78),
#'   Upper95  = c(1.91, 1.11, 3.14, 0.97, 1.33),
#'   stringsAsFactors = FALSE
#' )
#'
#' plot_oddsratio(
#'   x    = df_ci,
#'   or   = "OR",
#'   low  = "Lower95",
#'   high = "Upper95",
#'   ci   = 95,
#'   feature_col = "feature",
#'   plot_title  = "Logistic Regression: Odds Ratios"
#' )
#'
#' # ── Special-character column names ───────────────────────────────────────────
#' df_special <- data.frame(
#'   check.names = FALSE,
#'   "OR (adjusted)" = c(0.73, 1.88, 1.05),
#'   "Lower 95% CI"  = c(0.51, 1.20, 0.79),
#'   "Upper 95% CI"  = c(1.05, 2.94, 1.39)
#' )
#' rownames(df_special) <- c("Metabolite A", "Metabolite B", "Metabolite C")
#'
#' plot_oddsratio(
#'   x    = df_special,
#'   or   = "OR (adjusted)",
#'   low  = "Lower 95% CI",
#'   high = "Upper 95% CI"
#' )
#'
#' # ── Custom colors and theme ──────────────────────────────────────────────────
#' plot_oddsratio(
#'   x      = df_ci,
#'   or     = "OR",
#'   low    = "Lower95",
#'   high   = "Upper95",
#'   feature_col = "feature",
#'   colors = c("#2166AC", "#878787", "#D6604D"),
#'   theme  = "classic",
#'   global_font_size = 13,
#'   plot_title = "Forest Plot",
#'   x_label    = "Odds Ratio (95% CI)",
#'   y_label    = "Predictor"
#' )
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline
#'   scale_color_manual scale_x_log10 coord_cartesian labs theme_bw
#'   theme_classic theme_minimal theme_gray theme_dark theme element_text
#'   element_line element_blank margin expansion
#' @importFrom rlang sym
#'
#' @export
plot_oddsratio <- function(
    x,
    or,
    low             = NULL,
    high            = NULL,
    ci              = 95,
    feature_col     = NULL,
    colors          = c("blue", "grey", "red"),
    dot_size        = 3,
    line_size       = 0.8,
    theme           = "nature",
    global_font_size = 15,
    plot_title      = NULL,
    plot_subtitle   = NULL,
    x_label         = NULL,
    y_label         = NULL,
    font_family     = "Helvetica"
) {

  # ── 0. Internal helpers ──────────────────────────────────────────────────────

  .resolve_font <- function(family) {
    if (is.null(family) || !nzchar(family)) "sans" else family
  }

  .apply_nature_theme <- function(p, fs, ff) {
    p + ggplot2::theme_bw(base_size = fs, base_family = ff) +
      ggplot2::theme(
        panel.grid.major  = ggplot2::element_blank(),
        panel.grid.minor  = ggplot2::element_blank(),
        panel.border      = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
        axis.ticks        = ggplot2::element_line(color = "black", linewidth = 0.4),
        axis.text         = ggplot2::element_text(size = fs, family = ff, color = "black"),
        axis.title        = ggplot2::element_text(size = fs, family = ff, color = "black"),
        plot.title        = ggplot2::element_text(size = fs + 2, family = ff, face = "bold",
                                                  hjust = 0, margin = ggplot2::margin(b = 4)),
        plot.subtitle     = ggplot2::element_text(size = fs, family = ff, color = "grey30",
                                                  hjust = 0, margin = ggplot2::margin(b = 6)),
        legend.text       = ggplot2::element_text(size = fs, family = ff),
        legend.title      = ggplot2::element_text(size = fs, family = ff, face = "bold"),
        legend.position   = "bottom",
        legend.key        = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        strip.text        = ggplot2::element_text(size = fs, family = ff),
        plot.margin       = ggplot2::margin(10, 15, 10, 10)
      )
  }

  .apply_theme <- function(p, theme_name, fs, ff) {
    switch(
      theme_name,
      "nature"  = .apply_nature_theme(p, fs, ff),
      "bw"      = p + ggplot2::theme_bw(base_size = fs, base_family = ff) +
        ggplot2::theme(
          axis.text     = ggplot2::element_text(size = fs, family = ff),
          axis.title    = ggplot2::element_text(size = fs, family = ff),
          plot.title    = ggplot2::element_text(size = fs + 2, family = ff, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = fs, family = ff),
          legend.text   = ggplot2::element_text(size = fs, family = ff),
          legend.title  = ggplot2::element_text(size = fs, family = ff, face = "bold"),
          legend.position = "bottom"
        ),
      "classic" = p + ggplot2::theme_classic(base_size = fs, base_family = ff) +
        ggplot2::theme(
          axis.text     = ggplot2::element_text(size = fs, family = ff),
          axis.title    = ggplot2::element_text(size = fs, family = ff),
          plot.title    = ggplot2::element_text(size = fs + 2, family = ff, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = fs, family = ff),
          legend.text   = ggplot2::element_text(size = fs, family = ff),
          legend.title  = ggplot2::element_text(size = fs, family = ff, face = "bold"),
          legend.position = "bottom"
        ),
      "minimal" = p + ggplot2::theme_minimal(base_size = fs, base_family = ff) +
        ggplot2::theme(
          axis.text     = ggplot2::element_text(size = fs, family = ff),
          axis.title    = ggplot2::element_text(size = fs, family = ff),
          plot.title    = ggplot2::element_text(size = fs + 2, family = ff, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = fs, family = ff),
          legend.text   = ggplot2::element_text(size = fs, family = ff),
          legend.title  = ggplot2::element_text(size = fs, family = ff, face = "bold"),
          legend.position = "bottom"
        ),
      "gray"    = p + ggplot2::theme_gray(base_size = fs, base_family = ff) +
        ggplot2::theme(
          axis.text     = ggplot2::element_text(size = fs, family = ff),
          axis.title    = ggplot2::element_text(size = fs, family = ff),
          plot.title    = ggplot2::element_text(size = fs + 2, family = ff, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = fs, family = ff),
          legend.text   = ggplot2::element_text(size = fs, family = ff),
          legend.title  = ggplot2::element_text(size = fs, family = ff, face = "bold"),
          legend.position = "bottom"
        ),
      "dark"    = p + ggplot2::theme_dark(base_size = fs, base_family = ff) +
        ggplot2::theme(
          axis.text     = ggplot2::element_text(size = fs, family = ff),
          axis.title    = ggplot2::element_text(size = fs, family = ff),
          plot.title    = ggplot2::element_text(size = fs + 2, family = ff, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = fs, family = ff),
          legend.text   = ggplot2::element_text(size = fs, family = ff),
          legend.title  = ggplot2::element_text(size = fs, family = ff, face = "bold"),
          legend.position = "bottom"
        ),
      stop(
        sprintf(
          paste0(
            "plot_oddsratio: 'theme' must be one of \"nature\", \"bw\", \"classic\", ",
            "\"minimal\", \"gray\", or \"dark\". Got: \"%s\"."
          ),
          theme_name
        ),
        call. = FALSE
      )
    )
  }

  # ── 1. Validate 'x' ─────────────────────────────────────────────────────────

  if (missing(x) || is.null(x)) {
    stop("plot_oddsratio: 'x' is required but was not provided.", call. = FALSE)
  }

  if (!(is.data.frame(x) || is.matrix(x))) {
    stop(
      "plot_oddsratio: 'x' must be a data frame, matrix, or tibble. ",
      sprintf("Got an object of class: \"%s\".", paste(class(x), collapse = ", ")),
      call. = FALSE
    )
  }

  # Coerce to data frame to allow uniform downstream access
  x <- as.data.frame(x, check.names = FALSE)

  if (nrow(x) == 0L) {
    stop("plot_oddsratio: 'x' has zero rows. Nothing to plot.", call. = FALSE)
  }

  # ── 2. Validate 'or' ────────────────────────────────────────────────────────

  if (missing(or) || is.null(or)) {
    stop(
      "plot_oddsratio: 'or' is required. Provide the name of the column in 'x' ",
      "that contains the odds ratios.",
      call. = FALSE
    )
  }

  if (!is.character(or) || length(or) != 1L || !nzchar(or)) {
    stop(
      "plot_oddsratio: 'or' must be a single non-empty character string naming ",
      "a column in 'x'.",
      call. = FALSE
    )
  }

  if (!or %in% names(x)) {
    stop(
      sprintf(
        "plot_oddsratio: Column \"%s\" specified in 'or' was not found in 'x'.\n  Available columns: %s",
        or,
        paste(sprintf('"%s"', names(x)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  or_vals <- x[[or]]

  if (!is.numeric(or_vals)) {
    stop(
      sprintf(
        "plot_oddsratio: The 'or' column \"%s\" must be numeric. Got class: \"%s\".",
        or,
        class(or_vals)
      ),
      call. = FALSE
    )
  }

  n_or_na <- sum(is.na(or_vals))
  if (n_or_na == nrow(x)) {
    stop(
      "plot_oddsratio: All values in the 'or' column are NA. Nothing to plot.",
      call. = FALSE
    )
  }
  if (n_or_na > 0L) {
    warning(
      sprintf(
        "plot_oddsratio: %d row(s) in the 'or' column contain NA and will be excluded from the plot.",
        n_or_na
      ),
      call. = FALSE
    )
  }

  non_na_or <- or_vals[!is.na(or_vals)]
  if (any(non_na_or <= 0, na.rm = TRUE)) {
    n_neg <- sum(non_na_or <= 0, na.rm = TRUE)
    stop(
      sprintf(
        paste0(
          "plot_oddsratio: %d value(s) in the 'or' column are <= 0. ",
          "Odds ratios must be strictly positive (> 0). ",
          "Check whether your column contains raw ORs rather than log-ORs."
        ),
        n_neg
      ),
      call. = FALSE
    )
  }

  # ── 3. Validate 'low' and 'high' ────────────────────────────────────────────

  has_ci <- FALSE

  if (!is.null(low) || !is.null(high)) {

    # Must supply both or neither
    if (is.null(low) || is.null(high)) {
      stop(
        "plot_oddsratio: 'low' and 'high' must both be provided or both be NULL. ",
        if (is.null(low)) "'low' is missing." else "'high' is missing.",
        call. = FALSE
      )
    }

    # Type checks
    for (arg_name in c("low", "high")) {
      arg_val <- get(arg_name)
      if (!is.character(arg_val) || length(arg_val) != 1L || !nzchar(arg_val)) {
        stop(
          sprintf(
            "plot_oddsratio: '%s' must be a single non-empty character string naming a column in 'x'.",
            arg_name
          ),
          call. = FALSE
        )
      }
    }

    # Column existence checks
    for (arg_name in c("low", "high")) {
      col_nm <- get(arg_name)
      if (!col_nm %in% names(x)) {
        stop(
          sprintf(
            "plot_oddsratio: Column \"%s\" specified in '%s' was not found in 'x'.\n  Available columns: %s",
            col_nm, arg_name,
            paste(sprintf('"%s"', names(x)), collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }

    low_vals  <- x[[low]]
    high_vals <- x[[high]]

    # Numeric checks
    if (!is.numeric(low_vals)) {
      stop(
        sprintf(
          "plot_oddsratio: The 'low' column \"%s\" must be numeric. Got class: \"%s\".",
          low, class(low_vals)
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(high_vals)) {
      stop(
        sprintf(
          "plot_oddsratio: The 'high' column \"%s\" must be numeric. Got class: \"%s\".",
          high, class(high_vals)
        ),
        call. = FALSE
      )
    }

    # Logical checks: low should not exceed high
    valid_rows <- !is.na(low_vals) & !is.na(high_vals)
    if (any(valid_rows)) {
      inverted <- low_vals[valid_rows] > high_vals[valid_rows]
      if (any(inverted)) {
        warning(
          sprintf(
            paste0(
              "plot_oddsratio: %d row(s) have 'low' > 'high', which suggests the CI bounds ",
              "may be swapped. Those error bars will still be drawn."
            ),
            sum(inverted)
          ),
          call. = FALSE
        )
      }
    }

    has_ci <- TRUE
  }

  # ── 4. Validate 'ci' ────────────────────────────────────────────────────────

  if (!is.numeric(ci) || length(ci) != 1L || ci <= 0 || ci >= 100) {
    stop(
      "plot_oddsratio: 'ci' must be a single numeric value between 0 and 100 ",
      "(exclusive), e.g., 95 or 99.",
      call. = FALSE
    )
  }

  # ── 5. Validate 'colors' ────────────────────────────────────────────────────

  if (!is.character(colors) || length(colors) != 3L) {
    stop(
      "plot_oddsratio: 'colors' must be a character vector of length 3. ",
      "Element 1 = OR < 1, element 2 = OR = 1, element 3 = OR > 1.",
      call. = FALSE
    )
  }

  # Test colors are valid R colors
  invalid_colors <- character(0)
  for (col in colors) {
    tryCatch(
      grDevices::col2rgb(col),
      error = function(e) {
        invalid_colors <<- c(invalid_colors, col)
      }
    )
  }
  if (length(invalid_colors) > 0L) {
    stop(
      sprintf(
        "plot_oddsratio: The following value(s) in 'colors' are not recognized R colors: %s.\n  Use color names (e.g., \"red\") or hex codes (e.g., \"#D6604D\").",
        paste(sprintf('"%s"', invalid_colors), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ── 6. Validate scalar numerics ─────────────────────────────────────────────

  if (!is.numeric(dot_size)  || length(dot_size)  != 1L || dot_size  <= 0) {
    stop("plot_oddsratio: 'dot_size' must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(line_size) || length(line_size) != 1L || line_size <= 0) {
    stop("plot_oddsratio: 'line_size' must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(global_font_size) || length(global_font_size) != 1L || global_font_size <= 0) {
    stop("plot_oddsratio: 'global_font_size' must be a single positive numeric value.", call. = FALSE)
  }

  # ── 7. Validate 'feature_col' ────────────────────────────────────────────────

  if (!is.null(feature_col)) {
    if (!is.character(feature_col) || length(feature_col) != 1L || !nzchar(feature_col)) {
      stop(
        "plot_oddsratio: 'feature_col' must be a single non-empty character string naming a column in 'x'.",
        call. = FALSE
      )
    }
    if (!feature_col %in% names(x)) {
      stop(
        sprintf(
          "plot_oddsratio: Column \"%s\" specified in 'feature_col' was not found in 'x'.\n  Available columns: %s",
          feature_col,
          paste(sprintf('"%s"', names(x)), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # ── 8. Build the plotting data frame ─────────────────────────────────────────

  # Resolve feature labels (special-character safe)
  if (!is.null(feature_col)) {
    feature_labels <- as.character(x[[feature_col]])
  } else if (!is.null(rownames(x)) && !all(rownames(x) == as.character(seq_len(nrow(x))))) {
    feature_labels <- rownames(x)
  } else {
    feature_labels <- as.character(seq_len(nrow(x)))
  }

  tol <- .Machine$double.eps^0.5
  or_group <- dplyr::case_when(
    is.na(or_vals)                    ~ NA_character_,
    abs(or_vals - 1) <= tol           ~ "OR = 1",
    or_vals < 1                       ~ "OR < 1",
    TRUE                              ~ "OR > 1"
  )

  plt_df <- data.frame(
    .feature  = factor(feature_labels, levels = rev(feature_labels)),
    .or       = or_vals,
    .group    = factor(or_group, levels = c("OR < 1", "OR = 1", "OR > 1")),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  if (has_ci) {
    plt_df[[".low"]]  <- x[[low]]
    plt_df[[".high"]] <- x[[high]]
  }

  # Remove rows where OR is NA
  plt_df <- plt_df[!is.na(plt_df[[".or"]]), , drop = FALSE]

  # ── 9. Compute symmetric log-scale x limits centered on 1 ───────────────────

  all_values <- plt_df[[".or"]]
  if (has_ci) {
    all_values <- c(all_values, plt_df[[".low"]], plt_df[[".high"]])
  }
  all_values <- all_values[!is.na(all_values) & all_values > 0]

  # Max log-distance from 1 on either side
  log_max <- max(abs(log(all_values)), na.rm = TRUE)
  # Add 20% padding
  log_max_padded <- log_max * 1.20
  x_lims <- exp(c(-log_max_padded, log_max_padded))

  # ── 10. Resolve axis labels ──────────────────────────────────────────────────

  ci_label <- sprintf("%g%%", ci)

  if (is.null(x_label)) {
    x_label <- if (has_ci) {
      sprintf("Odds Ratio (%s CI)", ci_label)
    } else {
      "Odds Ratio"
    }
  }
  if (is.null(y_label)) y_label <- ""

  # ── 11. Color mapping ────────────────────────────────────────────────────────

  color_map <- stats::setNames(colors, c("OR < 1", "OR = 1", "OR > 1"))

  # ── 12. Build plot ───────────────────────────────────────────────────────────

  ff <- .resolve_font(font_family)

  p <- ggplot2::ggplot(
    data    = plt_df,
    mapping = ggplot2::aes(
      x     = .data[[".or"]],
      y     = .data[[".feature"]],
      color = .data[[".group"]]
    )
  ) +
    # Reference line at OR = 1
    ggplot2::geom_vline(
      xintercept = 1,
      linetype   = "dashed",
      color      = "grey40",
      linewidth  = 0.5
    )

  # Confidence interval error bars (drawn before dots so dots sit on top)
  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      mapping  = ggplot2::aes(
        xmin = .data[[".low"]],
        xmax = .data[[".high"]]
      ),
      height    = 0.25,
      linewidth = line_size,
      na.rm     = TRUE
    )
  }

  # OR dots
  p <- p + ggplot2::geom_point(size = dot_size, na.rm = TRUE)

  # Log x-axis, symmetric limits centered on 1
  p <- p + ggplot2::scale_x_log10(
    limits = x_lims,
    expand = ggplot2::expansion(0)
  )

  # Color scale — `limits` forces all three levels into the legend even when
  # no data row falls in that group (e.g., no variable has OR exactly = 1).
  p <- p + ggplot2::scale_color_manual(
    values = color_map,
    name   = "Direction",
    limits = c("OR < 1", "OR = 1", "OR > 1"),
    labels = c(
      "OR < 1" = "OR < 1 (decreased odds)",
      "OR = 1" = "OR = 1 (no effect)",
      "OR > 1" = "OR > 1 (increased odds)"
    )
  )

  # Labels
  p <- p + ggplot2::labs(
    title    = plot_title,
    subtitle = plot_subtitle,
    x        = x_label,
    y        = y_label
  )

  # Apply theme
  p <- .apply_theme(p, theme, global_font_size, ff)

  p
}