#' Plot PCA Scree Plot
#'
#' @description
#' Creates a scree plot showing variance explained or eigenvalues per principal component.
#' Helps determine the number of meaningful PCs to retain.
#'
#' @param pca_result List. Output from `run_pca()`.
#' @param title Character. Main plot title. Default: `"Scree Plot"`.
#' @param subtitle Character. Plot subtitle. Default: `NULL` (no subtitle).
#' @param caption Character. Plot caption (bottom-right). Default: `NULL` (no caption).
#' @param position Character. Horizontal alignment of title, subtitle, and caption.
#'   Options: `"left"`, `"center"`, `"right"`. Default: `"center"`.
#' @param what Character. What to plot: `"variance"` (variance explained %) or
#'   `"eigen"` (eigenvalues). Default: `"variance"`.
#' @param type Character. Plot geometry: `"bar"`, `"line"`, or `"both"`. Default: `"line"`.
#' @param max_pc Integer. Maximum number of PCs to display. Default: `15`.
#' @param show_cumulative Logical. Show cumulative variance explained as a dashed line
#'   with a secondary y-axis. Only available when `what = "variance"`. Default: `TRUE`.
#' @param show_values Logical. Display value labels on the plot using `ggrepel` to avoid
#'   overlap. Falls back to `geom_text` if `ggrepel` is not installed. Default: `TRUE`.
#' @param bar_color Character. Fill color for bars (when `type = "bar"` or `"both"`).
#'   Default: `"#0072B2"` (Okabe-Ito blue).
#' @param line_color Character. Color for the primary line and points (when `type = "line"`
#'   or `"both"`). Default: `"#0072B2"` (Okabe-Ito blue).
#' @param cumulative_color Character. Color for the cumulative variance line, points, and
#'   labels. Default: `"#D55E00"` (Okabe-Ito vermillion).
#' @param theme Character. ggplot2 theme to apply. Options: `"nature"`, `"minimal"`,
#'   `"classic"`, `"bw"`, `"light"`, `"dark"`. Default: `"nature"` (a clean,
#'   publication-ready style based on `theme_bw()` with no gridlines).
#' @param base_size Numeric. Base font size for the theme (pts). Default: `11`.
#' @param font_family Character. Font family for all text elements. Default: `"sans"`.
#' @param axis_title_size Numeric. Font size for axis titles. Default: `base_size + 1`.
#' @param axis_text_size Numeric. Font size for axis tick labels. Default: `base_size - 1`.
#' @param plot_title_size Numeric. Font size for the plot title. Default: `base_size + 3`.
#' @param zoom Numeric. Uniform scaling factor applied to all text sizes and point/line
#'   sizes. Values > 1 enlarge, values < 1 shrink. Default: `1`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A `ggplot2` object.
#'
#' @seealso \code{\link{run_pca}}
#' 
#' @details
#' **Scree Plot Interpretation:**
#'
#' The scree plot helps determine how many principal components to retain:
#'
#' - **Elbow Method**: Look for an "elbow" where the curve flattens. PCs before the
#'   elbow capture meaningful variance; those after capture mostly noise.
#' - **Kaiser Criterion**: Retain PCs with eigenvalues > 1 (when `what = "eigen"`).
#'   These PCs explain more variance than a single original variable.
#' - **Cumulative Variance**: Retain enough PCs to explain 70–90% of total variance
#'   (shown when `show_cumulative = TRUE`).
#'
#' **Plot Types:**
#'
#' - `"bar"`: Traditional bar chart; good for discrete comparison.
#' - `"line"`: Line plot; better for visualising trends and the elbow.
#' - `"both"`: Combines bars and a line for a comprehensive view.
#'
#' **Cumulative Variance Line:**
#'
#' When `show_cumulative = TRUE`, a dashed line showing cumulative variance is
#' overlaid on a shared y-axis scale. A secondary axis on the right is labelled
#' accordingly. Only available when `what = "variance"`.
#'
#' **When to Use Each Metric:**
#'
#' - **Variance explained**: More intuitive (percentages). Recommended for most users.
#' - **Eigenvalues**: Useful for the Kaiser criterion and advanced analysis.
#'
#' **Themes:**
#'
#' The `"nature"` theme (default) applies `theme_bw()` with all gridlines removed —
#' suitable for journal figures. Other options map directly to their ggplot2 equivalents,
#' also with gridlines suppressed.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Cattell, R.B. (1966). The scree test for the number of factors.
#' *Multivariate Behavioral Research*, 1(2), 245–276.
#' \doi{10.1207/s15327906mbr0102_10}
#'
#' Kaiser, H.F. (1960). The application of electronic computers to factor analysis.
#' *Educational and Psychological Measurement*, 20(1), 141–151.
#' \doi{10.1177/001316446002000116}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_point geom_text
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous sec_axis labs theme_bw theme_minimal
#' @importFrom ggplot2 theme_classic theme_light theme_dark theme element_text element_blank
#' @importFrom ggplot2 element_line element_rect margin
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run PCA first
#' pca_result <- run_pca(scaled_data, metadata, scale_method = "auto")
#'
#' # Basic scree plot (nature theme, Okabe-Ito colors, defaults)
#' plot_scree(pca_result)
#'
#' # Eigenvalue plot, bar + line geometry
#' plot_scree(pca_result,
#'            what     = "eigen",
#'            type     = "both",
#'            title    = "PCA Eigenvalues",
#'            position = "left")
#'
#' # Variance plot without cumulative line, classic theme
#' plot_scree(pca_result,
#'            show_cumulative = FALSE,
#'            max_pc          = 10,
#'            theme           = "classic")
#'
#' # Larger font, serif family, custom zoom
#' plot_scree(pca_result,
#'            base_size   = 13,
#'            font_family = "serif",
#'            zoom        = 1.2)
#' }
plot_scree <- function(
    pca_result,
    title             = "Scree Plot",
    subtitle          = NULL,
    caption           = NULL,
    position          = "center",
    what              = "variance",
    type              = "line",
    max_pc            = 15,
    show_cumulative   = TRUE,
    show_values       = TRUE,
    bar_color         = "#0072B2",
    line_color        = "#0072B2",
    cumulative_color  = "#D55E00",
    theme             = "nature",
    base_size         = 11,
    font_family       = "sans",
    axis_title_size   = NULL,
    axis_text_size    = NULL,
    plot_title_size   = NULL,
    zoom              = 1,
    verbose           = TRUE
) {

  msg <- function(...) if (verbose) message(...)

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  #if (!inherits(pca_result, "run_pca")) {
  #  stop(
  #    "'pca_result' must be output from run_pca().\n",
  #    "Solution: Run run_pca() first on your data before plotting."
  #  )
  #}

  what          <- tolower(what)
  type          <- tolower(type)
  position      <- tolower(position)
  theme_choice  <- tolower(theme)

  if (!what %in% c("variance", "eigen")) {
    stop(
      "'what' must be 'variance' or 'eigen'.\n",
      "Solution: Use what = 'variance' (default) or what = 'eigen'."
    )
  }
  if (!type %in% c("bar", "line", "both")) {
    stop(
      "'type' must be 'bar', 'line', or 'both'.\n",
      "Solution: Use type = 'line' (default), type = 'bar', or type = 'both'."
    )
  }
  if (!position %in% c("left", "center", "right")) {
    stop(
      "'position' must be 'left', 'center', or 'right'.\n",
      "Solution: Use position = 'center' (default), 'left', or 'right'."
    )
  }

  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  if (!theme_choice %in% valid_themes) {
    stop(
      sprintf("'theme' must be one of: %s.\n", paste(valid_themes, collapse = ", ")),
      "Solution: Use theme = 'nature' (default) or another listed option."
    )
  }

  if (!is.numeric(max_pc) || length(max_pc) != 1L || max_pc < 1L) {
    stop(
      "'max_pc' must be a positive integer.\n",
      sprintf("Solution: Set max_pc to a number between 1 and %d (total PCs available).",
              pca_result$n_pcs)
    )
  }
  if (!is.character(title) || length(title) != 1L) {
    stop(
      "'title' must be a single character string.\n",
      "Solution: Provide a title like: title = 'My Scree Plot'."
    )
  }
  if (!is.null(subtitle) && (!is.character(subtitle) || length(subtitle) != 1L)) {
    stop(
      "'subtitle' must be NULL or a single character string.\n",
      "Solution: Either omit subtitle or provide: subtitle = 'My subtitle'."
    )
  }
  if (!is.null(caption) && (!is.character(caption) || length(caption) != 1L)) {
    stop(
      "'caption' must be NULL or a single character string.\n",
      "Solution: Either omit caption or provide: caption = 'My caption'."
    )
  }
  if (!is.numeric(base_size) || length(base_size) != 1L || base_size <= 0) {
    stop(
      "'base_size' must be a positive numeric value.\n",
      "Solution: Use base_size = 11 (default) or another positive number."
    )
  }
  if (!is.character(font_family) || length(font_family) != 1L) {
    stop(
      "'font_family' must be a single character string.\n",
      "Solution: Use font_family = 'sans' (default), 'serif', or 'mono'."
    )
  }
  if (!is.numeric(zoom) || length(zoom) != 1L || zoom <= 0) {
    stop(
      "'zoom' must be a single positive numeric value.\n",
      "Solution: Use zoom = 1 (default), zoom = 1.5 to enlarge, or zoom = 0.8 to shrink."
    )
  }

  if (show_cumulative && what != "variance") {
    warning(
      "'show_cumulative' is only available when what = 'variance'. Setting show_cumulative = FALSE.\n",
      "Solution: Use what = 'variance' if you want to show cumulative variance."
    )
    show_cumulative <- FALSE
  }

  # ggrepel availability
  use_ggrepel <- if (show_values && !requireNamespace("ggrepel", quietly = TRUE)) {
    warning(
      "Package 'ggrepel' not available. Using geom_text instead (labels may overlap).\n",
      "Solution: Install it with: install.packages('ggrepel')"
    )
    FALSE
  } else {
    show_values
  }

  # ============================================================================
  # DATA PREPARATION
  # ============================================================================

  n_pcs    <- pca_result$n_pcs
  max_show <- min(as.integer(max_pc), n_pcs)

  if (max_pc > n_pcs) {
    warning(
      sprintf("max_pc = %d requested but only %d PCs available. Showing %d PCs.\n",
              max_pc, n_pcs, n_pcs),
      sprintf("Solution: Reduce max_pc to a value <= %d.", n_pcs)
    )
  }

  if (what == "variance") {
    y_values <- pca_result$variance_explained[seq_len(max_show)]
    y_label  <- "Variance Explained (%)"
  } else {
    y_values <- pca_result$eigenvalues[seq_len(max_show)]
    y_label  <- "Eigenvalue"
  }

  scree_data <- data.frame(
    PC_num   = seq_len(max_show),
    Value    = y_values,
    CumValue = if (what == "variance") pca_result$cumulative_variance[seq_len(max_show)] else NA_real_,
    stringsAsFactors = FALSE
  )

  msg(sprintf("Creating scree plot (%s, type: %s, max PC: %d)", what, type, max_show))

  # ============================================================================
  # RESOLVE SIZE DEFAULTS & ZOOM
  # ============================================================================

  axis_title_size <- if (is.null(axis_title_size)) base_size + 1 else axis_title_size
  axis_text_size  <- if (is.null(axis_text_size))  base_size - 1 else axis_text_size
  plot_title_size <- if (is.null(plot_title_size)) base_size + 3 else plot_title_size

  .z <- function(x) x * zoom

  # ============================================================================
  # THEME CONSTRUCTION
  # ============================================================================

  # All themes have gridlines removed via custom_theme below
  base_theme <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
      ggplot2::theme(
        panel.border     = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.6),
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold", family = font_family)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size = base_size,      base_family = font_family),
    "light"   = ggplot2::theme_light(base_size = base_size,   base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size = base_size,    base_family = font_family)
  )

  hjust_val <- switch(position, "left" = 0, "center" = 0.5, "right" = 1, 0.5)

  custom_theme <- ggplot2::theme(
    panel.grid.major      = ggplot2::element_blank(),
    panel.grid.minor      = ggplot2::element_blank(),
    plot.title            = ggplot2::element_text(
      hjust = hjust_val, size = .z(plot_title_size), face = "bold", family = font_family
    ),
    plot.subtitle         = ggplot2::element_text(
      hjust = hjust_val, size = .z(plot_title_size - 2), family = font_family
    ),
    plot.caption          = ggplot2::element_text(
      hjust = 1, size = .z(base_size - 2), color = "grey50", family = font_family
    ),
    axis.title.x          = ggplot2::element_text(size = .z(axis_title_size), family = font_family),
    axis.title.y          = ggplot2::element_text(size = .z(axis_title_size), family = font_family),
    axis.title.y.right    = if (show_cumulative && what == "variance") {
      ggplot2::element_text(
        size = .z(axis_title_size), color = cumulative_color, family = font_family
      )
    } else {
      ggplot2::element_blank()
    },
    axis.text.x           = ggplot2::element_text(size = .z(axis_text_size), family = font_family),
    axis.text.y           = if (show_values && !show_cumulative) {
      ggplot2::element_blank()
    } else {
      ggplot2::element_text(size = .z(axis_text_size), family = font_family)
    },
    axis.text.y.right     = if (show_cumulative && what == "variance") {
      ggplot2::element_text(
        size = .z(axis_text_size), color = cumulative_color, family = font_family
      )
    } else {
      ggplot2::element_blank()
    },
    axis.ticks.y          = if (show_values && !show_cumulative) {
      ggplot2::element_blank()
    } else {
      ggplot2::element_line()
    },
    plot.margin           = ggplot2::margin(10, 10, 10, 10)
  )

  # ============================================================================
  # BUILD PLOT
  # ============================================================================

  p <- ggplot2::ggplot(scree_data, ggplot2::aes(x = PC_num, y = Value))

  # Primary geometry
  if (type == "bar") {
    p <- p + ggplot2::geom_col(fill = bar_color, alpha = 0.75, width = 0.6)
  } else if (type == "line") {
    p <- p +
      ggplot2::geom_line(color = line_color, linewidth = .z(1.2), group = 1L) +
      ggplot2::geom_point(color = line_color, size = .z(3))
  } else {  # both
    p <- p +
      ggplot2::geom_col(fill = bar_color, alpha = 0.5, width = 0.6) +
      ggplot2::geom_line(color = line_color, linewidth = .z(1.2), group = 1L) +
      ggplot2::geom_point(color = line_color, size = .z(3))
  }

  # Primary value labels
  if (show_values) {
    fmt      <- if (what == "variance") "%.1f" else "%.2f"
    lbl_data <- data.frame(
      PC_num = scree_data$PC_num,
      Value  = scree_data$Value,
      Label  = sprintf(fmt, scree_data$Value)
    )
    nudge_y_val <- max(scree_data$Value, na.rm = TRUE) * 0.05

    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data               = lbl_data,
        ggplot2::aes(label = Label),
        size               = .z(base_size / 3.5),
        color              = "black",
        family             = font_family,
        box.padding        = 0.3,
        point.padding      = 0.2,
        segment.color      = "grey50",
        segment.linewidth  = 0.3,
        min.segment.length = 0,
        max.overlaps       = Inf,
        direction          = "y",
        nudge_y            = nudge_y_val
      )
    } else {
      p <- p + ggplot2::geom_text(
        data               = lbl_data,
        ggplot2::aes(label = Label),
        vjust              = -0.5,
        size               = .z(base_size / 3.5),
        color              = "black",
        family             = font_family
      )
    }
  }

  # Cumulative variance overlay
  if (show_cumulative && what == "variance") {
    cum_data <- data.frame(
      PC_num   = scree_data$PC_num,
      CumValue = scree_data$CumValue
    )
    nudge_cum <- max(cum_data$CumValue, na.rm = TRUE) * 0.02

    p <- p +
      ggplot2::geom_line(
        data               = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue),
        color              = cumulative_color,
        linewidth          = .z(1),
        linetype           = "dashed",
        group              = 1L
      ) +
      ggplot2::geom_point(
        data               = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue),
        color              = cumulative_color,
        size               = .z(2.5)
      )

    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data               = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue,
                     label = sprintf("%.1f", CumValue)),
        size               = .z(base_size / 4),
        color              = cumulative_color,
        family             = font_family,
        fontface           = "italic",
        box.padding        = 0.3,
        point.padding      = 0.2,
        segment.color      = "grey50",
        segment.linewidth  = 0.2,
        min.segment.length = 0,
        max.overlaps       = Inf,
        direction          = "y",
        nudge_y            = nudge_cum
      )
    } else {
      p <- p + ggplot2::geom_text(
        data               = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue,
                     label = sprintf("%.1f", CumValue)),
        vjust              = -0.5,
        size               = .z(base_size / 4),
        color              = cumulative_color,
        fontface           = "italic",
        family             = font_family
      )
    }

    p <- p +
      ggplot2::scale_y_continuous(
        name     = y_label,
        sec.axis = ggplot2::sec_axis(~ ., name = "Cumulative Variance (%)")
      )
  }

  # Axes, labels, theme
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(max_show),
      labels = as.character(seq_len(max_show))
    ) +
    ggplot2::labs(
      x        = "Principal Component",
      y        = if (!show_cumulative || what != "variance") y_label else NULL,
      title    = title,
      subtitle = subtitle,
      caption  = caption
    ) +
    base_theme +
    custom_theme

  msg("Scree plot created successfully.")
  return(p)
}