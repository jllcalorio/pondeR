# ============================================================
#  plot_multi_y.R
#  Part of the pondeR package
# ============================================================

# ---- Internal helpers ------------------------------------------------------

#' @keywords internal
.check_zoom <- function(zoom) {
  if (!is.numeric(zoom) || length(zoom) != 1L || is.na(zoom) || zoom <= 0) {
    stop(
      "`zoom` must be a single positive number (e.g. zoom = 1.5). ",
      "You supplied: ", deparse(zoom),
      call. = FALSE
    )
  }
  zoom
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @keywords internal
.validate_col <- function(df, col, arg_name) {
  if (!is.character(col) || length(col) != 1L) {
    stop("`", arg_name, "` must be a single character string (a column name). ",
         "You supplied an object of class: ", paste(class(col), collapse = "/"),
         call. = FALSE)
  }
  if (!col %in% names(df)) {
    stop("`", arg_name, "` = \"", col, "\" was not found in the data. ",
         "Available columns: ", paste(names(df), collapse = ", "),
         call. = FALSE)
  }
}

#' @keywords internal
.validate_cols <- function(df, cols, arg_name) {
  if (is.null(cols)) return(invisible(NULL))
  if (!is.character(cols)) {
    stop("`", arg_name, "` must be a character vector of column names. ",
         "You supplied an object of class: ", paste(class(cols), collapse = "/"),
         call. = FALSE)
  }
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("`", arg_name, "` contains column(s) not found in the data: ",
         paste(missing_cols, collapse = ", "), ". ",
         "Available columns: ", paste(names(df), collapse = ", "),
         call. = FALSE)
  }
}

#' @keywords internal
.validate_numeric_cols <- function(df, cols, arg_name) {
  if (is.null(cols)) return(invisible(NULL))
  non_num <- cols[!vapply(df[cols], is.numeric, logical(1L))]
  if (length(non_num) > 0L) {
    stop("`", arg_name, "` column(s) must be numeric, but the following are not: ",
         paste(non_num, collapse = ", "),
         call. = FALSE)
  }
}

#' @keywords internal
.validate_order <- function(order_vec, arg_name) {
  valid <- c("hist", "area", "line", "dot")
  bad   <- setdiff(order_vec, valid)
  if (length(bad) > 0L) {
    stop("`", arg_name, "` contains invalid values: ", paste(bad, collapse = ", "),
         ". Only combinations of c(\"hist\", \"area\", \"line\", \"dot\") are allowed.",
         call. = FALSE)
  }
  if (any(duplicated(order_vec))) {
    stop("`", arg_name, "` must not contain duplicated values.",
         call. = FALSE)
  }
}

# Okabe-Ito + extended colorblind-safe palette
.oi_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000",
  "#88CCEE", "#CC6677", "#DDCC77", "#117733",
  "#332288", "#AA4499"
)

#' @keywords internal
.resolve_colors <- function(user_colors, n, arg_name) {
  if (n == 0L) return(character(0L))
  if (is.null(user_colors) || length(user_colors) == 0L) {
    return(rep_len(.oi_palette, n))
  }
  if (!is.character(user_colors)) {
    stop("`", arg_name, "` must be a character vector of color names or hex codes.",
         call. = FALSE)
  }
  rep_len(user_colors, n)
}

# ---- Theme builder ---------------------------------------------------------

#' @keywords internal
.ponder_theme_multi_y <- function(
    theme_name    = "nature",
    base_size     = 20,
    fs_title      = NULL,
    fs_subtitle   = NULL,
    fs_caption    = NULL,
    fs_axis       = NULL,
    fs_legend     = NULL,
    fs_facet      = NULL,
    legend_pos    = "top",
    zoom          = 1
) {

  fs_title    <- (fs_title    %||% (base_size * 1.20)) * zoom
  fs_subtitle <- (fs_subtitle %||% (base_size * 0.90)) * zoom
  fs_caption  <- (fs_caption  %||% (base_size * 0.75)) * zoom
  fs_axis     <- (fs_axis     %||% (base_size * 1.00)) * zoom
  fs_legend   <- (fs_legend   %||% (base_size * 0.90)) * zoom
  fs_facet    <- (fs_facet    %||% (base_size * 0.85)) * zoom

  base_t <- ggplot2::theme_bw(base_size = base_size * zoom)

  nature_add <- ggplot2::theme(
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.border       = ggplot2::element_rect(colour = "grey80", fill = NA,
                                               linewidth = 0.5 * zoom),
    axis.ticks         = ggplot2::element_line(colour = "grey60",
                                               linewidth = 0.4 * zoom)
  )

  theme_add <- switch(
    theme_name,
    "nature"  = nature_add,
    "minimal" = ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "grey60",
                                               linewidth = 0.4 * zoom)
    ),
    "classic" = ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "black",
                                               linewidth = 0.5 * zoom)
    ),
    {
      warning("`theme` = \"", theme_name, "\" is not recognised. ",
              "Falling back to \"nature\". ",
              "Valid options: \"nature\", \"minimal\", \"classic\".",
              call. = FALSE)
      nature_add
    }
  )

  lp <- switch(
    legend_pos,
    "top"    = "top",
    "bottom" = "bottom",
    "left"   = "left",
    "right"  = "right",
    {
      warning("`legend_pos` = \"", legend_pos, "\" is not recognised. ",
              "Falling back to \"top\".",
              call. = FALSE)
      "top"
    }
  )

  base_t +
    theme_add +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(size = fs_title,
                                                face = "bold",
                                                hjust = 0),
      plot.subtitle     = ggplot2::element_text(size = fs_subtitle,
                                                hjust = 0,
                                                colour = "grey40"),
      plot.caption      = ggplot2::element_text(size = fs_caption,
                                                hjust = 1,
                                                colour = "grey50"),
      axis.title        = ggplot2::element_text(size = fs_axis),
      axis.text         = ggplot2::element_text(size = fs_axis * 0.9),
      axis.title.y.right = ggplot2::element_text(size = fs_axis,
                                                  angle = 90,
                                                  vjust = 0.5),
      legend.position   = lp,
      legend.text       = ggplot2::element_text(size = fs_legend),
      legend.title      = ggplot2::element_blank(),
      legend.key.size   = ggplot2::unit(fs_legend * 0.08, "cm"),
      strip.text        = ggplot2::element_text(size = fs_facet,
                                                face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "grey95",
                                                colour = "grey80")
    )
}

# ---- Main scaling helper for dual axis -------------------------------------

#' @keywords internal
.scale_to_primary <- function(primary_range, secondary_range) {
  slope     <- diff(primary_range)  / diff(secondary_range)
  intercept <- primary_range[1L]    - slope * secondary_range[1L]
  list(
    transform = function(y) slope * y + intercept,
    inverse   = function(y) (y - intercept) / slope,
    slope     = slope,
    intercept = intercept
  )
}

# ---- Long-format pivot helper ----------------------------------------------

#' @keywords internal
# Pivots selected columns to long format, keeping x_axis (and optionally
# facet) intact.  Rows where x_axis is NA are always dropped.  Rows where
# the y value is NA are kept as NA rows so ggplot2 can render partial series
# without dropping entire x positions.
.pivot_long <- function(df, cols, x_axis, facet = NULL, transform_fn = NULL) {
  if (is.null(cols) || length(cols) == 0L) return(NULL)

  keep_cols <- unique(c(x_axis, facet))

  # Drop rows where x_axis itself is NA — nothing to plot there
  df <- df[!is.na(df[[x_axis]]), , drop = FALSE]

  out <- do.call(rbind, lapply(cols, function(col) {
    sub <- df[, keep_cols, drop = FALSE]
    sub$.series <- col
    sub$.y      <- df[[col]]
    if (!is.null(transform_fn)) {
      sub$.y <- transform_fn(sub$.y)
    }
    sub
  }))

  out$.series <- factor(out$.series, levels = cols)
  out
}

# ============================================================
#  plot_multi_y  (main export)
# ============================================================

#' @title Plot Multiple Variables on Primary and Secondary Y Axes
#'
#' @description
#' Creates a dual-axis ggplot with a shared x-axis, allowing simultaneous
#' display of multiple variables on independent primary and secondary y scales.
#' Supports four geometry types (line, area, histogram, dot) on each axis,
#' optional faceting, colorblind-safe palettes, and full publication-ready
#' typography controls.
#'
#' @param x A \code{data.frame}, \code{tibble}, or \code{matrix} (converted
#'   internally). Column names may contain special characters.
#' @param x_axis A single character string giving the column of \code{x} to
#'   use on the x-axis. May be numeric, Date/POSIXct, or categorical.
#' @param main_as_line Character vector of column name(s) to plot as line
#'   geometry on the primary y-axis. \code{NULL} omits this geometry.
#' @param main_as_area Character vector of column name(s) to plot as area
#'   geometry on the primary y-axis. \code{NULL} omits this geometry.
#' @param main_as_hist Character vector of column name(s) to plot as histogram
#'   (bar) geometry on the primary y-axis. \code{NULL} omits this geometry.
#' @param main_as_dot Character vector of column name(s) to plot as point
#'   geometry on the primary y-axis. \code{NULL} omits this geometry.
#' @param secondary_as_line Character vector of column name(s) to plot as line
#'   geometry on the secondary y-axis. \code{NULL} omits this geometry.
#' @param secondary_as_area Character vector of column name(s) to plot as area
#'   geometry on the secondary y-axis. \code{NULL} omits this geometry.
#' @param secondary_as_hist Character vector of column name(s) to plot as
#'   histogram (bar) geometry on the secondary y-axis. \code{NULL} omits this
#'   geometry.
#' @param secondary_as_dot Character vector of column name(s) to plot as point
#'   geometry on the secondary y-axis. \code{NULL} omits this geometry.
#' @param order_main Character vector specifying the back-to-front layer order
#'   of primary geometries. Default \code{c("hist", "area", "line", "dot")}.
#' @param order_secondary Same as \code{order_main} but for secondary
#'   geometries.
#' @param color_main_line Character vector of R color names or hex codes for
#'   primary line geometry. Recycled if fewer colors than series.
#'   Defaults to the Okabe-Ito colorblind-safe palette.
#' @param color_main_area Character vector of colors for primary area geometry.
#' @param color_main_hist Character vector of colors for primary histogram
#'   geometry.
#' @param color_main_dot Character vector of colors for primary dot geometry.
#' @param color_secondary_line Character vector of colors for secondary line
#'   geometry.
#' @param color_secondary_area Character vector of colors for secondary area
#'   geometry.
#' @param color_secondary_hist Character vector of colors for secondary
#'   histogram geometry.
#' @param color_secondary_dot Character vector of colors for secondary dot
#'   geometry.
#' @param linewidth_main_line Numeric. Line width for primary line geometry.
#'   Default \code{0.8}.
#' @param linewidth_secondary_line Numeric. Line width for secondary line
#'   geometry. Default \code{0.8}.
#' @param size_main_dot Numeric. Point size for primary dot geometry.
#'   Default \code{2}.
#' @param size_secondary_dot Numeric. Point size for secondary dot geometry.
#'   Default \code{2}.
#' @param area_alpha_main Numeric in \eqn{[0, 1]}. Fill transparency for
#'   primary area and histogram geometries. Default \code{0.35}.
#' @param area_alpha_secondary Numeric in \eqn{[0, 1]}. Fill transparency for
#'   secondary area and histogram geometries. Default \code{0.35}.
#' @param facet A single character string naming a categorical column of
#'   \code{x} to use for \code{\link[ggplot2]{facet_wrap}}. \code{NULL}
#'   (default) disables faceting.
#' @param facet_ncol Integer. Number of columns when \code{facet} is set.
#'   \code{NULL} lets ggplot2 choose automatically.
#' @param legend_pos Character. Legend position: \code{"top"} (default),
#'   \code{"bottom"}, \code{"left"}, or \code{"right"}.
#' @param title Character. Plot title. \code{NULL} omits the title.
#' @param subtitle Character. Plot subtitle. \code{NULL} omits the subtitle.
#' @param caption Character. Plot caption. \code{NULL} omits the caption.
#' @param x_lab Character. x-axis label. Defaults to the value of
#'   \code{x_axis}.
#' @param y_lab_main Character. Primary y-axis label. Defaults to
#'   \code{"Primary axis"} when secondary variables are supplied, otherwise
#'   an empty string.
#' @param y_lab_secondary Character. Secondary y-axis label. Defaults to
#'   \code{"Secondary axis"}.
#' @param global_font_size Numeric. Base font size in points. Default
#'   \code{20}.
#' @param title_font_size Numeric. Override font size for the plot title.
#'   Defaults to \code{1.20 * global_font_size}.
#' @param subtitle_font_size Numeric. Override font size for the subtitle.
#'   Defaults to \code{0.90 * global_font_size}.
#' @param caption_font_size Numeric. Override font size for the caption.
#'   Defaults to \code{0.75 * global_font_size}.
#' @param axis_font_size Numeric. Override font size for axis titles and text.
#'   Defaults to \code{1.00 * global_font_size}.
#' @param legend_font_size Numeric. Override font size for legend text.
#'   Defaults to \code{0.90 * global_font_size}.
#' @param facet_font_size Numeric. Override font size for facet strip labels.
#'   Only relevant when \code{facet} is non-\code{NULL}. Defaults to
#'   \code{0.85 * global_font_size}.
#' @param theme Character. Visual theme. One of \code{"nature"} (default),
#'   \code{"minimal"}, or \code{"classic"}.
#' @param zoom Numeric. Multiplicative scaling factor applied to all resolved
#'   sizes at construction time. Default \code{1}.
#'
#' @details
#' \strong{Dual-axis implementation.}
#' ggplot2 implements secondary axes as linear rescalings of the primary axis
#' (\code{\link[ggplot2]{sec_axis}}). \code{plot_multi_y} computes the linear
#' transform mapping the secondary variables' observed range into the primary
#' range with a 5\% margin on each side, then back-transforms the tick labels.
#'
#' \strong{NA handling.}
#' Rows where \code{x_axis} is \code{NA} are always dropped before plotting.
#' Rows where a y-variable is \code{NA} are retained so that other series at
#' the same x position are still drawn; ggplot2 renders gaps in those series
#' naturally.  A warning is issued listing which columns contain \code{NA}
#' values so you remain aware of the gaps.
#'
#' \strong{Long-format pivoting.}
#' Each geometry group is pivoted to long format internally, with a
#' \code{.series} factor carrying the original column name. This is what
#' feeds the colour/fill aesthetics and drives the legend, ensuring the
#' manual colour scale always has matching levels.
#'
#' \strong{Layer ordering.}
#' \code{order_main} and \code{order_secondary} control the painter's order.
#' The default \code{c("hist", "area", "line", "dot")} places histograms at
#' the back, then areas, then lines, with dots in the foreground.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' \dontrun{
#' # --- Basic: area on primary, two lines on primary, faceted by Year ---
#' set.seed(42)
#' df <- data.frame(
#'   Year         = factor(rep(c("2021", "2022"), each = 20)),
#'   `Week No.`   = rep(1:20, 2),
#'   `Weekly cases` = c(rpois(20, 60), rpois(20, 80)),
#'   `RdRP GC/L`  = c(rep(NA, 5), runif(15, 1e4, 1e5),
#'                    rep(NA, 3), runif(17, 2e4, 2e5)),
#'   `N GC/L`     = c(rep(NA, 8), runif(12, 5e3, 8e4),
#'                    rep(NA, 4), runif(16, 1e4, 1e5)),
#'   check.names  = FALSE
#' )
#'
#' plot_multi_y(
#'   x             = df,
#'   x_axis        = "Week No.",
#'   main_as_area  = "Weekly cases",
#'   main_as_line  = c("RdRP GC/L", "N GC/L"),
#'   facet         = "Year",
#'   title         = NULL
#' )
#'
#' # --- Dual axis: cases on primary, viral load on secondary ---
#' plot_multi_y(
#'   x                  = df,
#'   x_axis             = "Week No.",
#'   main_as_area       = "Weekly cases",
#'   secondary_as_line  = c("RdRP GC/L", "N GC/L"),
#'   y_lab_main         = "Weekly cases",
#'   y_lab_secondary    = "GC/L",
#'   facet              = "Year"
#' )
#' }
#'
#' @author John Lennon L. Calorio
#' @export
plot_multi_y <- function(
    x,
    x_axis,
    main_as_line           = NULL,
    main_as_area           = NULL,
    main_as_hist           = NULL,
    main_as_dot            = NULL,
    secondary_as_line      = NULL,
    secondary_as_area      = NULL,
    secondary_as_hist      = NULL,
    secondary_as_dot       = NULL,
    order_main             = c("hist", "area", "line", "dot"),
    order_secondary        = c("hist", "area", "line", "dot"),
    color_main_line        = NULL,
    color_main_area        = NULL,
    color_main_hist        = NULL,
    color_main_dot         = NULL,
    color_secondary_line   = NULL,
    color_secondary_area   = NULL,
    color_secondary_hist   = NULL,
    color_secondary_dot    = NULL,
    linewidth_main_line      = 0.8,
    linewidth_secondary_line = 0.8,
    size_main_dot            = 2,
    size_secondary_dot       = 2,
    area_alpha_main          = 0.35,
    area_alpha_secondary     = 0.35,
    facet                    = NULL,
    facet_ncol               = NULL,
    legend_pos               = "top",
    title                    = NULL,
    subtitle                 = NULL,
    caption                  = NULL,
    x_lab                    = NULL,
    y_lab_main               = NULL,
    y_lab_secondary          = "Secondary axis",
    global_font_size         = 20,
    title_font_size          = NULL,
    subtitle_font_size       = NULL,
    caption_font_size        = NULL,
    axis_font_size           = NULL,
    legend_font_size         = NULL,
    facet_font_size          = NULL,
    theme                    = "nature",
    zoom                     = 1
) {

  # ---- 0. Dependencies -----------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("`ggplot2` is required. Install with: install.packages(\"ggplot2\")",
         call. = FALSE)

  # ---- 1. Coerce input -----------------------------------------------------
  if (inherits(x, "matrix")) x <- as.data.frame(x)
  if (!is.data.frame(x))
    stop("`x` must be a data.frame, tibble, or matrix. ",
         "You supplied: ", paste(class(x), collapse = "/"), call. = FALSE)
  if (nrow(x) == 0L) stop("`x` has zero rows. Nothing to plot.", call. = FALSE)
  if (ncol(x) == 0L) stop("`x` has zero columns.",              call. = FALSE)

  # ---- 2. Validate scalars -------------------------------------------------
  zoom <- .check_zoom(zoom)

  for (nm in c("global_font_size", "linewidth_main_line",
               "linewidth_secondary_line", "size_main_dot",
               "size_secondary_dot")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0)
      stop("`", nm, "` must be a single positive number.", call. = FALSE)
  }
  for (nm in c("area_alpha_main", "area_alpha_secondary")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1)
      stop("`", nm, "` must be a single number in [0, 1].", call. = FALSE)
  }
  for (nm in c("title_font_size", "subtitle_font_size", "caption_font_size",
               "axis_font_size", "legend_font_size", "facet_font_size")) {
    v <- get(nm)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0))
      stop("`", nm, "` must be a single positive number or NULL.", call. = FALSE)
  }

  # Apply zoom to sizes at construction time
  lw_main <- linewidth_main_line      * zoom
  lw_sec  <- linewidth_secondary_line * zoom
  sz_main <- size_main_dot            * zoom
  sz_sec  <- size_secondary_dot       * zoom

  # ---- 3. Validate column args ---------------------------------------------
  .validate_col(x, x_axis, "x_axis")

  all_main_cols <- c(main_as_line, main_as_area, main_as_hist, main_as_dot)
  all_sec_cols  <- c(secondary_as_line, secondary_as_area,
                     secondary_as_hist, secondary_as_dot)

  if (length(all_main_cols) == 0L && length(all_sec_cols) == 0L)
    stop("At least one of `main_as_*` or `secondary_as_*` must be non-NULL.",
         call. = FALSE)

  for (nm in c("main_as_line","main_as_area","main_as_hist","main_as_dot",
               "secondary_as_line","secondary_as_area",
               "secondary_as_hist","secondary_as_dot")) {
    .validate_cols(x, get(nm), nm)
  }
  .validate_numeric_cols(x, all_main_cols, "main_as_*")
  .validate_numeric_cols(x, all_sec_cols,  "secondary_as_*")

  if (!is.null(facet)) {
    .validate_col(x, facet, "facet")
    if (is.numeric(x[[facet]])) {
      warning("`facet` column \"", facet, "\" is numeric and will be coerced ",
              "to a factor. Convert it beforehand if this is unintended.",
              call. = FALSE)
      x[[facet]] <- factor(x[[facet]])
    }
  }

  .validate_order(order_main,      "order_main")
  .validate_order(order_secondary, "order_secondary")

  # ---- 4. Warn about NA y values (informational only) ----------------------
  y_cols_all <- c(all_main_cols, all_sec_cols)
  na_cols <- y_cols_all[vapply(y_cols_all,
                               function(col) anyNA(x[[col]]), logical(1L))]
  if (length(na_cols) > 0L) {
    warning(
      "The following column(s) contain NA values: ",
      paste(na_cols, collapse = ", "), ". ",
      "Rows where x_axis (\"", x_axis, "\") is non-missing are retained; ",
      "gaps will appear in the affected series at those x positions.",
      call. = FALSE
    )
  }

  # ---- 5. Drop rows where x_axis is NA ------------------------------------
  x_na <- is.na(x[[x_axis]])
  if (any(x_na)) {
    warning(sum(x_na), " row(s) where `x_axis` = \"", x_axis, "\" is NA ",
            "have been removed.", call. = FALSE)
    x <- x[!x_na, , drop = FALSE]
  }

  has_secondary <- length(all_sec_cols) > 0L

  # ---- 6. Dual-axis transform ----------------------------------------------
  tf <- NULL
  if (has_secondary) {
    sec_finite <- unlist(lapply(all_sec_cols,
                                function(col) x[[col]][is.finite(x[[col]])]))
    if (length(sec_finite) == 0L)
      stop("All secondary y columns contain only NA/Inf. ",
           "Cannot compute secondary axis range.", call. = FALSE)

    main_finite <- if (length(all_main_cols) > 0L) {
      unlist(lapply(all_main_cols,
                    function(col) x[[col]][is.finite(x[[col]])]))
    } else sec_finite

    if (length(main_finite) == 0L) main_finite <- sec_finite

    pad <- function(r, frac = 0.05) {
      span <- diff(r); if (span == 0) span <- abs(r[1L]) * 0.1 + 1
      c(r[1L] - frac * span, r[2L] + frac * span)
    }

    pr <- pad(range(main_finite, na.rm = TRUE))
    sr <- pad(range(sec_finite,  na.rm = TRUE))
    tf <- .scale_to_primary(pr, sr)
  }

  # ---- 7. Resolve color palettes -------------------------------------------
  # Colors are assigned per-column within each geom type, then merged into a
  # single named vector that exactly matches the .series factor levels used in
  # the long-format data.  This guarantees scale_colour/fill_manual never sees
  # a mismatch between names(values) and the data's colour values.

  color_lookup <- c(
    stats::setNames(.resolve_colors(color_main_line,        length(main_as_line),        "color_main_line"),        main_as_line),
    stats::setNames(.resolve_colors(color_main_area,        length(main_as_area),        "color_main_area"),        main_as_area),
    stats::setNames(.resolve_colors(color_main_hist,        length(main_as_hist),        "color_main_hist"),        main_as_hist),
    stats::setNames(.resolve_colors(color_main_dot,         length(main_as_dot),         "color_main_dot"),         main_as_dot),
    stats::setNames(.resolve_colors(color_secondary_line,   length(secondary_as_line),   "color_secondary_line"),   secondary_as_line),
    stats::setNames(.resolve_colors(color_secondary_area,   length(secondary_as_area),   "color_secondary_area"),   secondary_as_area),
    stats::setNames(.resolve_colors(color_secondary_hist,   length(secondary_as_hist),   "color_secondary_hist"),   secondary_as_hist),
    stats::setNames(.resolve_colors(color_secondary_dot,    length(secondary_as_dot),    "color_secondary_dot"),    secondary_as_dot)
  )
  # In case the same column appears in multiple geom types, last write wins
  color_lookup <- color_lookup[!duplicated(names(color_lookup), fromLast = TRUE)]

  # ---- 8. Long-format data builders ----------------------------------------

  make_long <- function(cols, transform_fn = NULL) {
    .pivot_long(x, cols, x_axis, facet, transform_fn)
  }

  ld_ml <- make_long(main_as_line)
  ld_ma <- make_long(main_as_area)
  ld_mh <- make_long(main_as_hist)
  ld_md <- make_long(main_as_dot)
  ld_sl <- make_long(secondary_as_line,  if (!is.null(tf)) tf$transform else NULL)
  ld_sa <- make_long(secondary_as_area,  if (!is.null(tf)) tf$transform else NULL)
  ld_sh <- make_long(secondary_as_hist,  if (!is.null(tf)) tf$transform else NULL)
  ld_sd <- make_long(secondary_as_dot,   if (!is.null(tf)) tf$transform else NULL)

  # ---- 9. Geom layer factory -----------------------------------------------
  # Each geom uses data = <long df>, aes(x, y = .y, colour/fill = .series).
  # This ensures the colour scale levels always match the data.

  make_geom <- function(long_df, geom_type, lw = NULL, sz = NULL, alpha = 0.35) {
    if (is.null(long_df)) return(NULL)

    base_aes <- ggplot2::aes(
      x      = .data[[x_axis]],
      y      = .data[[".y"]],
      colour = .data[[".series"]],
      fill   = .data[[".series"]],
      group  = .data[[".series"]]
    )

    switch(
      geom_type,
      "line" = ggplot2::geom_line(
        mapping   = base_aes,
        data      = long_df,
        linewidth = lw,
        na.rm     = FALSE,   # keep NAs → natural gaps in the line
        inherit.aes = FALSE
      ),
      "area" = ggplot2::geom_area(
        mapping   = base_aes,
        data      = long_df,
        alpha     = alpha,
        colour    = NA,
        na.rm     = FALSE,
        inherit.aes = FALSE
      ),
      "hist" = ggplot2::geom_col(
        mapping   = base_aes,
        data      = long_df,
        alpha     = alpha,
        colour    = NA,
        na.rm     = FALSE,
        inherit.aes = FALSE
      ),
      "dot"  = ggplot2::geom_point(
        mapping   = base_aes,
        data      = long_df,
        size      = sz,
        na.rm     = FALSE,
        inherit.aes = FALSE
      ),
      stop("Unknown geom_type: ", geom_type, call. = FALSE)
    )
  }

  # ---- 10. Ordered layer assembly ------------------------------------------

  layer_map <- list(
    main = list(
      line = list(df = ld_ml, lw = lw_main, sz = NULL,    alpha = area_alpha_main),
      area = list(df = ld_ma, lw = NULL,    sz = NULL,    alpha = area_alpha_main),
      hist = list(df = ld_mh, lw = NULL,    sz = NULL,    alpha = area_alpha_main),
      dot  = list(df = ld_md, lw = NULL,    sz = sz_main, alpha = area_alpha_main)
    ),
    secondary = list(
      line = list(df = ld_sl, lw = lw_sec,  sz = NULL,   alpha = area_alpha_secondary),
      area = list(df = ld_sa, lw = NULL,    sz = NULL,   alpha = area_alpha_secondary),
      hist = list(df = ld_sh, lw = NULL,    sz = NULL,   alpha = area_alpha_secondary),
      dot  = list(df = ld_sd, lw = NULL,    sz = sz_sec, alpha = area_alpha_secondary)
    )
  )

  collect_layers <- function(order_vec, axis) {
    layers <- list()
    for (gtype in order_vec) {
      m <- layer_map[[axis]][[gtype]]
      if (!is.null(m$df)) {
        layers <- c(layers, list(
          make_geom(m$df, gtype, lw = m$lw, sz = m$sz, alpha = m$alpha)
        ))
      }
    }
    layers
  }

  main_layers <- collect_layers(order_main,      "main")
  sec_layers  <- collect_layers(order_secondary, "secondary")

  # ---- 11. Initialise base plot and add layers -----------------------------
  p <- ggplot2::ggplot()

  for (lyr in main_layers) p <- p + lyr
  for (lyr in sec_layers)  p <- p + lyr

  # ---- 12. Colour / fill scales -------------------------------------------
  # color_lookup names == .series factor levels → no mismatch possible
  if (length(color_lookup) > 0L) {
    p <- p +
      ggplot2::scale_colour_manual(values = color_lookup, name = NULL) +
      ggplot2::scale_fill_manual(  values = color_lookup, name = NULL)
  }

  # ---- 13. Y-axis / secondary axis ----------------------------------------
  if (has_secondary && !is.null(tf)) {
    p <- p + ggplot2::scale_y_continuous(
      name     = y_lab_main %||% "Primary axis",
      sec.axis = ggplot2::sec_axis(
        transform = ~ tf$inverse(.),
        name      = y_lab_secondary,
        labels    = function(b) {
          mx <- max(abs(b), na.rm = TRUE)
          if (mx >= 1e4 || (!all(b == 0 | is.na(b)) && mx < 0.01)) {
            formatC(b, format = "e", digits = 2)
          } else {
            prettyNum(round(b, 3), big.mark = ",")
          }
        }
      )
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(name = y_lab_main %||% "")
  }

  # ---- 14. Labels ----------------------------------------------------------
  p <- p + ggplot2::labs(
    x        = x_lab %||% x_axis,
    title    = title,
    subtitle = subtitle,
    caption  = caption
  )

  # ---- 15. Faceting --------------------------------------------------------
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym(facet)),
      ncol = facet_ncol
    )
  }

  # ---- 16. Theme -----------------------------------------------------------
  p <- p + .ponder_theme_multi_y(
    theme_name       = theme,
    base_size        = global_font_size,
    fs_title         = title_font_size,
    fs_subtitle      = subtitle_font_size,
    fs_caption       = caption_font_size,
    fs_axis          = axis_font_size,
    fs_legend        = legend_font_size,
    fs_facet         = facet_font_size,
    legend_pos       = legend_pos,
    zoom             = zoom
  )

  p
}