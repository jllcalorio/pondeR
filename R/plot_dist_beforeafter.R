#' Plot Distribution Comparison Before and After Data Transformation
#'
#' @title Plot Distribution Comparison Before and After Data Transformation
#'
#' @description
#' Generates density and box plots comparing data distributions before and
#' after any transformation, normalization, scaling, or correction step.
#' Plots can be oriented by sample or by feature. Returns a named list of
#' individual \code{ggplot2} objects and a combined \code{gtable} grob.
#'
#' @param x A \code{data.frame}, \code{tibble}, or \code{matrix} representing
#'   the data \emph{before} transformation. Rows are samples; columns are
#'   features. Column names may contain special characters.
#' @param y A \code{data.frame}, \code{tibble}, or \code{matrix} representing
#'   the data \emph{after} transformation. Must have the same dimensions and
#'   column names as \code{x}.
#' @param group_by Character. Perspective for box plots. \code{"sample"}
#'   shows one box per sample; \code{"feature"} shows one box per feature.
#'   Default is \code{"sample"}.
#' @param plot_what Character vector or \code{NULL}. Subset of column names
#'   (features) from \code{x} and \code{y} to include in the plots. If
#'   \code{NULL}, all features are used. Default is \code{NULL}.
#' @param n_random Integer or \code{NULL}. Maximum number of samples (when
#'   \code{group_by = "sample"}) or features (when \code{group_by =
#'   "feature"}) to display in box plots. If \code{NULL} or if the data has
#'   fewer items than \code{n_random}, all are shown. Default is \code{30}.
#' @param x_label Character. Label describing the before state, used in plot
#'   titles. Default is \code{"Before"}.
#' @param y_label Character. Label describing the after state, used in plot
#'   titles. Default is \code{"After"}.
#' @param x_fill Character. Fill color for before density and box plots.
#'   Default is \code{"#56B4E9"} (Okabe-Ito sky blue).
#' @param y_fill Character. Fill color for after density and box plots.
#'   Default is \code{"#E69F00"} (Okabe-Ito orange).
#' @param point_alpha Numeric in \code{[0, 1]}. Transparency for density fill
#'   and box fill. Default is \code{0.6}.
#' @param theme Character. ggplot2 theme to apply. Options: \code{"nature"},
#'   \code{"minimal"}, \code{"classic"}, \code{"bw"}, \code{"light"},
#'   \code{"dark"}. Default is \code{"nature"}.
#' @param base_size Numeric. Base font size for the theme (pts). Default is
#'   \code{11}.
#' @param font_family Character. Font family for all text elements. Default
#'   is \code{"sans"}.
#' @param plot_cols Integer or \code{NULL}. Number of columns in the combined
#'   output grid. If \code{NULL}, defaults to \code{2} (density row +
#'   boxplot row). Default is \code{NULL}.
#' @param seed Numeric or \code{NULL}. Random seed for reproducible sampling
#'   when \code{n_random} is used. Default is \code{123}.
#' @param global_font_size Numeric or \code{NULL}. Base font size for
#'   proportional scaling of all text elements. Default is \code{NULL}.
#' @param title_size Numeric or \code{NULL}. Plot panel title font size.
#'   Default is \code{NULL}.
#' @param subtitle_size Numeric or \code{NULL}. Combined plot main title font
#'   size. Default is \code{NULL}.
#' @param xlab_size Numeric or \code{NULL}. X-axis label font size. Default
#'   is \code{NULL}.
#' @param ylab_size Numeric or \code{NULL}. Y-axis label font size. Default
#'   is \code{NULL}.
#'
#' @details
#' \code{x} and \code{y} must share identical row counts, column counts, and
#' column names. Row order is assumed to correspond between \code{x} and
#' \code{y}.
#'
#' When \code{plot_what} is specified, only those features are included in
#' both density and box plots. When \code{n_random} is also specified, random
#' sub-sampling is applied on top of the \code{plot_what} subset for box
#' plots only; density plots always use the full \code{plot_what} subset.
#'
#' Box plot items are ordered by median value for readability.
#'
#' Font sizes follow the same resolution priority as \code{plot_beforeafter}:
#' explicit argument > \code{global_font_size} scaling > \code{base_size}
#' derived default.
#'
#' The function renders the combined plot to the active graphics device and
#' also returns a named list so individual panels can be extracted,
#' re-styled, or exported independently.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{plot_density_before}}{Density plot of \code{x}.}
#'   \item{\code{plot_density_after}}{Density plot of \code{y}.}
#'   \item{\code{plot_box_before}}{Box plot of \code{x}.}
#'   \item{\code{plot_box_after}}{Box plot of \code{y}.}
#'   \item{\code{plot_combined}}{A \code{gtable} grob of all four panels,
#'     drawn to the active device.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n_samples  <- 40
#' n_features <- 20
#'
#' x <- as.data.frame(
#'   matrix(abs(rnorm(n_samples * n_features, mean = 1000, sd = 300)),
#'          nrow = n_samples,
#'          dimnames = list(
#'            paste0("S", seq_len(n_samples)),
#'            paste0("Feature_", seq_len(n_features))
#'          ))
#' )
#' y <- as.data.frame(scale(x))
#'
#' # Sample-wise comparison, default settings
#' result <- plot_dist_beforeafter(
#'   x       = x,
#'   y       = y,
#'   group_by = "sample"
#' )
#'
#' # Feature-wise, subset of features, classic theme
#' result <- plot_dist_beforeafter(
#'   x         = x,
#'   y         = y,
#'   group_by  = "feature",
#'   plot_what = paste0("Feature_", 1:10),
#'   theme     = "classic"
#' )
#'
#' # Extract individual panels
#' result$plot_density_before
#' result$plot_box_after
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_boxplot coord_flip labs
#'   theme_bw theme_minimal theme_classic theme_light theme_dark theme
#'   element_text element_blank element_line element_rect
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.newpage grid.draw textGrob gpar
#'
#' @export
plot_dist_beforeafter <- function(
    x,
    y,
    group_by         = "sample",
    plot_what        = NULL,
    n_random         = 30L,
    x_label          = "Before",
    y_label          = "After",
    x_fill           = "#56B4E9",
    y_fill           = "#E69F00",
    point_alpha      = 0.6,
    theme            = "nature",
    base_size        = 11,
    font_family      = "sans",
    plot_cols        = NULL,
    seed             = 123,
    global_font_size = NULL,
    title_size       = NULL,
    subtitle_size    = NULL,
    xlab_size        = NULL,
    ylab_size        = NULL
) {

  # -------------------------------------------------------------------------
  # Package availability
  # -------------------------------------------------------------------------
  for (pkg in c("ggplot2", "gridExtra", "grid")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Required package '", pkg, "' is not installed. ",
        "Install it with: install.packages('", pkg, "')"
      )
    }
  }

  # --- Integration with run_DIpreprocess ---
  if (inherits(x, "run_DIpreprocess")) {
    if (missing(metadata) || is.null(metadata)) {
      # Must use unmerged metadata because data_imputed is unmerged
      metadata <- x$metadata
    }
    
    # Extract states
    data_before <- x$data_imputed
    data_after  <- x$data_nonpls
    
    # Subset to common features to ensure exact alignment
    common_features <- intersect(colnames(data_before), colnames(data_after))
    if (length(common_features) == 0) {
      stop("No common features found between imputed and final non-PLS data.")
    }
    
    # Overwrite x and y
    x <- data_before[, common_features, drop = FALSE]
    y <- data_after[, common_features, drop = FALSE]
  }

  # -------------------------------------------------------------------------
  # Input validation: x and y
  # -------------------------------------------------------------------------
  .is_tabular <- function(z) is.data.frame(z) || is.matrix(z)
  if (!.is_tabular(x)) stop("'x' must be a data.frame, tibble, or matrix.")
  if (!.is_tabular(y)) stop("'y' must be a data.frame, tibble, or matrix.")

  if (!identical(dim(x), dim(y))) {
    stop(
      "'x' and 'y' must have the same dimensions. ",
      "Got: x = ", nrow(x), " x ", ncol(x), ", ",
      "y = ",      nrow(y), " x ", ncol(y), "."
    )
  }

  x_cols <- if (is.matrix(x)) colnames(x) else names(x)
  y_cols <- if (is.matrix(y)) colnames(y) else names(y)
  if (!identical(x_cols, y_cols)) {
    stop("'x' and 'y' must have identical column names in the same order.")
  }
  if (is.null(x_cols)) stop("'x' and 'y' must have named columns.")

  x <- as.data.frame(x, check.names = FALSE)
  y <- as.data.frame(y, check.names = FALSE)

  # -------------------------------------------------------------------------
  # Input validation: group_by
  # -------------------------------------------------------------------------
  group_by <- tolower(group_by)
  if (!group_by %in% c("sample", "feature")) {
    stop("'group_by' must be \"sample\" or \"feature\".")
  }

  # -------------------------------------------------------------------------
  # Input validation: plot_what
  # -------------------------------------------------------------------------
  if (!is.null(plot_what)) {
    if (!is.character(plot_what) || length(plot_what) < 1L) {
      stop("'plot_what' must be a non-empty character vector or NULL.")
    }
    missing_feats <- setdiff(plot_what, x_cols)
    if (length(missing_feats) > 0L) {
      stop(
        "The following feature(s) in 'plot_what' were not found in 'x'/'y': ",
        paste(missing_feats, collapse = ", "), "."
      )
    }
    x <- x[, plot_what, drop = FALSE]
    y <- y[, plot_what, drop = FALSE]
  }

  # -------------------------------------------------------------------------
  # Input validation: scalars and theme
  # -------------------------------------------------------------------------
  .check_positive_numeric <- function(val, name, allow_null = FALSE) {
    if (allow_null && is.null(val)) return(invisible(NULL))
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop("'", name, "' must be a single positive number.")
    }
  }

  .check_positive_numeric(point_alpha,      "point_alpha")
  .check_positive_numeric(base_size,        "base_size")
  .check_positive_numeric(global_font_size, "global_font_size", allow_null = TRUE)
  .check_positive_numeric(title_size,       "title_size",       allow_null = TRUE)
  .check_positive_numeric(subtitle_size,    "subtitle_size",    allow_null = TRUE)
  .check_positive_numeric(xlab_size,        "xlab_size",        allow_null = TRUE)
  .check_positive_numeric(ylab_size,        "ylab_size",        allow_null = TRUE)

  if (point_alpha > 1) {
    warning("'point_alpha' is > 1; it will be clamped to 1.")
    point_alpha <- 1
  }

  if (!is.null(n_random)) {
    if (!is.numeric(n_random) || length(n_random) != 1L ||
        is.na(n_random) || n_random < 1L || n_random != floor(n_random)) {
      stop("'n_random' must be a single positive integer or NULL.")
    }
    n_random <- as.integer(n_random)
  }

  if (!is.null(plot_cols)) {
    if (!is.numeric(plot_cols) || length(plot_cols) != 1L ||
        is.na(plot_cols) || plot_cols < 1L || plot_cols != floor(plot_cols)) {
      stop("'plot_cols' must be a single positive integer or NULL.")
    }
    plot_cols <- as.integer(plot_cols)
  }

  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  theme_choice <- tolower(theme)
  if (!theme_choice %in% valid_themes) {
    stop(
      "'theme' must be one of: ", paste(valid_themes, collapse = ", "), ". ",
      "Got: '", theme, "'."
    )
  }

  if (!is.character(font_family) || length(font_family) != 1L) {
    stop("'font_family' must be a single character string.")
  }
  if (!is.character(x_fill) || length(x_fill) != 1L) {
    stop("'x_fill' must be a single character color string.")
  }
  if (!is.character(y_fill) || length(y_fill) != 1L) {
    stop("'y_fill' must be a single character color string.")
  }

  if (identical(x, y)) {
    warning("'x' and 'y' are identical. No transformation may have been applied.")
  }

  # -------------------------------------------------------------------------
  # Seed
  # -------------------------------------------------------------------------
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("'seed' must be a single numeric value or NULL.")
    }
    set.seed(seed)
  }

  # -------------------------------------------------------------------------
  # Font size resolution
  # -------------------------------------------------------------------------
  .resolve_size <- function(explicit, scale_factor, base_default) {
    if (!is.null(explicit))     return(explicit)
    if (!is.null(scale_factor)) return(base_default * scale_factor)
    base_default
  }

  gfs_factor      <- if (!is.null(global_font_size)) global_font_size / 11 else NULL
  r_title_size    <- .resolve_size(title_size,    gfs_factor, base_size + 1)
  r_subtitle_size <- .resolve_size(subtitle_size, gfs_factor, base_size + 2)
  r_xlab_size     <- .resolve_size(xlab_size,     gfs_factor, base_size)
  r_ylab_size     <- .resolve_size(ylab_size,     gfs_factor, base_size)
  r_axis_text     <- .resolve_size(NULL,          gfs_factor, base_size - 1)

  # -------------------------------------------------------------------------
  # Theme construction
  # -------------------------------------------------------------------------
  base_theme <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.6)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size     = base_size, base_family = font_family),
    "light"   = ggplot2::theme_light(base_size  = base_size, base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size   = base_size, base_family = font_family)
  )

  custom_theme <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title       = ggplot2::element_text(size = r_title_size,  hjust = 0.5,
                                             face = "bold", family = font_family),
    axis.title.x     = ggplot2::element_text(size = r_xlab_size,   family = font_family),
    axis.title.y     = ggplot2::element_text(size = r_ylab_size,   family = font_family),
    axis.text        = ggplot2::element_text(size = r_axis_text,   family = font_family)
  )

  # -------------------------------------------------------------------------
  # Convert to long format (base R, no tidyr dependency)
  # -------------------------------------------------------------------------
  .to_long <- function(df, id_col_name, value_name) {
    n_rows <- nrow(df)
    n_cols <- ncol(df)
    ids    <- if (!is.null(rownames(df))) rownames(df) else as.character(seq_len(n_rows))
    data.frame(
      id    = rep(ids, times = n_cols),
      key   = rep(names(df), each  = n_rows),
      value = unlist(df, use.names = FALSE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  long_x <- .to_long(x)
  long_y <- .to_long(y)

  # -------------------------------------------------------------------------
  # Build density plots (all features, all samples)
  # -------------------------------------------------------------------------
  p_density_before <- ggplot2::ggplot(long_x, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = x_fill, color = x_fill, alpha = point_alpha) +
    ggplot2::labs(
      title = paste("Distribution —", x_label),
      x     = "Value",
      y     = "Density"
    ) +
    base_theme + custom_theme

  p_density_after <- ggplot2::ggplot(long_y, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = y_fill, color = y_fill, alpha = point_alpha) +
    ggplot2::labs(
      title = paste("Distribution —", y_label),
      x     = "Value",
      y     = "Density"
    ) +
    base_theme + custom_theme

  # -------------------------------------------------------------------------
  # Subsample for box plots
  # -------------------------------------------------------------------------
  # Items = samples (group_by = "sample") or features (group_by = "feature")
  if (group_by == "sample") {
    all_items    <- unique(long_x$id)
    item_col     <- "id"
    box_x_label  <- "Sample"
  } else {
    all_items    <- unique(long_x$key)
    item_col     <- "key"
    box_x_label  <- "Feature"
  }

  n_items <- length(all_items)
  if (!is.null(n_random) && n_random < n_items) {
    selected_items <- sample(all_items, n_random)
    n_shown        <- n_random
  } else {
    selected_items <- all_items
    n_shown        <- n_items
  }

  sub_x <- long_x[long_x[[item_col]] %in% selected_items, , drop = FALSE]
  sub_y <- long_y[long_y[[item_col]] %in% selected_items, , drop = FALSE]

  # Order items by median for readability
  med_order_x <- names(sort(tapply(sub_x$value, sub_x[[item_col]], median, na.rm = TRUE)))
  med_order_y <- names(sort(tapply(sub_y$value, sub_y[[item_col]], median, na.rm = TRUE)))

  sub_x[[item_col]] <- factor(sub_x[[item_col]], levels = med_order_x)
  sub_y[[item_col]] <- factor(sub_y[[item_col]], levels = med_order_y)

  box_title_before <- paste0(box_x_label, " Profiles — ", x_label,
                             " (n = ", n_shown, ")")
  box_title_after  <- paste0(box_x_label, " Profiles — ", y_label,
                             " (n = ", n_shown, ")")

  p_box_before <- ggplot2::ggplot(
    sub_x,
    ggplot2::aes(x = .data[[item_col]], y = value)
  ) +
    ggplot2::geom_boxplot(fill = x_fill, color = "grey30",
                          alpha = point_alpha, outlier.size = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = box_title_before, x = box_x_label, y = "Value") +
    base_theme + custom_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(6, r_axis_text - 2),
                                                       family = font_family))

  p_box_after <- ggplot2::ggplot(
    sub_y,
    ggplot2::aes(x = .data[[item_col]], y = value)
  ) +
    ggplot2::geom_boxplot(fill = y_fill, color = "grey30",
                          alpha = point_alpha, outlier.size = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = box_title_after, x = box_x_label, y = "Value") +
    base_theme + custom_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(6, r_axis_text - 2),
                                                       family = font_family))

  # -------------------------------------------------------------------------
  # Assemble combined plot
  # -------------------------------------------------------------------------
  main_title <- paste0(
    "Distribution Comparison: ", x_label, " vs ", y_label,
    " (", tools::toTitleCase(group_by), " Perspective)"
  )

  n_cols_out <- if (!is.null(plot_cols)) plot_cols else 2L

  combined_plot <- gridExtra::arrangeGrob(
    p_density_before, p_density_after,
    p_box_before,     p_box_after,
    ncol = n_cols_out,
    top  = grid::textGrob(
      main_title,
      gp = grid::gpar(
        fontsize   = r_subtitle_size,
        fontface   = "bold",
        fontfamily = font_family
      )
    )
  )

  grid::grid.newpage()
  grid::grid.draw(combined_plot)

  message("Distribution comparison plots created successfully.")

  return(invisible(list(
    plot_density_before = p_density_before,
    plot_density_after  = p_density_after,
    plot_box_before     = p_box_before,
    plot_box_after      = p_box_after,
    plot_combined       = combined_plot
  )))
}