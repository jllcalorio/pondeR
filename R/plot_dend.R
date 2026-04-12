# =============================================================================
# plot_dend.R
# Part of the pondeR package
# =============================================================================

#' @title Plot a Hierarchical Clustering Dendrogram
#'
#' @description
#' Produces a publication-ready dendrogram from a data frame, tibble, or
#' matrix using hierarchical clustering. Supports row-wise and column-wise
#' clustering with configurable distance metrics and agglomeration methods.
#' Returns a \pkg{ggplot2} object via \pkg{ggdendro}.
#'
#' @param x A data frame, tibble, or matrix of numeric values to cluster.
#'   Column names may contain special characters.
#' @param clustering_distance_rows A single character string specifying the
#'   distance metric for row clustering. One of \code{"euclidean"} (default),
#'   \code{"maximum"}, \code{"manhattan"}, \code{"canberra"},
#'   \code{"binary"}, or \code{"minkowski"}. See \code{\link[stats]{dist}}.
#' @param clustering_distance_cols A single character string specifying the
#'   distance metric for column clustering. Same options as
#'   \code{clustering_distance_rows}. Default \code{"euclidean"}.
#' @param agglomeration_method A single character string specifying the
#'   agglomeration method passed to \code{\link[stats]{hclust}}. One of
#'   \code{"ward.D"} (default), \code{"ward.D2"}, \code{"single"},
#'   \code{"complete"}, \code{"average"}, \code{"mcquitty"},
#'   \code{"median"}, or \code{"centroid"}.
#' @param orientation A single character string. One of \code{"rows"} (default;
#'   clusters rows), \code{"cols"} (clusters columns), or \code{"both"}
#'   (returns a named list of two \code{ggplot2} objects).
#' @param theme A single character string specifying the plot theme. One of
#'   \code{"nature"} (default), \code{"classic"}, \code{"minimal"}, or
#'   \code{"bw"}.
#' @param title A single character string for the plot title. Default
#'   \code{NULL} (no title).
#' @param subtitle A single character string for the plot subtitle. Default
#'   \code{NULL} (no subtitle).
#' @param global_font_size A single positive numeric specifying the base font
#'   size in points. Default \code{11}.
#' @param axis_title_size A single positive numeric for axis title font size.
#'   Defaults to \code{global_font_size}.
#' @param axis_text_size A single positive numeric for axis text font size.
#'   Defaults to \code{global_font_size}.
#' @param legend_title_size A single positive numeric for legend title font
#'   size. Defaults to \code{global_font_size}.
#' @param legend_text_size A single positive numeric for legend text font size.
#'   Defaults to \code{global_font_size}.
#' @param plot_title_size A single positive numeric for plot title font size.
#'   Defaults to \code{global_font_size + 2}.
#' @param plot_subtitle_size A single positive numeric for plot subtitle font
#'   size. Defaults to \code{global_font_size}.
#' @param label_size A single positive numeric controlling the font size of
#'   leaf labels. Defaults to \code{global_font_size * 0.8}.
#' @param line_width A single positive numeric for dendrogram branch width.
#'   Default \code{0.5}.
#' @param line_color A single character string for branch colour.
#'   Default \code{"black"}.
#' @param ... Additional arguments passed to \code{\link[stats]{dist}} or
#'   \code{\link[stats]{hclust}}.
#'
#' @details
#' Numeric columns are extracted and scaled (zero mean, unit variance) before
#' computing distances. Rows or columns containing all \code{NA} values are
#' silently dropped. Remaining \code{NA}s are imputed with column means before
#' scaling.
#'
#' When \code{orientation = "both"}, the returned object is a named list with
#' elements \code{rows} and \code{cols}.
#'
#' @return A \code{ggplot2} object, or a named list of two \code{ggplot2}
#'   objects (\code{rows} and \code{cols}) when \code{orientation = "both"}.
#'
#' @examples
#' ## Row dendrogram (default)
#' plot_dend(mtcars)
#'
#' ## Column dendrogram with Manhattan distance
#' plot_dend(mtcars, orientation = "cols",
#'           clustering_distance_cols = "manhattan")
#'
#' ## Both axes
#' dends <- plot_dend(mtcars, orientation = "both")
#' dends$rows
#' dends$cols
#'
#' ## Custom styling
#' plot_dend(mtcars, theme = "minimal", line_color = "#2166AC",
#'           line_width = 0.8, global_font_size = 10,
#'           title = "Sample Clustering")
#'
#' @author John Lennon L. Calorio
#' @export
plot_dend <- function(x,
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",
                       agglomeration_method     = "ward.D",
                       orientation              = "rows",
                       theme                    = "nature",
                       title                    = NULL,
                       subtitle                 = NULL,
                       global_font_size         = 11,
                       axis_title_size          = NULL,
                       axis_text_size           = NULL,
                       legend_title_size        = NULL,
                       legend_text_size         = NULL,
                       plot_title_size          = NULL,
                       plot_subtitle_size       = NULL,
                       label_size               = NULL,
                       line_width               = 0.5,
                       line_color               = "black",
                       ...) {

  # --- Dependencies -----------------------------------------------------------
  if (!requireNamespace("ggplot2",  quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')",
         call. = FALSE)
  if (!requireNamespace("ggdendro", quietly = TRUE))
    stop("Package 'ggdendro' is required. Install with: install.packages('ggdendro')",
         call. = FALSE)

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x))
    stop("`x` must be a data frame, tibble, or matrix.", call. = FALSE)
  if (is.matrix(x)) x <- as.data.frame(x)

  valid_dist <- c("euclidean", "maximum", "manhattan",
                  "canberra", "binary", "minkowski")
  if (!clustering_distance_rows %in% valid_dist)
    stop("`clustering_distance_rows` must be one of: ",
         paste0('"', valid_dist, '"', collapse = ", "), ".", call. = FALSE)
  if (!clustering_distance_cols %in% valid_dist)
    stop("`clustering_distance_cols` must be one of: ",
         paste0('"', valid_dist, '"', collapse = ", "), ".", call. = FALSE)

  valid_agg <- c("ward.D", "ward.D2", "single", "complete",
                 "average", "mcquitty", "median", "centroid")
  if (!agglomeration_method %in% valid_agg)
    stop("`agglomeration_method` must be one of: ",
         paste0('"', valid_agg, '"', collapse = ", "), ".", call. = FALSE)

  valid_orient <- c("rows", "cols", "both")
  if (!orientation %in% valid_orient)
    stop("`orientation` must be one of: ",
         paste0('"', valid_orient, '"', collapse = ", "), ".", call. = FALSE)

  valid_themes <- c("nature", "classic", "minimal", "bw")
  if (!theme %in% valid_themes)
    stop("`theme` must be one of: ",
         paste0('"', valid_themes, '"', collapse = ", "), ".", call. = FALSE)

  if (!is.null(title) && (!is.character(title) || length(title) != 1L))
    stop("`title` must be a single character string or NULL.", call. = FALSE)
  if (!is.null(subtitle) && (!is.character(subtitle) || length(subtitle) != 1L))
    stop("`subtitle` must be a single character string or NULL.", call. = FALSE)

  if (!is.numeric(global_font_size) || length(global_font_size) != 1L ||
      global_font_size <= 0)
    stop("`global_font_size` must be a single positive numeric.", call. = FALSE)

  .fs <- function(param, default) {
    if (!is.null(param)) {
      if (!is.numeric(param) || length(param) != 1L || param <= 0)
        stop(deparse(substitute(param)),
             " must be a single positive numeric.", call. = FALSE)
      param
    } else {
      default
    }
  }

  fs_axis_title    <- .fs(axis_title_size,    global_font_size)
  fs_axis_text     <- .fs(axis_text_size,     global_font_size)
  fs_legend_title  <- .fs(legend_title_size,  global_font_size)
  fs_legend_text   <- .fs(legend_text_size,   global_font_size)
  fs_plot_title    <- .fs(plot_title_size,    global_font_size + 2)
  fs_plot_subtitle <- .fs(plot_subtitle_size, global_font_size)
  fs_label         <- .fs(label_size,         global_font_size * 0.8)

  if (!is.numeric(line_width) || length(line_width) != 1L || line_width <= 0)
    stop("`line_width` must be a single positive numeric.", call. = FALSE)
  if (!is.character(line_color) || length(line_color) != 1L)
    stop("`line_color` must be a single character string.", call. = FALSE)

  # --- Prepare numeric matrix -------------------------------------------------
  num_cols <- vapply(x, is.numeric, logical(1L))
  if (!any(num_cols))
    stop("No numeric columns found in `x`.", call. = FALSE)
  mat <- as.matrix(x[, num_cols, drop = FALSE])

  row_ok <- rowSums(!is.na(mat)) > 0L
  col_ok <- colSums(!is.na(mat)) > 0L
  mat    <- mat[row_ok, col_ok, drop = FALSE]

  if (nrow(mat) < 2L)
    stop("At least 2 non-missing rows are required for row clustering.",
         call. = FALSE)
  if (ncol(mat) < 2L)
    stop("At least 2 non-missing columns are required for column clustering.",
         call. = FALSE)

  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx))
      mat[na_idx, j] <- mean(mat[!na_idx, j], na.rm = TRUE)
  }
  mat_scaled <- scale(mat)

  # --- Helper: build one dendrogram ggplot ------------------------------------
  .build_dend_gg <- function(data_mat, dist_method, flip,
                              plot_title, plot_subtitle) {

    extra      <- list(...)
    dist_extra <- extra[names(extra) %in% names(formals(stats::dist))]
    hc_extra   <- extra[names(extra) %in% names(formals(stats::hclust))]

    d   <- do.call(stats::dist,
                   c(list(x = data_mat, method = dist_method), dist_extra))
    hc  <- do.call(stats::hclust,
                   c(list(d = d, method = agglomeration_method), hc_extra))
    ddr <- ggdendro::dendro_data(hc, type = "rectangle")
    seg <- ggdendro::segment(ddr)
    lab <- ggdendro::label(ddr)

    # Always build ggplot first, then add theme and theme() overrides
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data      = seg,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = line_width,
        colour    = line_color
      ) +
      ggplot2::geom_text(
        data    = lab,
        ggplot2::aes(x = x, y = y, label = label),
        hjust   = 1,
        size    = fs_label / ggplot2::.pt,
        nudge_y = -max(seg$y, na.rm = TRUE) * 0.02
      ) +
      ggplot2::labs(title = plot_title, subtitle = plot_subtitle)

    if (flip) p <- p + ggplot2::coord_flip(clip = "off")

    base_theme_fn <- switch(theme,
      nature  = ggplot2::theme_bw,
      classic = ggplot2::theme_classic,
      minimal = ggplot2::theme_minimal,
      bw      = ggplot2::theme_bw
    )

    p +
      base_theme_fn(base_size = global_font_size) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_blank(),
        axis.title       = ggplot2::element_blank(),
        axis.text        = ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(
          size = fs_plot_title, face = "bold"),
        plot.subtitle    = ggplot2::element_text(size = fs_plot_subtitle),
        legend.title     = ggplot2::element_text(size = fs_legend_title),
        legend.text      = ggplot2::element_text(size = fs_legend_text)
      )
  }

  # --- Build and return -------------------------------------------------------
  if (orientation == "rows")
    return(.build_dend_gg(mat_scaled, clustering_distance_rows,
                          flip = TRUE,  title, subtitle))

  if (orientation == "cols")
    return(.build_dend_gg(t(mat_scaled), clustering_distance_cols,
                          flip = FALSE, title, subtitle))

  list(
    rows = .build_dend_gg(mat_scaled,    clustering_distance_rows,
                          flip = TRUE,  title, subtitle),
    cols = .build_dend_gg(t(mat_scaled), clustering_distance_cols,
                          flip = FALSE, title, subtitle)
  )
}