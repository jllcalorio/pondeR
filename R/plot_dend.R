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
#'   \code{"nature"} (default), \code{"classic"}, \code{"plos"},
#'   \code{"minimal"}, \code{"bw"}, \code{"light"}, \code{"dark"}, or
#'   \code{"gray"}.
#' @param title A single character string for the plot title. Default
#'   \code{NULL} (no title).
#' @param subtitle A single character string for the plot subtitle. Default
#'   \code{NULL} (no subtitle).
#' @param base_size Base font size in points; every text element scales from
#'   it. Default \code{11}, which suits a default 7 x 5 inch Quarto figure.
#' @param title_size,subtitle_size,label_size Font sizes in points for the
#'   title, subtitle and leaf labels. Default \code{NULL} derives each from
#'   \code{base_size}.
#' @param leaf_label_angle A single numeric specifying the rotation angle (in
#'   degrees) for leaf label text. Default \code{NULL}. When \code{NULL}, the
#'   angle defaults to 0 if \code{orientation = "rows"} (horizontal labels) or
#'   90 if \code{orientation = "cols"} (vertical labels).
#' @param leaf_labels A character vector of custom leaf names. Default
#'   \code{NULL}. Must have one name per row for \code{orientation = "rows"} and
#'   one per column for \code{orientation = "cols"}; not used with
#'   \code{orientation = "both"}.
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
#' @references
#' Ward JH (1963). Hierarchical grouping to optimize an objective function.
#'   \emph{Journal of the American Statistical Association}, 58(301), 236-244.
#'   \doi{10.1080/01621459.1963.10500845}
#'
#' Murtagh F, Contreras P (2012). Algorithms for hierarchical clustering: an
#'   overview. \emph{WIREs Data Mining and Knowledge Discovery}, 2(1), 86-97.
#'   \doi{10.1002/widm.53}
#'
#' @examplesIf requireNamespace("ggdendro", quietly = TRUE)
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
#'           line_width = 0.8, base_size = 10,
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
                       base_size                = 11,
                       title_size               = NULL,
                       subtitle_size            = NULL,
                       label_size               = NULL,
                       leaf_label_angle         = NULL,
                       leaf_labels              = NULL,
                       line_width               = 0.5,
                       line_color               = "black",
                       ...) {

  # --- Dependencies -----------------------------------------------------------
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

  if (!is.null(title) && (!is.character(title) || length(title) != 1L))
    stop("`title` must be a single character string or NULL.", call. = FALSE)
  if (!is.null(subtitle) && (!is.character(subtitle) || length(subtitle) != 1L))
    stop("`subtitle` must be a single character string or NULL.", call. = FALSE)

  sizes <- .text_sizes(base_size, title_size, subtitle_size, label = label_size)
  if (!is.numeric(line_width) || length(line_width) != 1L || line_width <= 0)
    stop("`line_width` must be a single positive numeric.", call. = FALSE)
  if (!is.character(line_color) || length(line_color) != 1L)
    stop("`line_color` must be a single character string.", call. = FALSE)

  # Determine leaf label angle based on orientation if not provided
  if (!is.null(leaf_label_angle) &&
      (!is.numeric(leaf_label_angle) || length(leaf_label_angle) != 1L))
    stop("`leaf_label_angle` must be a single numeric value or NULL.",
         call. = FALSE)

  # Validate leaf_labels if provided
  if (!is.null(leaf_labels)) {
    if (!is.character(leaf_labels))
      stop("`leaf_labels` must be a character vector or NULL.", call. = FALSE)
    if (orientation == "both")
      stop("`leaf_labels` needs `orientation = \"rows\"` or `\"cols\"`.", call. = FALSE)
  }

  # --- Prepare numeric matrix -------------------------------------------------
  num_cols <- vapply(x, is.numeric, logical(1L))
  if (!any(num_cols))
    stop("No numeric columns found in `x`.", call. = FALSE)
  mat <- as.matrix(x[, num_cols, drop = FALSE])

  # Apply custom row names if provided (before any row filtering)
  if (!is.null(leaf_labels)) {
    n_leaf <- if (orientation == "rows") nrow(mat) else ncol(mat)
    if (length(leaf_labels) != n_leaf)
      stop("`leaf_labels` must have ", n_leaf, " element(s), one per ",
           if (orientation == "rows") "row" else "numeric column", " of `x`.",
           call. = FALSE)
    if (orientation == "rows") rownames(mat) <- leaf_labels
    else colnames(mat) <- leaf_labels
  }

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

    # Leaf names are axis labels, so ggplot2 reserves the room they need
    angle     <- leaf_label_angle %||% (if (flip) 0 else 90)
    leaf_text <- ggplot2::element_text(
      size  = sizes$label,
      angle = angle,
      hjust = if (flip || angle != 0) 1 else 0.5,
      vjust = if (flip || angle != 0) 0.5 else 1
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data      = seg,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = line_width,
        colour    = line_color
      ) +
      ggplot2::scale_x_continuous(breaks = lab$x, labels = as.character(lab$label),
                                  expand = ggplot2::expansion(add = 0.5)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(title = plot_title, subtitle = plot_subtitle, x = NULL, y = NULL)

    if (flip) p <- p + ggplot2::coord_flip()

    p +
      .ponder_theme(theme, sizes) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_blank(),
        axis.text.x      = if (flip) ggplot2::element_blank() else leaf_text,
        axis.text.y      = if (flip) leaf_text else ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_blank()
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
