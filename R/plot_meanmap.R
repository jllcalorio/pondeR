# =============================================================================
# plot_meanmap.R
# Part of the pondeR package
# =============================================================================

#' @title Mean Intensity Heatmap with Hierarchical Clustering Dendrogram
#'
#' @description
#' Produces a publication-ready mean intensity heatmap where rows are
#' variables (e.g. metabolites) clustered by hierarchical clustering and
#' columns are group means derived from a grouping variable in the metadata.
#' A coloured dendrogram is aligned to the left of the heatmap, with optional
#' cluster bracket annotations and per-row significance markers on the right.
#'
#' This function is intended for visualising group-level mean differences
#' across variables, for example, metabolite intensities across disease
#' stages.
#'
#' @param x A numeric matrix or data frame whose \strong{columns} are the
#'   variables (e.g. metabolites) and whose \strong{rows} are
#'   observations/samples. Must have column names.
#' @param metadata A data frame with the same number of rows as \code{x} and
#'   valid column names. Must contain the column specified in \code{group}.
#' @param group A single character string naming the column in \code{metadata}
#'   that defines sample groups. Group means are computed per level of this
#'   variable and used as heatmap columns.
#' @param group_order A character vector specifying the display order of group
#'   levels (left to right). Must contain exactly the same levels as
#'   \code{metadata[[group]]}. When \code{NULL} (default), levels are ordered
#'   as they appear in \code{factor(metadata[[group]])}.
#' @param sig_vars A character vector of variable (column) names to mark with
#'   a significance marker on the right side of the heatmap. Default
#'   \code{NULL} (no markers).
#' @param sig_marker A single character string used as the significance marker
#'   symbol. Default \code{"*"}.
#' @param sig_color A single character string (colour name or hex) for the
#'   significance marker. Default \code{"red"}.
#' @param clust_method A single character string passed to the
#'   \code{method} argument of \code{\link[stats]{hclust}}. Default
#'   \code{"ward.D2"}.
#' @param dist_method A single character string passed to the \code{method}
#'   argument of \code{\link[stats]{dist}}. Default \code{"euclidean"}.
#' @param scale_vars A single character string controlling row-wise scaling
#'   applied before clustering and display. One of \code{"none"} (default;
#'   use values as supplied), \code{"zscore"} (subtract row mean, divide by
#'   row SD), or \code{"minmax"} (scale each row to (0, 1) inclusive).
#' @param n_clusters A positive integer or \code{NULL}. When not \code{NULL},
#'   the dendrogram is cut into \code{n_clusters} clusters using
#'   \code{\link[stats]{cutree}}, each cluster is coloured with a distinct
#'   colour-blind-friendly palette, and cluster bracket + label annotations
#'   are added to the left of the dendrogram. Default \code{NULL}.
#' @param cluster_labels A character vector of labels for each cluster, in
#'   order from the top of the plot to the bottom (i.e. in dendrogram leaf
#'   order). Must have length equal to \code{n_clusters} when supplied.
#'   Default \code{NULL}, which auto-generates labels
#'   \code{"Cluster I"}, \code{"Cluster II"}, etc.
#' @param cluster_colors A character vector of colours for each cluster.
#'   Must have length equal to \code{n_clusters} when supplied. Default
#'   \code{NULL}, which uses the Okabe-Ito colour-blind-safe palette.
#' @param heatmap_colors A character vector of length \eqn{\ge 2} specifying
#'   the colour gradient for the heatmap fill, from low to high values.
#'   Default \code{c("#008000", "#FFFF00", "#FF0000")} (green–yellow–red).
#' @param show_values Logical. Whether to print the rounded mean value inside
#'   each heatmap cell. Default \code{FALSE}.
#' @param value_digits A single non-negative integer. Number of decimal places
#'   for in-cell value labels. Default \code{2}. Ignored when
#'   \code{show_values = FALSE}.
#' @param value_color A single character string for in-cell text colour.
#'   Default \code{"black"}.
#' @param global_font_size Numeric. Base font size in points. Default
#'   \code{11}.
#' @param row_font_size Numeric or \code{NULL}. Font size for row (variable)
#'   labels on the right y-axis. Derived from \code{global_font_size} when
#'   \code{NULL} (default).
#' @param col_font_size Numeric or \code{NULL}. Font size for column (group)
#'   labels on the x-axis. Derived from \code{global_font_size} when
#'   \code{NULL} (default).
#' @param legend_font_size Numeric or \code{NULL}. Font size for legend text.
#'   Derived from \code{global_font_size} when \code{NULL} (default).
#' @param cluster_label_size Numeric or \code{NULL}. Font size for cluster
#'   bracket labels. Derived from \code{global_font_size} when \code{NULL}
#'   (default). Ignored when \code{n_clusters} is \code{NULL}.
#' @param sig_marker_size Numeric or \code{NULL}. Font size for significance
#'   markers. Derived from \code{global_font_size} when \code{NULL} (default).
#' @param dend_width Numeric in (0, 1). Relative width of the dendrogram panel
#'   as a fraction of the total plot width. Default \code{0.2}.
#' @param legend_position A single character string for legend placement.
#'   One of \code{"right"} (default), \code{"bottom"}, \code{"left"},
#'   \code{"top"}, or \code{"none"}.
#' @param na_color A single character string for the fill colour of \code{NA}
#'   cells. Default \code{"grey80"}.
#'
#' @details
#' \strong{Algorithm overview}
#'
#' \enumerate{
#'   \item \strong{Group mean computation}: For each level of
#'         \code{metadata[[group]]}, the column means of \code{x} are computed
#'         (ignoring \code{NA}s), producing a \eqn{p \times g} matrix where
#'         \eqn{p} is the number of variables and \eqn{g} the number of groups.
#'   \item \strong{Optional scaling}: When \code{scale_vars} is
#'         \code{"zscore"} or \code{"minmax"}, each \strong{row} (variable) of
#'         the mean matrix is scaled independently before clustering and
#'         display.
#'   \item \strong{Clustering}: Euclidean (or user-specified) distances are
#'         computed on the rows of the (possibly scaled) mean matrix and
#'         hierarchical clustering is applied.
#'   \item \strong{Reordering}: Rows are permuted to the leaf order of the
#'         dendrogram.
#'   \item \strong{Dendrogram}: Built via \pkg{ggdendro} with optional cluster
#'         colouring. A bracket-and-label annotation layer is added when
#'         \code{n_clusters} is supplied.
#'   \item \strong{Heatmap}: Built with \code{ggplot2::geom_tile} using a
#'         user-controlled colour gradient. Significance markers and in-cell
#'         values are added as optional annotation layers.
#'   \item \strong{Assembly}: \code{patchwork::wrap_plots()} aligns the two
#'         panels with zero inter-panel spacing.
#' }
#'
#' \strong{Scaling recommendation}: If variables are on comparable scales
#' (e.g. all z-scored or log-transformed), \code{scale_vars = "none"} is
#' appropriate. For raw intensities on very different scales, use
#' \code{scale_vars = "zscore"} to make the colour gradient interpretable
#' across all rows.
#'
#' @return A \code{patchwork} object. The following attributes are attached:
#'   \describe{
#'     \item{\code{mean_matrix}}{The \eqn{p \times g} matrix of group means
#'       (unscaled), rows in dendrogram leaf order.}
#'     \item{\code{scaled_matrix}}{The scaled version of \code{mean_matrix}
#'       used for clustering and display (identical to \code{mean_matrix}
#'       when \code{scale_vars = "none"}).}
#'     \item{\code{hclust_obj}}{The \code{hclust} object.}
#'     \item{\code{leaf_order}}{Integer vector of row indices in dendrogram
#'       leaf order.}
#'     \item{\code{cluster_membership}}{Named integer vector of cluster
#'       assignments per variable (only set when \code{n_clusters} is not
#'       \code{NULL}).}
#'   }
#'
#' @examples
#' \donttest{
#' ## ── Basic example: iris mean intensity by species ────────────────────────
#' data(iris)
#' num_part <- iris[, 1:4]
#' meta_part <- data.frame(Species = iris$Species)
#'
#' p <- plot_meanmap(
#'   x         = num_part,
#'   metadata  = meta_part,
#'   group     = "Species",
#'   scale_vars = "zscore"
#' )
#' print(p)
#'
#' ## ── With cluster colouring and significance markers ──────────────────────
#' p2 <- plot_meanmap(
#'   x                = num_part,
#'   metadata         = meta_part,
#'   group            = "Species",
#'   scale_vars       = "zscore",
#'   n_clusters       = 2,
#'   cluster_labels   = c("Cluster I", "Cluster II"),
#'   sig_vars         = c("Petal.Length", "Petal.Width"),
#'   sig_color        = "red",
#'   global_font_size = 11
#' )
#' print(p2)
#'
#' ## ── Custom group order and colour palette ────────────────────────────────
#' p3 <- plot_meanmap(
#'   x           = num_part,
#'   metadata    = meta_part,
#'   group       = "Species",
#'   group_order = c("virginica", "versicolor", "setosa"),
#'   scale_vars  = "zscore",
#'   n_clusters  = 2,
#'   heatmap_colors = c("#2166AC", "#F7F7F7", "#B2182B")
#' )
#' print(p3)
#' }
#'
#' @author John Lennon L. Calorio
#' @export
plot_meanmap <- function(
    x,
    metadata          = NULL,
    group,
    group_order       = NULL,
    sig_vars          = NULL,
    sig_marker        = "*",
    sig_color         = "red",
    clust_method      = "ward.D2",
    dist_method       = "euclidean",
    scale_vars        = "none",
    n_clusters        = NULL,
    cluster_labels    = NULL,
    cluster_colors    = NULL,
    heatmap_colors    = c("#008000", "#FFFF00", "#FF0000"),
    show_values       = FALSE,
    value_digits      = 2L,
    value_color       = "black",
    global_font_size  = 11,
    row_font_size     = NULL,
    col_font_size     = NULL,
    legend_font_size  = NULL,
    cluster_label_size = NULL,
    sig_marker_size   = NULL,
    dend_width        = 0.2,
    legend_position   = "right",
    na_color          = "grey80"
) {

  # ── 0. Dependencies ──────────────────────────────────────────────────────────
  .need <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf(
        "Package '%s' is required but not installed. Install with: install.packages('%s')",
        pkg, pkg
      ), call. = FALSE)
  }
  .need("ggplot2")
  .need("patchwork")
  .need("ggdendro")

  # ── 1. Input validation ──────────────────────────────────────────────────────

  ## 1a. x
  if (!is.matrix(x) && !is.data.frame(x))
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`x` must contain only numeric values.", call. = FALSE)
  if (is.null(colnames(x)))
    stop("`x` must have column names.", call. = FALSE)
  if (ncol(x) < 1L)
    stop("`x` must have at least one column.", call. = FALSE)
  if (nrow(x) < 1L)
    stop("`x` must have at least one row.", call. = FALSE)

  ## 1b. metadata
  if (is.null(metadata) || !is.data.frame(metadata))
    stop("`metadata` must be a data frame.", call. = FALSE)
  if (nrow(metadata) != nrow(x))
    stop(sprintf(
      "`metadata` must have the same number of rows as `x` (%d), but has %d.",
      nrow(x), nrow(metadata)
    ), call. = FALSE)

  ## 1c. group
  if (!is.character(group) || length(group) != 1L)
    stop("`group` must be a single character string.", call. = FALSE)
  if (!group %in% colnames(metadata))
    stop(sprintf(
      "`group` column '%s' not found in `metadata`. Available columns: %s.",
      group, paste(colnames(metadata), collapse = ", ")
    ), call. = FALSE)
  grp_vec    <- as.character(metadata[[group]])
  grp_levels <- if (!is.null(group_order)) {
    if (!is.character(group_order))
      stop("`group_order` must be a character vector.", call. = FALSE)
    missing_lvl <- setdiff(unique(grp_vec), group_order)
    extra_lvl   <- setdiff(group_order, unique(grp_vec))
    if (length(missing_lvl) > 0L)
      stop(sprintf(
        "`group_order` is missing level(s): %s.",
        paste(missing_lvl, collapse = ", ")
      ), call. = FALSE)
    if (length(extra_lvl) > 0L)
      stop(sprintf(
        "`group_order` contains level(s) not present in `metadata[[\"%s\"]]`: %s.",
        group, paste(extra_lvl, collapse = ", ")
      ), call. = FALSE)
    group_order
  } else {
    levels(factor(grp_vec))
  }
  if (length(grp_levels) < 2L)
    stop("`group` must have at least 2 distinct levels.", call. = FALSE)

  ## 1d. scale_vars
  scale_vars <- match.arg(scale_vars, c("none", "zscore", "minmax"))

  ## 1e. clustering
  valid_clust <- c("ward.D", "ward.D2", "single", "complete", "average",
                   "mcquitty", "median", "centroid")
  clust_method <- match.arg(clust_method, valid_clust)
  valid_dist   <- c("euclidean", "maximum", "manhattan", "canberra",
                    "binary", "minkowski")
  dist_method  <- match.arg(dist_method, valid_dist)

  ## 1f. n_clusters
  if (!is.null(n_clusters)) {
    if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
        is.na(n_clusters) || n_clusters < 1L)
      stop("`n_clusters` must be a single positive integer.", call. = FALSE)
    n_clusters <- as.integer(n_clusters)
  }

  ## 1g. cluster_labels / cluster_colors
  if (!is.null(cluster_labels)) {
    if (is.null(n_clusters))
      stop("`cluster_labels` requires `n_clusters` to be set.", call. = FALSE)
    if (!is.character(cluster_labels) || length(cluster_labels) != n_clusters)
      stop(sprintf(
        "`cluster_labels` must be a character vector of length %d (= `n_clusters`).",
        n_clusters
      ), call. = FALSE)
  }
  if (!is.null(cluster_colors)) {
    if (is.null(n_clusters))
      stop("`cluster_colors` requires `n_clusters` to be set.", call. = FALSE)
    if (!is.character(cluster_colors) || length(cluster_colors) != n_clusters)
      stop(sprintf(
        "`cluster_colors` must be a character vector of length %d (= `n_clusters`).",
        n_clusters
      ), call. = FALSE)
    tryCatch(grDevices::col2rgb(cluster_colors),
             error = function(e) stop(
               "Invalid colour(s) in `cluster_colors`.", call. = FALSE))
  }

  ## 1h. heatmap_colors
  if (!is.character(heatmap_colors) || length(heatmap_colors) < 2L)
    stop("`heatmap_colors` must be a character vector of length >= 2.", call. = FALSE)
  tryCatch(grDevices::col2rgb(heatmap_colors),
           error = function(e) stop(
             "Invalid colour(s) in `heatmap_colors`.", call. = FALSE))

  ## 1i. sig_vars
  if (!is.null(sig_vars)) {
    if (!is.character(sig_vars))
      stop("`sig_vars` must be a character vector of column names.", call. = FALSE)
    unknown <- setdiff(sig_vars, colnames(x))
    if (length(unknown) > 0L)
      warning(sprintf(
        "%d name(s) in `sig_vars` not found in `colnames(x)` and will be ignored: %s.",
        length(unknown), paste(unknown, collapse = ", ")
      ), call. = FALSE)
    sig_vars <- intersect(sig_vars, colnames(x))
  }

  ## 1j. numerics
  .pos_num <- function(val, nm) {
    if (!is.null(val) &&
        (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0))
      stop(sprintf("`%s` must be a single positive numeric value or NULL.", nm),
           call. = FALSE)
  }
  if (!is.numeric(global_font_size) || length(global_font_size) != 1L ||
      global_font_size <= 0)
    stop("`global_font_size` must be a single positive numeric.", call. = FALSE)
  .pos_num(row_font_size,     "row_font_size")
  .pos_num(col_font_size,     "col_font_size")
  .pos_num(legend_font_size,  "legend_font_size")
  .pos_num(cluster_label_size,"cluster_label_size")
  .pos_num(sig_marker_size,   "sig_marker_size")
  if (!is.numeric(dend_width) || length(dend_width) != 1L ||
      dend_width <= 0 || dend_width >= 1)
    stop("`dend_width` must be a single numeric in (0, 1).", call. = FALSE)

  ## 1k. legend_position
  legend_position <- match.arg(
    legend_position, c("right", "bottom", "left", "top", "none")
  )

  ## 1l. logical flags
  if (!is.logical(show_values) || length(show_values) != 1L || is.na(show_values))
    stop("`show_values` must be a single non-NA logical.", call. = FALSE)

  # ── 2. Derived font sizes ────────────────────────────────────────────────────
  `%||%` <- function(a, b) if (is.null(a)) b else a
  fs_row   <- row_font_size      %||% (global_font_size * 0.85)
  fs_col   <- col_font_size      %||% (global_font_size * 0.85)
  fs_leg   <- legend_font_size   %||% (global_font_size * 0.80)
  fs_clust <- cluster_label_size %||% (global_font_size * 0.85)
  fs_sig   <- sig_marker_size    %||% (global_font_size * 0.90)

  # ── 3. Compute group means (variables × groups) ──────────────────────────────
  # x is observations × variables; transpose to variables × observations first
  x_t      <- t(x)    # p × n
  mean_mat <- vapply(grp_levels, function(lvl) {
    idx <- which(grp_vec == lvl)
    rowMeans(x_t[, idx, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(x_t)))                  # result: p × g
  rownames(mean_mat) <- colnames(x)
  colnames(mean_mat) <- grp_levels

  # ── 4. Row scaling ───────────────────────────────────────────────────────────
  scaled_mat <- switch(scale_vars,
    none = mean_mat,
    zscore = {
      t(apply(mean_mat, 1L, function(r) {
        sd_r <- sd(r, na.rm = TRUE)
        if (is.na(sd_r) || sd_r == 0) r - mean(r, na.rm = TRUE)
        else (r - mean(r, na.rm = TRUE)) / sd_r
      }))
    },
    minmax = {
      t(apply(mean_mat, 1L, function(r) {
        mn <- min(r, na.rm = TRUE); mx <- max(r, na.rm = TRUE)
        if (is.na(mn) || mn == mx) rep(0, length(r)) else (r - mn) / (mx - mn)
      }))
    }
  )
  rownames(scaled_mat) <- rownames(mean_mat)
  colnames(scaled_mat) <- colnames(mean_mat)

  # Remove all-NA rows before clustering
  ok_rows    <- rowSums(!is.na(scaled_mat)) > 0L
  if (!any(ok_rows))
    stop("All rows of the scaled mean matrix are entirely NA. Cannot cluster.",
         call. = FALSE)
  if (!all(ok_rows))
    warning(sprintf(
      "%d variable(s) are entirely NA after scaling and will be excluded from clustering.",
      sum(!ok_rows)
    ), call. = FALSE)
  scaled_clust <- scaled_mat[ok_rows, , drop = FALSE]

  # Impute remaining NAs with row mean for distance computation
  scaled_imp <- t(apply(scaled_clust, 1L, function(r) {
    r[is.na(r)] <- mean(r, na.rm = TRUE)
    r
  }))
  rownames(scaled_imp) <- rownames(scaled_clust)

  # ── 5. Hierarchical clustering ────────────────────────────────────────────────
  if (nrow(scaled_imp) < 2L)
    stop("At least 2 non-missing variables are required for clustering.",
         call. = FALSE)

  dist_obj <- stats::dist(scaled_imp, method = dist_method)
  hc       <- tryCatch(
    stats::hclust(dist_obj, method = clust_method),
    error = function(e) stop(
      sprintf("Hierarchical clustering failed: %s", conditionMessage(e)),
      call. = FALSE
    )
  )
  leaf_ord <- hc$order   # variable indices in leaf order (bottom → top in ggdendro)

  # ── 6. Cluster membership ────────────────────────────────────────────────────
  clust_membership <- NULL
  n_clust_actual   <- 1L
  clust_id         <- rep(1L, nrow(scaled_imp))

  if (!is.null(n_clusters)) {
    if (n_clusters > nrow(scaled_imp)) {
      warning(sprintf(
        "`n_clusters` (%d) exceeds the number of variables (%d); using %d cluster(s).",
        n_clusters, nrow(scaled_imp), nrow(scaled_imp)
      ), call. = FALSE)
      n_clusters <- nrow(scaled_imp)
    }
    clust_id         <- stats::cutree(hc, k = n_clusters)
    n_clust_actual   <- length(unique(clust_id))
    clust_membership <- clust_id
  }

  # Colour palette (Okabe-Ito, colour-blind safe)
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")
  palette_use <- if (!is.null(cluster_colors)) {
    rep_len(cluster_colors, n_clust_actual)
  } else {
    rep_len(okabe_ito, n_clust_actual)
  }

  # ── 7. Reorder matrices ───────────────────────────────────────────────────────
  var_order    <- rownames(scaled_imp)[leaf_ord]   # top → bottom in final plot
  mean_ord     <- mean_mat[var_order,   , drop = FALSE]
  scaled_ord   <- scaled_mat[var_order, , drop = FALSE]

  # ── 8. Cluster bracket geometry ──────────────────────────────────────────────
  # For each cluster, determine the contiguous span of rows in leaf order.
  # Rows are indexed 1 (bottom of dendrogram) to p (top); ggplot y-axis is
  # reversed so row 1 appears at the top of the heatmap.
  n_vars <- length(var_order)

  bracket_df <- NULL
  clust_label_df <- NULL

  if (!is.null(n_clusters)) {
    # Map each variable in leaf order to its cluster id
    clust_by_leaf <- clust_id[leaf_ord]   # length n_vars, ordered top→bottom
    # (ggdendro default: leaf 1 is at x=1 which we coord_flip to y=1=top)

    # Build one bracket per cluster: find first and last position in leaf order
    bracket_list <- lapply(seq_len(n_clust_actual), function(k) {
      pos <- which(clust_by_leaf == k)
      if (length(pos) == 0L) return(NULL)
      data.frame(
        clust    = k,
        ymin     = min(pos) - 0.5,
        ymax     = max(pos) + 0.5,
        color    = palette_use[k],
        stringsAsFactors = FALSE
      )
    })
    bracket_df <- do.call(rbind, Filter(Negate(is.null), bracket_list))

    # Auto-generate cluster labels if not supplied
    roman <- c("I","II","III","IV","V","VI","VII","VIII","IX","X")
    auto_labels <- paste0("Cluster ", roman[seq_len(n_clust_actual)])
    labels_use  <- cluster_labels %||% auto_labels

    clust_label_df <- data.frame(
      clust  = seq_len(n_clust_actual),
      y      = (bracket_df$ymin + bracket_df$ymax) / 2,
      label  = labels_use[seq_len(nrow(bracket_df))],
      color  = palette_use,
      stringsAsFactors = FALSE
    )
  }

  # ── 9. Build dendrogram panel ─────────────────────────────────────────────────
  dend_data <- ggdendro::dendro_data(hc, type = "rectangle")
  segs      <- ggdendro::segment(dend_data)
  labs      <- ggdendro::label(dend_data)

  # Colour segments by cluster when requested
  if (!is.null(n_clusters) && n_clust_actual > 1L) {
    # Each leaf position → cluster colour
    leaf_to_color <- stats::setNames(
      palette_use[clust_id[hc$order]],
      seq_along(hc$order)
    )
    segs$seg_col <- "#444444"
    leaf_mask    <- segs$yend == 0
    if (any(leaf_mask)) {
      nearest <- pmax(1L, pmin(round(segs$xend[leaf_mask]),
                                length(leaf_to_color)))
      segs$seg_col[leaf_mask] <- leaf_to_color[as.character(nearest)]
    }
  } else {
    segs$seg_col <- "#444444"
  }

  # x-axis limits for dendrogram: extend left for bracket + label margin
  dend_x_expand <- if (!is.null(n_clusters)) 0.5 else 0.05

  dend_plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data    = segs,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                              colour = seg_col),
      linewidth = 0.5
    ) +
    ggplot2::scale_colour_identity()

  # Cluster brackets: vertical coloured lines on the left of the dendrogram
  if (!is.null(bracket_df)) {
    bracket_x <- -max(segs$y, na.rm = TRUE) * 0.08   # just left of tree root
    tick_x    <-  max(segs$y, na.rm = TRUE) * 0.04

    dend_plot <- dend_plot +
      # Vertical bracket bar
      ggplot2::geom_segment(
        data    = bracket_df,
        mapping = ggplot2::aes(x = ymin, xend = ymax,
                                y = bracket_x, yend = bracket_x,
                                colour = color),
        linewidth = 2.5,
        lineend   = "butt",
        inherit.aes = FALSE
      ) +
      # Short horizontal ticks at bracket ends
      ggplot2::geom_segment(
        data    = bracket_df,
        mapping = ggplot2::aes(x = ymin, xend = ymin,
                                y = bracket_x, yend = bracket_x - tick_x,
                                colour = color),
        linewidth = 1.2, inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data    = bracket_df,
        mapping = ggplot2::aes(x = ymax, xend = ymax,
                                y = bracket_x, yend = bracket_x - tick_x,
                                colour = color),
        linewidth = 1.2, inherit.aes = FALSE
      ) +
      # Cluster labels
      ggplot2::geom_text(
        data    = clust_label_df,
        mapping = ggplot2::aes(x = y, y = bracket_x - tick_x * 2.5,
                                label = label, colour = color),
        size        = fs_clust / ggplot2::.pt,
        hjust       = 0.5,
        inherit.aes = FALSE
      )
  }

  dend_plot <- dend_plot +
    ggplot2::scale_x_continuous(
      breaks = seq_len(nrow(labs)),
      labels = labs$label,
      expand = c(0, 0.5)
    ) +
    ggplot2::scale_y_reverse(
      expand = ggplot2::expansion(mult = c(dend_x_expand, 0.05))
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_void(base_size = global_font_size) +
    ggplot2::theme(
      plot.margin     = ggplot2::margin(0, 0, 0, 0),
      axis.text       = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "none"
    )

  # ── 10. Build heatmap panel ───────────────────────────────────────────────────
  # Melt scaled_ord to long format (base R)
  nr <- nrow(scaled_ord)
  nc <- ncol(scaled_ord)
  heat_long <- data.frame(
    variable = rep(rownames(scaled_ord), times = nc),
    group    = rep(colnames(scaled_ord), each  = nr),
    value    = as.vector(scaled_ord),
    stringsAsFactors = FALSE
  )

  # Factor levels: variables in leaf order (top of plot = last level for y)
  heat_long$variable <- factor(heat_long$variable,
                                levels = rev(var_order))
  heat_long$group    <- factor(heat_long$group, levels = grp_levels)

  # Significance flag per row
  heat_long$is_sig <- heat_long$variable %in% sig_vars

  # Cell value labels
  if (show_values) {
    heat_long$label <- ifelse(
      is.na(heat_long$value), "",
      formatC(round(heat_long$value, value_digits),
              format = "f", digits = value_digits)
    )
  }

  # Colour scale midpoint
  val_range  <- range(scaled_ord, na.rm = TRUE)
  clr_mid    <- if (length(heatmap_colors) >= 3L) {
    mean(val_range)
  } else {
    NULL
  }

  heat_plot <- ggplot2::ggplot(
    heat_long,
    ggplot2::aes(x = group, y = variable, fill = value)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    {
      if (length(heatmap_colors) == 2L) {
        ggplot2::scale_fill_gradient(
          low      = heatmap_colors[1L],
          high     = heatmap_colors[2L],
          na.value = na_color,
          name     = NULL,
          guide    = ggplot2::guide_colorbar(
            barwidth  = ggplot2::unit(0.5, "lines"),
            barheight = ggplot2::unit(4,   "lines"),
            title     = NULL,
            label.theme = ggplot2::element_text(size = fs_leg)
          )
        )
      } else {
        ggplot2::scale_fill_gradientn(
          colours  = heatmap_colors,
          na.value = na_color,
          name     = NULL,
          guide    = ggplot2::guide_colorbar(
            barwidth  = ggplot2::unit(0.5, "lines"),
            barheight = ggplot2::unit(4,   "lines"),
            title     = NULL,
            label.theme = ggplot2::element_text(size = fs_leg)
          )
        )
      }
    }

  # In-cell values
  if (show_values) {
    heat_plot <- heat_plot +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        size   = (global_font_size * 0.6) / ggplot2::.pt,
        colour = value_color,
        na.rm  = TRUE
      )
  }

  # Significance markers on the right (outside the heatmap tiles)
  if (length(sig_vars) > 0L) {
    sig_df <- data.frame(
      variable = factor(sig_vars, levels = levels(heat_long$variable)),
      x_pos    = as.numeric(factor(grp_levels,
                                    levels = grp_levels))[length(grp_levels)] + 0.65,
      stringsAsFactors = FALSE
    )
    sig_df <- sig_df[!is.na(sig_df$variable), , drop = FALSE]

    heat_plot <- heat_plot +
      ggplot2::geom_text(
        data        = sig_df,
        mapping     = ggplot2::aes(x = x_pos, y = variable, label = sig_marker),
        colour      = sig_color,
        size        = fs_sig / ggplot2::.pt,
        inherit.aes = FALSE,
        hjust       = 0.5
      ) +
      # Expand x-axis right side to accommodate the marker column
      ggplot2::scale_x_discrete(
        position = "bottom",
        expand   = ggplot2::expansion(add = c(0.5, 1.2))
      )
  } else {
    heat_plot <- heat_plot +
      ggplot2::scale_x_discrete(position = "bottom",
                                  expand = ggplot2::expansion(add = c(0.5, 0.5)))
  }

  heat_plot <- heat_plot +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = global_font_size) +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text(
        angle = 0, hjust = 0.5, vjust = 1, size = fs_col
      ),
      axis.text.y    = ggplot2::element_text(
        hjust = 0, size = fs_row
      ),
      axis.title     = ggplot2::element_blank(),
      axis.ticks     = ggplot2::element_blank(),
      panel.grid     = ggplot2::element_blank(),
      panel.border   = ggplot2::element_blank(),
      plot.margin    = ggplot2::margin(0, 0, 0, 0),
      legend.position = legend_position,
      legend.text    = ggplot2::element_text(size = fs_leg)
    )

  # ── 11. Assemble ──────────────────────────────────────────────────────────────
  combined <- patchwork::wrap_plots(
    dend_plot, heat_plot,
    ncol   = 2L,
    widths = c(dend_width, 1 - dend_width)
  ) & ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  # ── 12. Attach result attributes ─────────────────────────────────────────────
  attr(combined, "mean_matrix")        <- mean_ord
  attr(combined, "scaled_matrix")      <- scaled_ord
  attr(combined, "hclust_obj")         <- hc
  attr(combined, "leaf_order")         <- leaf_ord
  attr(combined, "cluster_membership") <- clust_membership

  combined
}