# =============================================================================
# plot_heatmap.R
# Part of the pondeR package
# =============================================================================

#' @title Plot a Correlation Heatmap
#'
#' @description
#' Produces a publication-ready correlation heatmap with optional hierarchical
#' clustering dendrograms, significance masking, triangle display, and
#' in-cell coefficient annotations. Internally calls \code{\link{run_correl}}
#' and returns a \pkg{ggplot2} object.
#'
#' @param x A data frame, tibble, or matrix. Column names may contain special
#'   characters.
#' @param method Correlation method passed to \code{\link{run_correl}}. One of
#'   \code{"auto"} (default), \code{"pearson"}, \code{"spearman"},
#'   \code{"kendall"}, \code{"pointbiserial"}, \code{"phi"}, \code{"cramerv"},
#'   or \code{"poly"}.
#' @param out A single character string: \code{"all"} (default; plot all
#'   correlations) or \code{"masked"} (set non-significant cells to
#'   \code{NA}).
#' @param sig_threshold Significance threshold passed to
#'   \code{\link{run_correl}}. Default \code{0.05}.
#' @param remove Level-exclusion list passed to \code{\link{run_correl}}.
#' @param theme A single character string: \code{"nature"} (default),
#'   \code{"classic"}, \code{"minimal"}, or \code{"bw"}.
#' @param plot_title A single character string for the plot title. Default
#'   \code{NULL}.
#' @param plot_subtitle A single character string for the plot subtitle. Default
#'   \code{NULL}.
#' @param show_rlab Logical. Show row labels. Default \code{TRUE}.
#' @param show_clab Logical. Show column labels. Default \code{TRUE}.
#' @param top_n A single positive integer. Plot only the \code{top_n} pairs
#'   with the highest absolute correlation. Mutually exclusive with
#'   \code{bottom_n}. Default \code{NULL}.
#' @param bottom_n A single positive integer. Plot only the \code{bottom_n}
#'   pairs with the lowest absolute correlation. Mutually exclusive with
#'   \code{top_n}. Default \code{NULL}.
#' @param color A character vector of length 3 specifying fill colours for
#'   negative, neutral, and positive correlation respectively. Default
#'   \code{c("#2166AC", "#F7F7F7", "#B2182B")}.
#' @param clab_angle A single numeric (0–360) controlling the rotation angle
#'   of column (x-axis) labels in degrees. Default \code{45}.
#' @param show_triangle A single character string: \code{"all"} (default),
#'   \code{"lower"} (retain lower triangle only), or \code{"upper"} (retain
#'   upper triangle only).
#' @param show_coef Logical. Annotate each cell with its rounded correlation
#'   coefficient. Default \code{FALSE}.
#' @param round A single non-negative integer. Number of decimal places for
#'   coefficient annotation. Default \code{2}. When rounding produces exactly
#'   \code{1} or \code{-1} for a non-unity value, scientific notation is used
#'   and the font is reduced. Ignored when \code{show_coef = FALSE}.
#' @param mark_x Logical. Overlay an \code{"\u00D7"} on non-significant cells
#'   when \code{show_coef = TRUE}. Default \code{TRUE}.
#' @param mark_x_size A single positive numeric acting as a multiplier for the 
#'   size of the "×" mark. Default \code{1.0} (fills ~85% of the cell). 
#'   Ignored when \code{mark_x = FALSE}.
#' @param mark_x_alpha A single numeric (0–1) controlling the transparency of 
#'   the "×" mark. Default \code{0.6}. Ignored when \code{mark_x = FALSE}.
#' @param show_dend Logical. Whether to attach hierarchical clustering
#'   dendrograms to the heatmap margins. Requires \pkg{patchwork}.
#'   Default \code{FALSE}.
#' @param dend_row_width A single positive numeric controlling the relative
#'   width of the row dendrogram panel when \code{show_dend = TRUE}.
#'   Default \code{0.2}.
#' @param dend_col_height A single positive numeric controlling the relative
#'   height of the column dendrogram panel when \code{show_dend = TRUE}.
#'   Default \code{0.2}.
#' @param global_font_size A single positive numeric specifying the base font
#'   size in points. Default \code{15}.
#' @param axis_title_size Font size for axis titles. Defaults to
#'   \code{global_font_size}.
#' @param axis_text_size Font size for axis tick labels. Defaults to
#'   \code{global_font_size}.
#' @param legend_title_size Font size for the legend title. Defaults to
#'   \code{global_font_size}.
#' @param legend_text_size Font size for legend text. Defaults to
#'   \code{global_font_size}.
#' @param coef_text_size Font size for in-cell coefficient annotations.
#'   Defaults to \code{global_font_size * 0.7}. Ignored when
#'   \code{show_coef = FALSE}.
#' @param plot_title_size Font size for the plot title. Defaults to
#'   \code{global_font_size + 2}.
#' @param plot_subtitle_size Font size for the plot subtitle. Defaults to
#'   \code{global_font_size}.
#' @param ... Additional arguments passed to \code{\link{run_correl}} and,
#'   when \code{show_dend = TRUE}, to \code{\link{plot_dend}}.
#'
#' @details
#' Hierarchical clustering uses Ward's D2 linkage on \code{1 - |r|} as the
#' dissimilarity to reorder variables before plotting.
#'
#' When \code{show_dend = TRUE}, dendrograms are composed with the heatmap
#' using \pkg{patchwork}. The row dendrogram is placed to the left and the
#' column dendrogram above the heatmap. Arguments accepted by
#' \code{\link{plot_dend}} (e.g. \code{agglomeration_method},
#' \code{line_color}) are forwarded automatically via \code{...}.
#'
#' @return A \code{ggplot2} object, or a \pkg{patchwork} composite when
#'   \code{show_dend = TRUE}.
#'
#' @examples
#' ## Basic heatmap
#' plot_heatmap(mtcars)
#'
#' ## Masked, lower triangle only
#' plot_heatmap(mtcars, out = "masked", show_triangle = "lower")
#'
#' ## Show coefficients
#' plot_heatmap(mtcars, show_coef = TRUE, mark_x = TRUE)
#'
#' ## With dendrograms
#' plot_heatmap(mtcars, show_dend = TRUE)
#'
#' @author John Lennon L. Calorio
#' @export
plot_heatmap <- function(x,
                          method             = "auto",
                          out                = "all",
                          sig_threshold      = 0.05,
                          remove             = NULL,
                          theme              = "nature",
                          plot_title         = NULL,
                          plot_subtitle      = NULL,
                          show_rlab          = TRUE,
                          show_clab          = TRUE,
                          top_n              = NULL,
                          bottom_n           = NULL,
                          color              = c("#2166AC", "#F7F7F7", "#B2182B"),
                          clab_angle         = 45,
                          show_triangle      = "all",
                          show_coef          = FALSE,
                          round              = 2L,
                          mark_x             = TRUE,
                          mark_x_size        = 1.0,
                          mark_x_alpha       = 0.6,
                          show_dend          = FALSE,
                          dend_row_width     = 0.2,
                          dend_col_height    = 0.2,
                          global_font_size   = 15,
                          axis_title_size    = NULL,
                          axis_text_size     = NULL,
                          legend_title_size  = NULL,
                          legend_text_size   = NULL,
                          coef_text_size     = NULL,
                          plot_title_size    = NULL,
                          plot_subtitle_size = NULL,
                          ...) {

  # --- Dependencies -----------------------------------------------------------
  if (!requireNamespace("ggplot2",  quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')",
         call. = FALSE)
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required. Install with: install.packages('reshape2')",
         call. = FALSE)
  if (show_dend && !requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for `show_dend = TRUE`. ",
         "Install with: install.packages('patchwork')", call. = FALSE)

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x))
    stop("`x` must be a data frame, tibble, or matrix.", call. = FALSE)
  if (is.matrix(x)) x <- as.data.frame(x)

  .check_str1 <- function(param, valid, nm) {
    if (!is.character(param) || length(param) != 1L || !param %in% valid)
      stop("`", nm, "` must be one of: ",
           paste0('"', valid, '"', collapse = ", "), ".", call. = FALSE)
  }
  .check_str1(out,           c("all", "masked"),                  "out")
  .check_str1(show_triangle, c("all", "lower", "upper"),          "show_triangle")
  .check_str1(theme,         c("nature", "classic", "minimal", "bw"), "theme")

  if (!is.numeric(mark_x_size) || length(mark_x_size) != 1L || mark_x_size <= 0)
    stop("`mark_x_size` must be a single positive numeric.", call. = FALSE)
  if (!is.numeric(mark_x_alpha) || length(mark_x_alpha) != 1L || 
      mark_x_alpha < 0 || mark_x_alpha > 1)
    stop("`mark_x_alpha` must be a numeric between 0 and 1.", call. = FALSE)
  
  if (!is.null(plot_title) && (!is.character(plot_title) || length(plot_title) != 1L))
    stop("`plot_title` must be a single character string or NULL.", call. = FALSE)
  if (!is.null(plot_subtitle) && (!is.character(plot_subtitle) || length(plot_subtitle) != 1L))
    stop("`plot_subtitle` must be a single character string or NULL.", call. = FALSE)

  if (!is.character(color) || length(color) != 3L)
    stop("`color` must be a character vector of length 3.", call. = FALSE)

  if (!is.numeric(clab_angle) || length(clab_angle) != 1L ||
      clab_angle < 0 || clab_angle > 360)
    stop("`clab_angle` must be a single numeric between 0 and 360.", call. = FALSE)

  if (!is.null(top_n) && !is.null(bottom_n))
    stop("`top_n` and `bottom_n` cannot be specified simultaneously.", call. = FALSE)
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1L)
      stop("`top_n` must be a single positive integer.", call. = FALSE)
    top_n <- as.integer(top_n)
  }
  if (!is.null(bottom_n)) {
    if (!is.numeric(bottom_n) || length(bottom_n) != 1L || bottom_n < 1L)
      stop("`bottom_n` must be a single positive integer.", call. = FALSE)
    bottom_n <- as.integer(bottom_n)
  }

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
  fs_coef          <- .fs(coef_text_size,     global_font_size * 0.7)
  fs_plot_title    <- .fs(plot_title_size,    global_font_size + 2)
  fs_plot_subtitle <- .fs(plot_subtitle_size, global_font_size)

  # --- Separate ... args for run_correl vs plot_dend --------------------------
  dots          <- list(...)
  dend_formals  <- c(names(formals(plot_dend)),
                     names(formals(stats::dist)),
                     names(formals(stats::hclust)))
  correl_formals <- names(formals(run_correl))

  dend_args   <- dots[names(dots) %in% dend_formals]
  correl_args <- dots[names(dots) %in% correl_formals]

  # --- Run correlation --------------------------------------------------------
  cor_res <- do.call(
    run_correl,
    c(list(x = x, remove = remove, method = method,
           sig_threshold = sig_threshold),
      correl_args)
  )

  tbl_use  <- if (out == "masked") cor_res$correlations_masked else cor_res$correlations
  pval_tbl <- cor_res$p_values

  # --- top_n / bottom_n filter ------------------------------------------------
  if (!is.null(top_n) || !is.null(bottom_n)) {
    abs_est <- abs(cor_res$correlations$estimate)
    keep_idx <- if (!is.null(top_n)) {
      order(abs_est, decreasing = TRUE)[seq_len(min(top_n, length(abs_est)))]
    } else {
      order(abs_est, decreasing = FALSE)[seq_len(min(bottom_n, length(abs_est)))]
    }
    vars_keep <- unique(c(tbl_use$var1[keep_idx], tbl_use$var2[keep_idx]))
    tbl_use  <- tbl_use[tbl_use$var1 %in% vars_keep & tbl_use$var2 %in% vars_keep, ]
    pval_tbl <- pval_tbl[pval_tbl$var1 %in% vars_keep & pval_tbl$var2 %in% vars_keep, ]
  }

  # --- Build symmetric matrix -------------------------------------------------
  all_vars <- unique(c(tbl_use$var1, tbl_use$var2))
  nv       <- length(all_vars)

  mat <- matrix(NA_real_, nv, nv, dimnames = list(all_vars, all_vars))
  diag(mat) <- 1
  for (i in seq_len(nrow(tbl_use))) {
    v1 <- tbl_use$var1[i]; v2 <- tbl_use$var2[i]; est <- tbl_use$estimate[i]
    mat[v1, v2] <- est
    mat[v2, v1] <- est
  }

  pmat <- matrix(NA_real_, nv, nv, dimnames = list(all_vars, all_vars))
  diag(pmat) <- 0
  for (i in seq_len(nrow(pval_tbl))) {
    v1 <- pval_tbl$var1[i]; v2 <- pval_tbl$var2[i]
    if (v1 %in% all_vars && v2 %in% all_vars) {
      pmat[v1, v2] <- pval_tbl$p_value[i]
      pmat[v2, v1] <- pval_tbl$p_value[i]
    }
  }

  # --- Cluster and reorder ----------------------------------------------------
  col_ok    <- colSums(!is.na(mat)) > 0L
  mat_clean <- mat[col_ok, col_ok, drop = FALSE]

  if (nrow(mat_clean) > 1L) {
    dist_mat <- stats::dist(1 - abs(mat_clean), method = "euclidean")
    hc       <- stats::hclust(dist_mat, method = "ward.D2")
    ord      <- hc$labels[hc$order]
  } else {
    ord <- rownames(mat_clean)
  }
  mat_ord  <- mat[ord, ord]
  pmat_ord <- pmat[ord, ord]

  # --- Apply triangle mask ----------------------------------------------------
  # "lower" keeps the lower triangle; upper.tri positions are set to NA
  # "upper" keeps the upper triangle; lower.tri positions are set to NA
  if (show_triangle == "lower") {
    mat_ord[upper.tri(mat_ord, diag = FALSE)] <- NA_real_
  } else if (show_triangle == "upper") {
    mat_ord[lower.tri(mat_ord, diag = FALSE)] <- NA_real_
  }

  # --- Melt to long -----------------------------------------------------------
  df_long <- reshape2::melt(mat_ord, varnames = c("Var1", "Var2"),
                            value.name = "estimate", na.rm = FALSE)
  df_long$Var1 <- factor(df_long$Var1, levels = ord)
  df_long$Var2 <- factor(df_long$Var2, levels = rev(ord))

  pmat_long    <- reshape2::melt(pmat_ord, varnames = c("Var1", "Var2"),
                                 value.name = "p_value", na.rm = FALSE)
  df_long$p_value <- pmat_long$p_value[
    match(paste(df_long$Var1, df_long$Var2),
          paste(pmat_long$Var1, pmat_long$Var2))
  ]

  # --- Coefficient labels -----------------------------------------------------
  if (show_coef) {
    round_val     <- max(0L, as.integer(round))
    df_long$label <- vapply(df_long$estimate, function(v) {
      if (is.na(v)) return("")
      rv <- base::round(v, round_val)
      if (abs(rv) == 1 && abs(v) != 1) formatC(v, format = "e", digits = 2)
      else formatC(rv, format = "f", digits = round_val)
    }, character(1L))
    df_long$use_small <- !is.na(df_long$estimate) &
      abs(base::round(df_long$estimate, round_val)) == 1 &
      abs(df_long$estimate) != 1
  }

  # --- Build ggplot -----------------------------------------------------------
  # Determine axis label hjust based on angle
  clab_hjust <- if (clab_angle == 0) 0.5 else 1

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Var1, y = Var2, fill = estimate)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4, na.rm = TRUE) +
    ggplot2::scale_fill_gradient2(
      low      = color[1L],
      mid      = color[2L],
      high     = color[3L],
      midpoint = 0,
      limits   = c(-1, 1),
      na.value = "grey95",
      name     = "r"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x        = NULL,
      y        = NULL,
      title    = plot_title,
      subtitle = plot_subtitle
    )

  # Coefficient annotations
  if (show_coef) {
    df_coef        <- df_long[!is.na(df_long$estimate) & df_long$label != "", ]
    df_coef_normal <- df_coef[!df_coef$use_small, ]
    df_coef_small  <- df_coef[ df_coef$use_small,  ]

    if (nrow(df_coef_normal) > 0L)
      p <- p + ggplot2::geom_text(
        data = df_coef_normal,
        ggplot2::aes(label = label),
        size   = fs_coef / ggplot2::.pt,
        colour = "black"
      )
    if (nrow(df_coef_small) > 0L)
      p <- p + ggplot2::geom_text(
        data = df_coef_small,
        ggplot2::aes(label = label),
        size   = (fs_coef * 0.7) / ggplot2::.pt,
        colour = "black"
      )

    if (show_coef && mark_x) {
    df_ns <- df_long[!is.na(df_long$p_value) & 
                     df_long$p_value >= sig_threshold & 
                     !is.na(df_long$estimate), ]
    
    if (nrow(df_ns) > 0L) {
      # Calculate base size to fill the tile
      n_vars      <- length(levels(df_long$Var1))
      base_x_size <- (150 / n_vars) * 0.85
      
      # Apply the user-defined multiplier
      final_x_size <- base_x_size * mark_x_size

      p <- p + ggplot2::geom_text(
        data = df_ns,
        ggplot2::aes(label = "\u00D7"),
        size   = final_x_size,
        colour = "grey40",
        alpha  = mark_x_alpha # Uses the new transparency parameter
      )
    }
  }
  }

  # Base theme — call the function, not the object, so base_size is applied
  base_theme_fn <- switch(theme,
    nature  = ggplot2::theme_bw,
    classic = ggplot2::theme_classic,
    minimal = ggplot2::theme_minimal,
    bw      = ggplot2::theme_bw
  )

  p <- p +
    base_theme_fn(base_size = global_font_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.title.x     = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_blank(),
      axis.text.x = if (show_clab)
        ggplot2::element_text(angle = clab_angle, hjust = clab_hjust,
                              size  = fs_axis_text)
      else
        ggplot2::element_blank(),
      axis.text.y = if (show_rlab)
        ggplot2::element_text(size = fs_axis_text)
      else
        ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(size = fs_legend_title),
      legend.text   = ggplot2::element_text(size = fs_legend_text),
      plot.title    = ggplot2::element_text(size = fs_plot_title,    face = "bold"),
      plot.subtitle = ggplot2::element_text(size = fs_plot_subtitle)
    )

  # --- Attach dendrograms if requested ----------------------------------------
  if (!show_dend) return(p)

  dend_both <- do.call(
    plot_dend,
    c(list(x           = x,
           orientation = "both",
           theme       = theme,
           global_font_size = global_font_size),
      dend_args)
  )

  # Strip titles/subtitles from dendrograms (they live on the heatmap)
  dend_row <- dend_both$rows +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_blank())
  dend_col <- dend_both$cols +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_blank())

  # Compose with patchwork:
  #   [col dend  ]
  #   [row dend | heatmap]
  empty_panel <- ggplot2::ggplot() + ggplot2::theme_void()

  top_row    <- patchwork::wrap_plots(empty_panel, dend_col,
                                      widths = c(dend_row_width, 1))
  bottom_row <- patchwork::wrap_plots(dend_row, p,
                                      widths = c(dend_row_width, 1))

  patchwork::wrap_plots(top_row, bottom_row,
                        ncol    = 1,
                        heights = c(dend_col_height, 1))
}