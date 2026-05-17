#' @title Plot Results from a Linear Mixed-Effects Model Analysis
#'
#' @description
#' Produces a volcano plot from the \code{combined_fixed_effects} component of
#' a \code{run_mixedmodel} object, with log2 fold-change (or raw estimates) on
#' the x-axis and \eqn{-\log_{10}(p_{\text{adj}})} on the y-axis. When the
#' model formula contained more than one fixed effect term, panels are
#' automatically faceted by \code{term} so each predictor is shown separately.
#' Points are colored by significance and direction, and labels can be added to
#' the top-\eqn{n} features per panel or to a user-specified set of features.
#'
#' @param x A \code{run_mixedmodel} object returned by \code{\link{run_mixedmodel}}.
#' @param type Character string. Currently \code{"volcano"} is supported. Reserved
#'   for future types (\code{"effect"}, \code{"heatmap"}). Default: \code{"volcano"}.
#' @param term A character vector of fixed-effect term name(s) to include. When
#'   \code{NULL} (default), all non-intercept terms are included. Useful for
#'   focusing on a single predictor when the formula has many.
#' @param use_padj Logical. If \code{TRUE} (default), the y-axis reflects
#'   \code{p_adj} (Benjamini-Hochberg or whichever method was passed to
#'   \code{run_mixedmodel}). If \code{FALSE}, raw \code{p_value} is used.
#' @param p_threshold Numeric. Horizontal dashed significance line drawn at
#'   \eqn{-\log_{10}(\text{p\_threshold})}. Default: \code{0.05}.
#' @param estimate_threshold Numeric. Vertical dashed lines drawn at
#'   \eqn{\pm}\code{estimate_threshold}. Set to \code{0} to suppress. Default:
#'   \code{0} (no vertical lines).
#' @param label_features A character vector of specific feature names to label,
#'   regardless of significance rank. Combined with \code{label_top} (union).
#'   Default: \code{NULL}.
#' @param minmax A numeric vector of length 2, e.g. \code{c(-0.02, 0.02)}.
#'   Features with estimate \eqn{\le} \code{minmax[1]} or \eqn{\ge}
#'   \code{minmax[2]} are annotated. Overrides automatic IQR-based annotation.
#'   Default: \code{NULL}.
#' @param iqr_k Numeric. IQR fence multiplier used for automatic annotation
#'   when both \code{minmax} and \code{label_features} are \code{NULL}.
#'   Features beyond \eqn{Q1 - k \times IQR} or \eqn{Q3 + k \times IQR} are
#'   labelled. Set to \code{Inf} to disable automatic annotation entirely.
#'   Default: \code{1.5}.
#' @param label_size Numeric. Font size for point labels (in ggplot2 pt units).
#'   Default: \code{3}.
#' @param label_max_overlap Integer. Passed to \code{ggrepel::geom_text_repel}'s
#'   \code{max.overlaps}. Default: \code{20}.
#' @param label_box Logical. If \code{TRUE}, draws a border box around labels
#'   via \code{ggrepel::geom_label_repel}. Default: \code{FALSE}.
#' @param label_color Character. Color of label text. \code{"match"} inherits
#'   the point's significance color. Any valid color string overrides.
#'   Default: \code{"match"}.
#' @param point_size Numeric. Base point size for non-significant features.
#'   Default: \code{2}.
#' @param point_size_sig Numeric. Point size for significant features.
#'   Default: \code{2.5}.
#' @param alpha Numeric in (0, 1]. Transparency for non-significant points.
#'   Default: \code{0.6}.
#' @param label_sig_only Logical. If \code{TRUE} (default), only features
#'   passing \code{p_threshold} are eligible for \code{minmax} or IQR-based
#'   auto-annotation. Features in \code{label_features} are always annotated
#'   regardless of this setting. Default: \code{TRUE}.
#' @param color_up Character. Color for significantly up-regulated (positive
#'   estimate + significant) features. Default: \code{"#D85A30"} (Okabe-Ito
#'   vermillion).
#' @param color_down Character. Color for significantly down-regulated
#'   (negative estimate + significant) features. Default: \code{"#0072B2"}
#'   (Okabe-Ito blue).
#' @param color_ns Character. Color for non-significant features.
#'   Default: \code{"#888780"}.
#' @param x_label Character. x-axis label. Default: \code{"Estimate (fixed effect)"}.
#' @param y_label Character. y-axis label. Default:
#'   \code{expression(-log[10](p[adj]))} or \code{expression(-log[10](p))}
#'   depending on \code{use_padj}.
#' @param title Character. Plot title. Default: \code{NULL} (no title).
#' @param subtitle Character. Plot subtitle. Default: \code{NULL}.
#' @param caption Character. Plot caption. \code{NULL} produces an automatic
#'   caption reporting the p-value threshold and adjustment method.
#'   \code{""} suppresses the caption entirely. Default: \code{NULL}.
#' @param facet_scales Character. Passed to \code{ggplot2::facet_wrap}'s
#'   \code{scales} argument when multiple terms are plotted. One of
#'   \code{"fixed"}, \code{"free_x"}, \code{"free_y"}, \code{"free"}.
#'   Default: \code{"fixed"}.
#' @param facet_ncol Integer or \code{NULL}. Number of columns in the facet
#'   grid. \code{NULL} lets \code{ggplot2} decide. Default: \code{NULL}.
#' @param theme Character string. One of \code{"nature"} (default, built on
#'   \code{theme_bw}), \code{"minimal"}, \code{"classic"}, \code{"bw"}.
#' @param base_size Numeric. Base font size passed to the theme. Default:
#'   \code{12}.
#' @param font_family Character. Font family. Default: \code{""} (system
#'   default).
#' @param zoom Numeric scalar. Multiplies all point sizes, label sizes, line
#'   widths, and base font size proportionally for export scaling. Default:
#'   \code{1}.
#'
#' @details
#' \strong{Significance classes}: Each feature is assigned one of three classes
#' per term panel:
#' \itemize{
#'   \item \strong{Up} — estimate > \code{estimate_threshold} AND p
#'     (adjusted or raw) < \code{p_threshold}.
#'   \item \strong{Down} — estimate < \eqn{-}\code{estimate_threshold} AND p <
#'     \code{p_threshold}.
#'   \item \strong{NS} — does not meet either criterion.
#' }
#' When \code{estimate_threshold = 0} (default), the up/down split is purely on
#' the sign of the estimate — useful for a single continuous predictor (e.g.,
#' age, time as numeric) where there is no conventional fold-change cutoff.
#'
# \strong{Labeling strategy}: Auto-labels (\code{label_top}) are selected
# within each facet panel independently, ranked by ascending p-value. If
# \code{label_sig_only = TRUE}, only features crossing \code{p_threshold}
# are eligible. Manual labels (\code{label_features}) are always drawn even
# if non-significant. Duplicate labels (a feature in both \code{label_top}
# and \code{label_features}) are deduplicated.
#'
#' \strong{ggrepel}: If \pkg{ggrepel} is installed (strongly recommended for
#' dense plots), it is used automatically. Install with
#' \code{install.packages("ggrepel")}. If absent, \code{geom_text} is used
#' with a slight upward nudge.
#'
#' @return A \code{ggplot} object (single-term) or a patchwork-composed plot
#'   (multi-term with \code{facet_wrap}). The object can be passed to
#'   \code{save_plots()} or modified with standard \code{ggplot2} additions.
#'
#' @examples
#' \dontrun{
#' result <- run_mixedmodel(
#'   x        = x,
#'   metadata = metadata,
#'   formula  = ".feature ~ time + age + (1 | subject_id)"
#' )
#'
#' # Default: all terms, auto-label top 10 significant features
#' plot_mixedmodel(result)
#'
#' # Focus on one term, label specific features
#' plot_mixedmodel(
#'   result,
#'   term           = "timepost",
#'   label_features = c("glucose", "insulin"),
#'   color_up       = "#CC79A7",
#'   alpha          = 0.3,
#'   alpha_sig      = 0.9,
#'   title          = "Time effect (pre → post)"
#' )
#'
#' # No auto-labels, only highlight manually specified features
#' plot_mixedmodel(
#'   result,
#'   label_features = c("glucose", "insulin", "PA O-28:1"),
#'   label_box      = TRUE,
#'   label_color    = "black"
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text
#'   scale_color_manual scale_alpha_manual scale_size_manual labs facet_wrap
#'   theme_bw theme_minimal theme_classic theme element_text element_rect
#'   element_line element_blank margin unit
#' @importFrom stats setNames
#'
#' @export
plot_mixedmodel <- function(
    x,
    type               = "volcano",
    term               = NULL,
    use_padj           = TRUE,
    p_threshold        = 0.05,
    estimate_threshold = 0,
    label_features     = NULL,
    minmax             = NULL,
    iqr_k              = 1.5,
    label_size         = 3,
    label_max_overlap  = 20L,
    label_box          = FALSE,
    label_color        = "match",
    point_size         = 2,
    point_size_sig     = 2.5,
    alpha              = 0.6,
    label_sig_only     = TRUE,
    color_up           = "#D85A30",
    color_down         = "#0072B2",
    color_ns           = "#888780",
    x_label            = "Estimate (fixed effect)",
    y_label            = NULL,
    title              = NULL,
    subtitle           = NULL,
    caption            = NULL,
    facet_scales       = "fixed",
    facet_ncol         = NULL,
    theme              = "nature",
    base_size          = 12,
    font_family        = "",
    zoom               = 1
) {

  # ---- soft dependency: ggplot2 --------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")

  # ---- input class ---------------------------------------------------------
  if (!inherits(x, "run_mixedmodel"))
    stop("'x' must be a 'run_mixedmodel' object returned by run_mixedmodel().")

  # ---- type ----------------------------------------------------------------
  type <- match.arg(type, choices = "volcano")   # expand as types are added

  # ---- minmax ----------------------------------------------------------------
  if (!is.null(minmax)) {
    if (!is.numeric(minmax) || length(minmax) != 2L || anyNA(minmax))
      stop("'minmax' must be a numeric vector of length 2, e.g. c(-0.02, 0.02).")
    if (minmax[1L] >= minmax[2L])
      stop("'minmax[1]' must be strictly less than 'minmax[2]'.")
  }

  # ---- iqr_k ----------------------------------------------------------------
  if (!is.numeric(iqr_k) || length(iqr_k) != 1L || iqr_k < 0)
    stop("'iqr_k' must be a single non-negative numeric value.")

  # ---- zoom ----------------------------------------------------------------
  if (!is.numeric(zoom) || length(zoom) != 1L || zoom <= 0)
    stop("'zoom' must be a single positive numeric.")
  base_size      <- base_size      * zoom
  label_size     <- label_size     * zoom
  point_size     <- point_size     * zoom
  point_size_sig <- point_size_sig * zoom

  # ---- combined_fixed_effects ----------------------------------------------
  cfe <- x$combined_fixed_effects
  if (is.null(cfe) || !is.data.frame(cfe) || nrow(cfe) == 0L)
    stop(
      "No combined_fixed_effects table found. ",
      "Ensure at least one model converged in run_mixedmodel()."
    )

  required_cols <- c("feature", "term", "estimate", "p_value")
  missing_cols  <- setdiff(required_cols, names(cfe))
  if (length(missing_cols))
    stop(sprintf(
      "combined_fixed_effects is missing required columns: %s.",
      paste(missing_cols, collapse = ", ")
    ))

  # ---- p column selection --------------------------------------------------
  p_col <- if (use_padj && "p_adj" %in% names(cfe)) "p_adj" else "p_value"
  if (use_padj && p_col == "p_value")
    warning(
      "'use_padj = TRUE' but 'p_adj' column not found in combined_fixed_effects; ",
      "falling back to 'p_value'."
    )

  # ---- filter to requested terms, drop intercept ---------------------------
  cfe <- cfe[cfe$term != "(Intercept)", , drop = FALSE]
  if (nrow(cfe) == 0L)
    stop("No non-intercept terms remain in combined_fixed_effects after removing '(Intercept)'.")

  if (!is.null(term)) {
    if (!is.character(term))
      stop("'term' must be a character vector of fixed-effect term names.")
    bad_terms <- setdiff(term, unique(cfe$term))
    if (length(bad_terms))
      warning(sprintf(
        "The following terms were not found in combined_fixed_effects and will be ignored: %s.\n  Available terms: %s.",
        paste(bad_terms, collapse = ", "),
        paste(unique(cfe$term), collapse = ", ")
      ))
    cfe <- cfe[cfe$term %in% term, , drop = FALSE]
    if (nrow(cfe) == 0L)
      stop("No rows remain after filtering to the requested term(s).")
  }

  # ---- validate scalar numerics --------------------------------------------
  .check_scalar_num <- function(val, nm) {
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val))
      stop(sprintf("'%s' must be a single finite numeric value.", nm))
  }
  .check_scalar_num(p_threshold, "p_threshold")
  .check_scalar_num(estimate_threshold, "estimate_threshold")
  .check_scalar_num(alpha, "alpha")
  if (p_threshold <= 0 || p_threshold >= 1)
    stop("'p_threshold' must be strictly between 0 and 1.")

  # ---- build plot data -----------------------------------------------------
  p_vals <- cfe[[p_col]]
  p_vals[p_vals <= 0] <- .Machine$double.eps  # guard against log(0)

  cfe$.neg_log10_p <- -log10(p_vals)
  cfe$.sig <- ifelse(
    !is.finite(cfe$.neg_log10_p) | is.na(p_vals),
    "NS",
    ifelse(
      p_vals < p_threshold & cfe$estimate >  estimate_threshold, "Up",
      ifelse(
        p_vals < p_threshold & cfe$estimate < -estimate_threshold, "Down",
        "NS"
      )
    )
  )

  # If estimate_threshold == 0, sign alone determines direction
  if (estimate_threshold == 0) {
    up_mask   <- p_vals < p_threshold & cfe$estimate >= 0
    down_mask <- p_vals < p_threshold & cfe$estimate <  0
    cfe$.sig  <- ifelse(up_mask, "Up", ifelse(down_mask, "Down", "NS"))
  }

  cfe$.sig <- factor(cfe$.sig, levels = c("Up", "Down", "NS"))

  # ---- annotation selection (by estimate extremity) -----------------------
  label_features <- if (!is.null(label_features)) as.character(label_features) else character(0)

  terms_in_plot <- unique(as.character(cfe$term))
  label_rows    <- logical(nrow(cfe))

  for (trm in terms_in_plot) {
    idx_term <- which(cfe$term == trm)
    ests     <- cfe$estimate[idx_term]
    p_vals_t <- cfe[[p_col]][idx_term]

    # Significance mask for this panel (used when label_sig_only = TRUE)
    sig_mask <- !is.na(p_vals_t) & p_vals_t < p_threshold

    # 1. Manual label_features — always annotated regardless of label_sig_only
    manual_mask <- cfe$feature[idx_term] %in% label_features
    label_rows[idx_term[manual_mask]] <- TRUE

    # Candidate pool respects label_sig_only (excludes manual — already handled)
    candidate_idx <- if (label_sig_only) idx_term[sig_mask] else idx_term
    cand_ests     <- cfe$estimate[candidate_idx]

    # 2. minmax — explicit thresholds on the estimate
    if (!is.null(minmax)) {
      mm_mask <- cand_ests <= minmax[1L] | cand_ests >= minmax[2L]
      label_rows[candidate_idx[mm_mask]] <- TRUE

    # 3. Automatic IQR fence (only when minmax is NULL and no label_features given)
    } else if (length(label_features) == 0L && is.finite(iqr_k)) {
      if (length(cand_ests) > 0L) {
        q1  <- stats::quantile(cand_ests, 0.25, na.rm = TRUE)
        q3  <- stats::quantile(cand_ests, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        if (iqr > 0) {
          fence_mask <- cand_ests < (q1 - iqr_k * iqr) | cand_ests > (q3 + iqr_k * iqr)
          label_rows[candidate_idx[fence_mask]] <- TRUE
        }
      }
    }
  }

  cfe$.label     <- ifelse(label_rows, as.character(cfe$feature), "")
  cfe$.lbl_color <- ifelse(
    cfe$.sig == "Up",   color_up,
    ifelse(cfe$.sig == "Down", color_down, color_ns)
  )

  # ---- y-axis label --------------------------------------------------------
  if (is.null(y_label)) {
    y_label <- if (p_col == "p_adj")
      expression(-log[10](p[adj]))
    else
      expression(-log[10](p))
  }

  # ---- auto caption --------------------------------------------------------
  if (is.null(caption)) {
    adj_method <- x$parameters$p_adjust_method %||% "BH"
    p_src      <- if (p_col == "p_adj") paste0("p.adj (", adj_method, ")") else "p (unadjusted)"
    caption    <- sprintf(
      "Significance threshold: %s < %.3g | Dashed line at -log\u2081\u2080(%.3g) = %.2f",
      p_src, p_threshold, p_threshold, -log10(p_threshold)
    )
  }

  # ---- theme builder -------------------------------------------------------
  .build_theme <- function(theme_nm, bs, ff) {
    base <- switch(
      theme_nm,
      nature  = ggplot2::theme_bw(base_size = bs, base_family = ff),
      minimal = ggplot2::theme_minimal(base_size = bs, base_family = ff),
      classic = ggplot2::theme_classic(base_size = bs, base_family = ff),
      bw      = ggplot2::theme_bw(base_size = bs, base_family = ff),
      stop(sprintf(
        "'theme' must be one of: 'nature', 'minimal', 'classic', 'bw'. Got: '%s'.", theme_nm
      ))
    )
    if (theme_nm == "nature") {
      base <- base + ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
        strip.text       = ggplot2::element_text(face = "bold", size = bs * 0.9),
        axis.title       = ggplot2::element_text(size = bs),
        axis.text        = ggplot2::element_text(size = bs * 0.85),
        legend.position  = "none",
        plot.caption     = ggplot2::element_text(
          size = bs * 0.72, color = "grey50", hjust = 0,
          margin = ggplot2::margin(t = 6)
        ),
        plot.title       = ggplot2::element_text(face = "bold", size = bs * 1.1),
        plot.subtitle    = ggplot2::element_text(size = bs * 0.9, color = "grey40")
      )
    }
    base
  }

  thm <- .build_theme(theme, base_size, font_family)

  # ---- color / alpha / size scales -----------------------------------------
  sig_colors <- c(Up = color_up, Down = color_down, NS = color_ns)
  sig_alphas <- c(Up = alpha, Down = alpha, NS = alpha)
  sig_sizes  <- c(Up = point_size_sig, Down = point_size_sig, NS = point_size)

  # ---- build base plot -----------------------------------------------------
  p <- ggplot2::ggplot(
    cfe,
    ggplot2::aes(
      x     = estimate,
      y     = .neg_log10_p,
      color = .sig,
      alpha = .sig,
      size  = .sig
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_threshold),
      linetype   = "dashed",
      linewidth  = 0.45 * zoom,
      color      = "grey40"
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = sig_colors, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = sig_alphas, drop = FALSE) +
    ggplot2::scale_size_manual(values  = sig_sizes,  drop = FALSE) +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = title,
      subtitle = subtitle,
      caption  = if (nzchar(caption)) caption else NULL
    ) +
    thm

  # Vertical threshold lines (only if estimate_threshold > 0)
  if (estimate_threshold > 0) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = c(-estimate_threshold, estimate_threshold),
        linetype   = "dashed",
        linewidth  = 0.45 * zoom,
        color      = "grey40"
      )
  }

  # ---- facet (multi-term) --------------------------------------------------
  if (length(terms_in_plot) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~ term,
      scales = facet_scales,
      ncol   = facet_ncol
    )
  }

  # ---- labels --------------------------------------------------------------
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  label_df    <- cfe[cfe$.label != "", , drop = FALSE]

  if (nrow(label_df) > 0L) {

    if (has_ggrepel) {
      repel_fn <- if (label_box) ggrepel::geom_label_repel else ggrepel::geom_text_repel

      # Resolve label color mapping ----------------------------------------
      if (identical(label_color, "match")) {
        # Per-point color lives inside aes(); no bare color argument
        p <- p + repel_fn(
          data         = label_df,
          mapping      = ggplot2::aes(
            x     = estimate,
            y     = .neg_log10_p,
            label = .label,
            color = .sig          # <-- mapped, not manual
          ),
          size          = label_size,
          max.overlaps  = label_max_overlap,
          segment.size  = 0.3 * zoom,
          segment.color = "grey50",
          show.legend   = FALSE,
          inherit.aes   = FALSE
        ) +
          # Reuse the same color scale already on the plot
          ggplot2::scale_color_manual(values = sig_colors, drop = FALSE)
      } else {
        # Single fixed color — safe to pass as bare argument
        p <- p + repel_fn(
          data         = label_df,
          mapping      = ggplot2::aes(
            x     = estimate,
            y     = .neg_log10_p,
            label = .label
          ),
          color         = label_color,
          size          = label_size,
          max.overlaps  = label_max_overlap,
          segment.size  = 0.3 * zoom,
          segment.color = "grey50",
          show.legend   = FALSE,
          inherit.aes   = FALSE
        )
      }

    } else {
      message(
        "Install 'ggrepel' for better label placement: install.packages('ggrepel').\n",
        "Falling back to geom_text."
      )

      lbl_color_vec <- if (identical(label_color, "match")) label_df$.lbl_color else label_color

      p <- p + ggplot2::geom_text(
        data        = label_df,
        mapping     = ggplot2::aes(
          x     = estimate,
          y     = .neg_log10_p,
          label = .label
        ),
        color       = lbl_color_vec,
        size        = label_size,
        nudge_y     = 0.15,
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
  }

  p   # visible return: auto-prints when unassigned, silent when assigned
}