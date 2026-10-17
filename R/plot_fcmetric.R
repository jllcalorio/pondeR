# =============================================================================
#  plot_fcmetric
# =============================================================================

#' @title Fold Change vs. Feature Metric Plot (MA-style)
#'
#' @description
#' Plots the fold change of each feature (y-axis) against a per-feature
#' metric (x-axis): mean abundance (the classic MA plot), effect size,
#' variance, standard deviation, or VIP score. Unlike
#' \code{\link{plot_volcano}()} it needs no p-values, so it remains useful
#' when \code{run_diff()} found nothing significant. With
#' \code{metric = "abundance"} it shows whether large fold changes come from
#' low-intensity "noise" features whose baseline is near zero. Styling
#' arguments match \code{\link{plot_volcano}()}.
#'
#' @param x A \code{run_foldchange()} result, or a named list holding it plus
#'   the object needed by \code{metric}, e.g.
#'   \code{list(fc = fc, diff = diff, pls = pls)}; \strong{or} a
#'   \code{data.frame}/\code{tibble}/named \code{matrix} with one row per
#'   feature. One plot is made per comparison.
#'
#'   The data stage (normalized, transformed, scaled, ...) is chosen upstream,
#'   e.g. \code{run_foldchange(from = "normalized")}. For a \code{data.frame},
#'   supply columns computed from whichever stage you want.
#' @param y A single character string naming the column of \strong{raw} fold
#'   changes (positive numeric). Default \code{NULL} (\code{"fold_change"} for
#'   \code{run_foldchange()} input).
#' @param z A single character string naming the column holding the x-axis
#'   metric. Default \code{NULL} (derived from \code{x} and \code{metric} for
#'   \code{run_foldchange()} input; required for a \code{data.frame}, where
#'   \code{metric} then only sets the axis label and default scale).
#' @param group An optional column name in \code{x}; one plot is made per level.
#'   Default \code{NULL}.
#' @param features A single character string naming the feature identifier
#'   column. Required when \code{annotate} is used. Default \code{NULL}.
#' @param metric A single string, the x-axis metric. See \emph{Metrics}.
#'   Default \code{"abundance"}.
#' @param log2 Logical. If \code{TRUE} (default) the y-axis is
#'   \eqn{\log_2}(fold change); if \code{FALSE} it is the raw fold change.
#' @param log10_x Logical or \code{NULL}. If \code{TRUE} the x-axis is
#'   \eqn{\log_{10}} of the metric. Default \code{NULL} chooses automatically:
#'   \code{TRUE} for \code{"abundance"}, \code{"variance"} and \code{"sd"};
#'   \code{FALSE} for \code{"effect_size"} and \code{"vip"}. If the metric has
#'   non-positive values it is set to \code{FALSE} with a warning.
#' @param x_cut Optional numeric vector of reference values (on the original
#'   metric scale, not log10) drawn as vertical dashed lines, e.g.
#'   \code{1} for VIP or \code{c(0.2, 0.5, 0.8)} for Cohen's d. Default
#'   \code{NULL}.
#' @param up,down Raw fold-change thresholds (\code{up > 1}, \code{0 < down < 1}).
#'   Features beyond them are coloured Up/Down. Inherited from a
#'   \code{run_foldchange()} object when missing; otherwise \code{1.5}/\code{0.5}.
#' @param annotate Optional (named) character vector of feature identifiers
#'   to label; see \code{\link{plot_volcano}}. Default \code{NULL}.
#' @param use_ggrepel,up_color,down_color,ns_color,h_color,v_color,size,alpha,theme,plot_title,plot_subtitle,xlab,ylab,global_font_size,title_size,subtitle_size,xlab_size,ylab_size,axis_text_size,label_size,legend_size,legend_title,legend_nrow,...
#'   As in \code{\link{plot_volcano}}. Here \code{h_color} colours the
#'   horizontal \code{up}/\code{down} fold-change lines and \code{v_color} the
#'   vertical \code{x_cut} lines. \code{plot_title = NULL} gives a title naming
#'   the metric (plus the comparison when there are several).
#'
#' @section Metrics:
#' \describe{
#'   \item{\code{"abundance"}}{Mean of the two group means of each comparison.
#'     Needs only \code{run_foldchange()}.}
#'   \item{\code{"variance"}, \code{"sd"}}{Pooled within-group variance / SD of
#'     the two groups in each comparison,
#'     \eqn{s_p^2 = ((n_1-1)s_1^2 + (n_2-1)s_2^2)/(n_1+n_2-2)}, on the same
#'     data scale as the fold change. Needs only \code{run_foldchange()}.
#'     Raw-intensity variance typically grows with abundance; consider
#'     transformed data for a variance-stabilised view.}
#'   \item{\code{"effect_size"}}{Per-feature estimate from
#'     \code{run_diff()} (multi-outcome; e.g. Cohen's d, rank-biserial,
#'     eta-squared depending on the test chosen). Supply
#'     \code{list(fc = ..., diff = ...)}. The estimate is one value per
#'     feature, so with more than two groups it is omnibus and repeats across
#'     comparisons. Estimates from different effect-size metrics are not
#'     directly comparable; a warning is given when they are mixed.}
#'   \item{\code{"vip"}}{VIP score from \code{run_pls()} (PLS-DA/OPLS-DA;
#'     sPLS-DA yields none). Supply \code{list(fc = ..., pls = ...)}. One
#'     value per feature, repeated across comparisons.}
#' }
#'
#' @return A list of class \code{plot_fcmetric} with \code{plots} (named list
#'   of \code{ggplot} objects), \code{classified} (data with the metric in
#'   \code{metric_value} and \code{.regulation}) and \code{params}.
#'
#' @examples
#' \dontrun{
#' fc   <- run_foldchange(x, metadata = meta, group = "Group")
#' diff <- run_diff(x, metadata = meta, outcome = colnames(x), group = "Group",
#'                  summary_table = TRUE)
#' pls  <- run_pls(x, metadata = meta, group = "Group", method = "plsda")
#'
#' plot_fcmetric(fc)                                      # MA plot
#' plot_fcmetric(fc, metric = "sd")
#' plot_fcmetric(list(fc = fc, diff = diff), metric = "effect_size",
#'               x_cut = c(0.2, 0.5, 0.8))
#' plot_fcmetric(list(fc = fc, pls = pls), metric = "vip", x_cut = 1)
#' }
#'
#' @author John Lennon L. Calorio
#' @seealso \code{\link{plot_volcano}}, \code{\link{run_foldchange}},
#'   \code{\link{run_diff}}, \code{\link{run_pls}}
#' @export
plot_fcmetric <- function(
    x,
    y                = NULL,
    z                = NULL,
    group            = NULL,
    features         = NULL,
    metric           = c("abundance", "effect_size", "variance", "sd", "vip"),
    log2             = TRUE,
    log10_x          = NULL,
    x_cut            = NULL,
    up               = 1.5,
    down             = 0.5,
    annotate         = NULL,
    use_ggrepel      = TRUE,
    up_color         = "red",
    down_color       = "blue",
    ns_color         = "grey60",
    h_color          = NULL,
    v_color          = NULL,
    size             = 5,
    alpha            = 0.6,
    theme            = "nature",
    plot_title       = NULL,
    plot_subtitle    = NULL,
    xlab             = NULL,
    ylab             = NULL,
    global_font_size = 25,
    title_size       = NULL,
    subtitle_size    = NULL,
    xlab_size        = NULL,
    ylab_size        = NULL,
    axis_text_size   = NULL,
    label_size       = NULL,
    legend_size      = NULL,
    legend_title     = NULL,
    legend_nrow      = NULL,
    ...
) {

  metric <- match.arg(metric)
  if (is.null(log10_x)) log10_x <- metric %in% c("abundance", "variance", "sd")

  # 1. Locate source objects (run_foldchange directly, or inside a list)
  objs  <- if (inherits(x, "run_foldchange")) list(x)
           else if (is.list(x) && !is.data.frame(x)) x else list()
  .find <- function(cls) for (o in objs) if (inherits(o, cls)) return(o)
  fc_obj   <- .find("run_foldchange")
  diff_obj <- .find("run_diff")
  pls_obj  <- .find("run_pls")

  es_metric_name <- NULL
  if (!is.null(fc_obj)) {
    d <- fc_obj$summary_table
    feats <- unique(d$feature)

    if (metric %in% c("variance", "sd") && !"sd_num" %in% names(d))
      stop("`x` lacks per-group SDs; re-run `run_foldchange()` with the current version.",
           call. = FALSE)
    if (metric == "effect_size" && is.null(diff_obj))
      stop("`metric = \"effect_size\"` needs a multi-outcome `run_diff()` result, ",
           "e.g. `x = list(fc = fc, diff = diff)`.", call. = FALSE)
    if (metric == "vip" && is.null(pls_obj))
      stop("`metric = \"vip\"` needs a `run_pls()` result, ",
           "e.g. `x = list(fc = fc, pls = pls)`.", call. = FALSE)

    d$metric_value <- switch(metric,
      abundance = rowMeans(cbind(d$mean_num, d$mean_den)),
      variance  = , sd = {
        v <- ((d$n_num - 1) * d$sd_num^2 + (d$n_den - 1) * d$sd_den^2) /
             (d$n_num + d$n_den - 2)
        if (metric == "sd") sqrt(v) else v
      },
      effect_size = {
        es <- lapply(feats, function(f) {
          e <- if (is.list(diff_obj[[f]])) diff_obj[[f]]$effect_size else NULL
          list(est = e$estimate %||% NA_real_, met = e$metric %||% NA_character_)
        })
        est <- vapply(es, `[[`, numeric(1L), "est")
        met <- vapply(es, `[[`, character(1L), "met")
        if (all(is.na(est)))
          stop("No effect sizes found in `diff` for these features. ",
               "`run_diff()` must be multi-outcome (one outcome per feature).",
               call. = FALSE)
        mets <- unique(stats::na.omit(met))
        if (length(mets) > 1L)
          warning("Effect sizes use different metrics (", paste(mets, collapse = ", "),
                  "); they are not directly comparable.", call. = FALSE)
        else es_metric_name <- mets
        est[match(d$feature, feats)]
      },
      vip = {
        v <- pls_obj$vip_scores
        if (length(v) == 0L)
          stop("`pls` has no VIP scores (sPLS-DA yields none, and a model ",
               "without a predictive component yields none).", call. = FALSE)
        out <- unname(v[d$feature])
        if (all(is.na(out)))
          stop("No feature names in `fc` match the VIP scores in `pls`.", call. = FALSE)
        out
      }
    )
    if (is.null(y))        y        <- "fold_change"
    if (is.null(z))        z        <- "metric_value"
    if (is.null(features)) features <- "feature"
    if (is.null(group) && length(unique(d$comparison)) > 1L) group <- "comparison"
    if (missing(up))   up   <- fc_obj$params$up   %||% up
    if (missing(down)) down <- fc_obj$params$down %||% down
    x <- d
  } else if (is.matrix(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else if (!is.data.frame(x)) {
    stop("`x` must be a run_foldchange result (or list containing one), ",
         "data.frame, tibble, or named matrix.", call. = FALSE)
  } else x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(x) == 0L) stop("`x` has zero rows.", call. = FALSE)

  # 2. Validate
  .chk_col <- function(nm, arg, numeric = TRUE) {
    if (!is.character(nm) || length(nm) != 1L || is.na(nm))
      stop("`", arg, "` must be a single non-NA character string.", call. = FALSE)
    if (!nm %in% colnames(x))
      stop("`", arg, "` = \"", nm, "\" is not a column in `x`.", call. = FALSE)
    if (numeric && !is.numeric(x[[nm]]))
      stop("Column `", arg, "` (\"", nm, "\") must be numeric.", call. = FALSE)
  }
  if (is.null(y) || is.null(z))
    stop("`y` and `z` must be specified when `x` is not a run_foldchange result.",
         call. = FALSE)
  .chk_col(y, "y"); .chk_col(z, "z")
  if (!is.null(group))    .chk_col(group, "group", FALSE)
  if (!is.null(features)) .chk_col(features, "features", FALSE)

  if (!is.numeric(up) || length(up) != 1L || is.na(up) || up <= 1)
    stop("`up` must be a single numeric value greater than 1.", call. = FALSE)
  if (!is.numeric(down) || length(down) != 1L || is.na(down) || down <= 0 || down >= 1)
    stop("`down` must be a single numeric value in (0, 1).", call. = FALSE)
  for (nm in c("log2", "log10_x", "use_ggrepel"))
    if (!is.logical(get(nm)) || length(get(nm)) != 1L || is.na(get(nm)))
      stop("`", nm, "` must be TRUE or FALSE.", call. = FALSE)
  for (nm in c("size", "global_font_size"))
    if (!is.numeric(get(nm)) || length(get(nm)) != 1L || get(nm) <= 0)
      stop("`", nm, "` must be a single positive numeric value.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha > 1)
    stop("`alpha` must be a single numeric value in (0, 1].", call. = FALSE)
  if (!is.null(x_cut) && (!is.numeric(x_cut) || anyNA(x_cut)))
    stop("`x_cut` must be a numeric vector without NA.", call. = FALSE)

  valid_themes <- c("nature", "minimal", "bw", "classic", "plos", "gray", "light", "dark")
  if (!is.character(theme) || length(theme) != 1L || !theme %in% valid_themes)
    stop("`theme` must be one of: ", paste(valid_themes, collapse = ", "), ".",
         call. = FALSE)
  if (!is.null(legend_nrow) &&
      (!is.numeric(legend_nrow) || length(legend_nrow) != 1L ||
       is.na(legend_nrow) || legend_nrow < 1 || legend_nrow %% 1 != 0))
    stop("`legend_nrow` must be a single positive integer.", call. = FALSE)
  if (!is.null(legend_title) && (!is.character(legend_title) || length(legend_title) != 1L))
    stop("`legend_title` must be a single character string or NULL.", call. = FALSE)

  if (any(!is.na(x[[y]]) & x[[y]] <= 0))
    stop("Column `y` (\"", y, "\") contains non-positive values. ",
         "Fold changes must be positive.", call. = FALSE)
  if (log10_x && (any(!is.na(x[[z]]) & x[[z]] <= 0) ||
                  (!is.null(x_cut) && any(x_cut <= 0)))) {
    warning("The metric (or `x_cut`) has non-positive values; using `log10_x = FALSE`.",
            call. = FALSE)
    log10_x <- FALSE
  }

  # 3. Annotation lookup
  ann_lookup <- NULL
  if (!is.null(annotate)) {
    if (is.null(features))
      stop("`features` must be specified when `annotate` is supplied.", call. = FALSE)
    if (!is.character(annotate) || length(annotate) == 0L)
      stop("`annotate` must be a non-empty character vector.", call. = FALSE)
    ids  <- unname(annotate)
    miss <- setdiff(ids, as.character(x[[features]]))
    if (length(miss))
      stop("The following `annotate` identifiers were not found in `features`: ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    ann_lookup <- stats::setNames(
      if (!is.null(names(annotate)) && any(nzchar(names(annotate)))) names(annotate) else ids,
      ids)
    if (use_ggrepel && !requireNamespace("ggrepel", quietly = TRUE))
      stop("Package 'ggrepel' is required when `use_ggrepel = TRUE`.", call. = FALSE)
  }

  # 4. Derived columns
  tf <- function(v) if (log10_x) log10(v) else v
  x$.x_plot <- tf(x[[z]])
  x$.y_plot <- if (log2) log2(x[[y]]) else x[[y]]
  up_line   <- if (log2) log2(up)   else up
  down_line <- if (log2) log2(down) else down
  ref_line  <- if (log2) 0 else 1
  x$.regulation <- ifelse(is.na(x[[y]]), "NS",
                   ifelse(x[[y]] >= up, "Up", ifelse(x[[y]] <= down, "Down", "NS")))
  if (!is.null(ann_lookup)) {
    f <- as.character(x[[features]])
    x$.ann_lbl <- ifelse(f %in% names(ann_lookup), ann_lookup[f], NA_character_)
  }

  metric_nm <- switch(metric,
    abundance   = "Mean Abundance",
    effect_size = if (length(es_metric_name)) paste0("Effect Size, ", es_metric_name)
                  else "Effect Size",
    variance    = "Variance",
    sd          = "Standard Deviation",
    vip         = "VIP Score")
  x_lbl <- if (!is.null(xlab)) xlab else
    if (log10_x) as.expression(bquote(log[10]~.(paste0("(", metric_nm, ")")))) else metric_nm
  y_lbl <- if (!is.null(ylab)) ylab else
    if (log2) expression(log[2]~"(Fold Change)") else "Fold Change"

  # 5. Theme / colours (same look as plot_volcano)
  gfs <- global_font_size
  .rsz <- function(usr, prop) if (!is.null(usr)) usr else gfs * prop
  fs_label  <- .rsz(label_size, 0.75)
  fs_legend <- .rsz(legend_size, 0.85)
  h_col <- if (!is.null(h_color)) h_color else ns_color
  v_col <- if (!is.null(v_color)) v_color else ns_color
  res_title <- if (is.null(legend_title)) "" else legend_title

  th <- switch(theme,
    nature = , plos = , classic = ggplot2::theme_classic(base_size = gfs),
    minimal = ggplot2::theme_minimal(base_size = gfs),
    bw      = ggplot2::theme_bw(base_size = gfs),
    gray    = ggplot2::theme_gray(base_size = gfs),
    light   = ggplot2::theme_light(base_size = gfs),
    dark    = ggplot2::theme_dark(base_size = gfs)
  ) + ggplot2::theme(
    plot.title        = ggplot2::element_text(size = .rsz(title_size, 1.2), face = "bold"),
    plot.subtitle     = ggplot2::element_text(size = .rsz(subtitle_size, 0.9), colour = "grey40"),
    axis.title.x      = ggplot2::element_text(size = .rsz(xlab_size, 1)),
    axis.title.y      = ggplot2::element_text(size = .rsz(ylab_size, 1)),
    axis.text         = ggplot2::element_text(size = .rsz(axis_text_size, 0.85)),
    legend.position   = "bottom",
    legend.justification = "center",
    legend.direction  = "horizontal",
    legend.text       = ggplot2::element_text(size = fs_legend),
    legend.title      = if (is.null(legend_title)) ggplot2::element_blank() else
      ggplot2::element_text(size = fs_legend, face = "bold")
  )

  guide_args <- Filter(Negate(is.null), list(
    override.aes = list(size = size * 1.5, alpha = 1),
    nrow = if (!is.null(legend_nrow)) as.integer(legend_nrow),
    title = res_title))
  col_scale <- ggplot2::scale_colour_manual(
    values = c(Up = up_color[1L], Down = down_color[1L], NS = ns_color),
    breaks = c("Up", "Down", "NS"),
    name   = res_title,
    guide  = do.call(ggplot2::guide_legend, guide_args))

  # 6. Plot builder
  .build_plot <- function(df, title_str) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[".x_plot"]], y = .data[[".y_plot"]],
                                          colour = .data[[".regulation"]])) +
      ggplot2::geom_point(size = size, alpha = alpha, na.rm = TRUE) +
      col_scale +
      ggplot2::geom_hline(yintercept = ref_line, colour = h_col, linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = c(down_line, up_line), linetype = "dashed",
                          colour = h_col, linewidth = 0.6) +
      ggplot2::labs(title = title_str, subtitle = plot_subtitle, x = x_lbl, y = y_lbl) +
      th
    if (!is.null(x_cut))
      p <- p + ggplot2::geom_vline(xintercept = tf(x_cut), linetype = "dashed",
                                   colour = v_col, linewidth = 0.6)

    if (!is.null(ann_lookup)) {
      ann_df <- df[!is.na(df$.ann_lbl), , drop = FALSE]
      if (nrow(ann_df) > 0L) {
        fn     <- if (use_ggrepel) ggrepel::geom_text_repel else ggplot2::geom_text
        extra  <- list(...)
        extra  <- extra[names(extra) %in% names(formals(fn))]
        defs   <- list(data = ann_df, mapping = ggplot2::aes(label = .data[[".ann_lbl"]]),
                       size = fs_label, na.rm = TRUE, show.legend = FALSE)
        if (use_ggrepel)
          defs <- c(defs, list(max.overlaps = Inf, min.segment.length = 0,
                               segment.color = "black", segment.size = 0.4))
        p <- p + do.call(fn, utils::modifyList(defs, extra))
      }
    }
    p
  }

  # 7. One plot per group level
  base_title <- if (!is.null(plot_title)) plot_title else
    paste0("Fold Change vs. ", sub(",.*", "", metric_nm))
  if (is.null(group)) {
    plots <- stats::setNames(list(.build_plot(x, base_title)), base_title)
  } else {
    lv    <- unique(as.character(x[[group]]))
    ttl   <- paste0(base_title, ": ", lv)
    plots <- stats::setNames(
      Map(function(l, t) .build_plot(x[as.character(x[[group]]) == l, , drop = FALSE], t),
          lv, ttl), ttl)
  }

  drop_cols <- intersect(c(".x_plot", ".y_plot", ".ann_lbl"), colnames(x))
  structure(
    list(
      plots      = plots,
      classified = x[, setdiff(colnames(x), drop_cols), drop = FALSE],
      params     = list(y = y, z = z, group = group, features = features,
                        metric = metric, log2 = log2, log10_x = log10_x,
                        x_cut = x_cut, up = up, down = down, annotate = annotate,
                        theme = theme)
    ),
    class = c("plot_fcmetric", "list")
  )
}


#' @title Print Method for \code{plot_fcmetric} Objects
#' @description Compact console summary of a \code{plot_fcmetric} result.
#' @param x An object as output of \code{plot_fcmetric}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @author John Lennon L. Calorio
#' @export
print.plot_fcmetric <- function(x, ...) {
  r <- x$classified$.regulation
  cat("-- plot_fcmetric results -----------------------------------\n")
  cat(sprintf("  Metric     : %s\n", x$params$metric))
  cat(sprintf("  Features   : %d  (Up: %d | Down: %d | NS: %d)\n",
              nrow(x$classified), sum(r == "Up"), sum(r == "Down"), sum(r == "NS")))
  cat(sprintf("  Thresholds : FC >= %.4f (up) | <= %.4f (down)\n",
              x$params$up, x$params$down))
  cat(sprintf("  Plots      : %s\n", paste(names(x$plots), collapse = ", ")))
  invisible(x)
}
