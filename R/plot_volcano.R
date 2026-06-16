# =============================================================================
#  plot_volcano
# =============================================================================

#' @title Volcano Plot for Fold Change and p-value Data
#'
#' @description
#' Generates publication-quality volcano plots that visualise differential
#' features as a scatter of log2 fold change (x-axis) versus statistical
#' significance (y-axis). Designed to integrate directly with the output of
#' \code{run_foldchange()} and \code{run_diff()}, but accepts any
#' \code{data.frame} with the required columns. When \code{group} is supplied,
#' all comparisons are overlaid on a single plot with distinct colours per
#' group, saving publication space.
#'
#' @param x A \code{data.frame}, \code{tibble}, or named \code{matrix} with
#'   one row per feature (must contain the columns named by \code{y} and
#'   \code{z}), \strong{or} a named \code{list} containing the output of
#'   \code{\link{run_foldchange}()} and/or a multi-outcome
#'   \code{\link{run_diff}()} result (with \code{summary_table = TRUE}).
#'   Column names may contain special characters.
#'
#'   \strong{List input - supported combinations:}
#'   \describe{
#'     \item{\code{list(fc = <run_foldchange>, diff = <run_diff>)}}{
#'       Both objects are supplied. Fold changes are sourced from
#'       \code{run_foldchange()}'s \code{$summary_table}
#'       and p-values are sourced from
#'       \code{run_diff()}'s \code{$summary_table}, matched on feature name.
#'       When \code{run_foldchange()} produced more than one pairwise
#'       comparison, the \code{comparison} column is automatically used as
#'       \code{group}, overlaying all comparisons on a single plot with
#'       distinct colours. This is the recommended usage.}
#'     \item{\code{list(fc = <run_foldchange>)}}{
#'       Fold-change object only. Log2 fold changes are available but p-values are unavailable; all features are
#'       classified as non-significant and the significance threshold line is
#'       not meaningful. A warning is issued.}
#'     \item{\code{list(diff = <run_diff>)}}{
#'       Differential-analysis object only. Fold changes are set to 1
#'       (\eqn{\log_2 \text{FC} = 0}) for all features; the x-axis is not
#'       meaningful. A warning is issued.}
#'   }
#'
#'   When a list is supplied, the arguments \code{y}, \code{z},
#'   \code{features}, and \code{group} are set automatically and do not need
#'   to be specified by the caller (though they can be overridden).
#' @param y A single character string naming the column in \code{x} that
#'   contains \strong{raw} fold change values (numeric). The function
#'   automatically applies a log2 transformation internally. The \code{up}
#'   and \code{down} thresholds are also specified as raw fold changes
#'   (e.g., \code{up = 1.5}).
#' @param z A single character string naming the column in \code{x} that
#'   contains raw p-values (numeric, in \eqn{[0, 1]}).
#' @param group An optional single character string naming the column in
#'   \code{x} that identifies the comparison group each row belongs to (e.g.
#'   \code{"Severe_vs_Healthy"}). When supplied, all groups are overlaid on
#'   one plot and coloured by group. Default \code{NULL}.
#' @param features A single character string naming the column in \code{x}
#'   that contains the feature identifiers. Required when \code{annotate} is
#'   not \code{NULL}.
#' @param up A single numeric value greater than 1. Features with a fold
#'   change \eqn{\ge \text{up}} \emph{and} \eqn{p < \text{pval}} are
#'   coloured as up-regulated. Default \code{1.5} (i.e., a 1.5-fold increase).
#' @param down A single numeric value in \eqn{(0, 1)}. Features with a fold
#'   change \eqn{\le \text{down}} \emph{and} \eqn{p < \text{pval}} are
#'   coloured as down-regulated. For example, \code{down = 0.667} corresponds
#'   to a ~1.5-fold decrease. Default \code{0.5}.
#' @param pval A single numeric value in \eqn{(0, 1]}. Default \code{0.05}.
#' @param annotate An optional vector of feature identifiers to label on the
#'   plot. Values must be found in the column specified by \code{features}.
#'   May be a plain character vector (labels shown as-is) or a named character
#'   vector created with \code{setNames(feature_ids, display_names)}, in which
#'   case the display names are shown on the plot. A line segment always
#'   connects each label to its corresponding point. Default \code{NULL}.
#' @param annotate2 An optional numeric vector of length 2 providing fold-change
#'   thresholds for automatic annotation. Features with fold change \eqn{\le}
#'   \code{annotate2[1]} or \eqn{\ge} \code{annotate2[2]} will be labelled.
#'   The values are specified as raw fold changes (not log2). To annotate only
#'   one side, use \code{NA} for the other, e.g., \code{c(0.6, NA)}. Features
#'   selected by \code{annotate2} must \strong{also} meet the significance
#'   threshold defined by \code{pval} (i.e., \eqn{p < \text{pval}}).
#'   Requires \code{features} to be specified.
#'   Labels added via \code{annotate2} use the feature identifiers
#'   found in the \code{features} column. If \code{annotate} is also provided,
#'   both sets of features are labelled. Default \code{NULL}.
#' @param use_ggrepel Logical. If \code{TRUE} (default), uses
#'   \code{ggrepel::geom_text_repel()} for non-overlapping labels with
#'   connector segments. If \code{FALSE}, uses \code{ggplot2::geom_text()}.
#' @param up_color Character vector of colours for up-regulated points. When
#'   \code{group} is \code{NULL} only the first element is used (default
#'   \code{"red"}). When \code{group} is supplied, one colour per group level
#'   is used; if fewer colours than levels are provided they are recycled with
#'   a warning. Default \code{c("red", "purple", "yellow")}.
#' @param down_color Character vector of colours for down-regulated points.
#'   Same recycling behaviour as \code{up_color}. Default
#'   \code{c("blue", "green", "cyan")}.
#' @param ns_color Colour for non-significant points. Default \code{"grey60"}.
#' @param neglog10 Logical. If \code{TRUE} (default), the y-axis shows
#'   \eqn{-\log_{10}(p)}; if \code{FALSE}, raw p-values are plotted with an
#'   inverted y-axis.
#' @param h_color Colour of the horizontal significance threshold line.
#'   Default \code{NULL} inherits \code{ns_color}.
#' @param v_color Colour(s) of the vertical fold-change threshold lines. A
#'   single string colours both lines; a two-element vector colours the
#'   \code{down} line first and the \code{up} line second. Default \code{NULL}
#'   inherits \code{ns_color}.
#' @param size A single positive numeric value controlling the point size.
#'   Default \code{5}.
#' @param alpha A single numeric value in \eqn{(0, 1]} controlling point
#'   transparency. Overlapping points appear darker due to additive blending.
#'   Default \code{0.6}.
#' @param theme A single string specifying the \pkg{ggplot2} theme. One of
#'   \code{"nature"} (default), \code{"minimal"}, \code{"bw"},
#'   \code{"classic"}, \code{"plos"}, \code{"gray"}, \code{"light"},
#'   \code{"dark"}.
#' @param plot_title Character string or \code{NULL}. Default \code{NULL}
#'   generates an automatic title.
#' @param plot_subtitle Character string or \code{NULL}. Default \code{NULL}.
#' @param xlab Character string or \code{NULL}. Default \code{NULL} uses
#'   \code{"log2 Fold Change"}.
#' @param ylab Character string or \code{NULL}. Default \code{NULL} uses
#'   \code{"-log10(p-value)"} or \code{"p-value"}.
#' @param global_font_size A positive numeric value setting the base font size
#'   for all text elements. Default \code{25}.
#' @param title_size Override for title font size. Default \code{NULL}.
#' @param subtitle_size Override for subtitle font size. Default \code{NULL}.
#' @param xlab_size Override for x-axis label font size. Default \code{NULL}.
#' @param ylab_size Override for y-axis label font size. Default \code{NULL}.
#' @param axis_text_size Override for axis tick label font size. Default
#'   \code{NULL}.
#' @param label_size Font size for point annotation labels. Default \code{NULL}
#'   (derived from \code{global_font_size} at \eqn{0.75\times}).
#' @param legend_size Font size for legend text. Respects
#'   \code{global_font_size} (\eqn{0.85\times}) when \code{NULL}. Default
#'   \code{NULL}.
#' @param legend_title A single character string for the legend title. Default
#'   \code{NULL} (no legend title shown).
#' @param legend_nrow A single positive integer controlling the number of rows
#'   in the bottom legend. Default \code{NULL} (automatic).
#' @param ... Additional arguments forwarded to
#'   \code{ggrepel::geom_text_repel()} or \code{ggplot2::geom_text()}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{plots}}{A named list of \code{ggplot} objects.}
#'   \item{\code{classified}}{A \code{data.frame} identical to \code{x} with
#'     an additional column \code{.regulation} (\code{"Up"}, \code{"Down"},
#'     or \code{"NS"}).}
#'   \item{\code{params}}{A list of all resolved parameter values.}
#' }
#'
#' @examples
#' \dontrun{
#' ## -----------------------------------------------------------------------
#' ## Example 1 - minimal synthetic data, no group
#' ## -----------------------------------------------------------------------
#' set.seed(42)
#' n   <- 200
#' dat <- data.frame(
#'   feature     = paste0("feat_", seq_len(n)),
#'   fold_change = 2^rnorm(n, 0, 2),
#'   pvalue      = runif(n, 0, 0.5)
#' )
#'
#' res <- plot_volcano(x = dat, y = "fold_change", z = "pvalue")
#' print(res$plots[[1]])
#'
#' ## -----------------------------------------------------------------------
#' ## Example 2 - asymmetric thresholds as raw fold changes, with labels
#' ## -----------------------------------------------------------------------
#' ann <- setNames(c("feat_1", "feat_2", "feat_3"),
#'                 c("Marker A", "Marker B", "Marker C"))
#'
#' res2 <- plot_volcano(
#'   x        = dat,
#'   y        = "fold_change",
#'   z        = "pvalue",
#'   features = "feature",
#'   up       = 1.5,
#'   down     = 1/1.5,
#'   annotate = ann
#' )
#' print(res2$plots[[1]])
#'
#' ## -----------------------------------------------------------------------
#' ## Example 3 - multi-group overlay
#' ## -----------------------------------------------------------------------
#' set.seed(7)
#' dat_multi <- data.frame(
#'   feature    = rep(paste0("feat_", seq_len(50)), 3),
#'   fold_change = 2^rnorm(150, 0, 1.5),
#'   pvalue     = runif(150, 0, 0.3),
#'   comparison = rep(c("Severe_vs_Healthy", "Nonsevere_vs_Healthy",
#'                      "Severe_vs_Nonsevere"), each = 50)
#' )
#'
#' res3 <- plot_volcano(
#'   x           = dat_multi,
#'   y           = "fold_change",
#'   z           = "pvalue",
#'   group       = "comparison",
#'   features    = "feature",
#'   annotate    = setNames(c("feat_1", "feat_5"), c("Marker A", "Marker B")),
#'   alpha       = 0.5,
#'   legend_nrow = 2
#' )
#' print(res3$plots[[1]])
#'
#' ## -----------------------------------------------------------------------
#' ## Example 4 - threshold-based annotation (annotate2)
#' ## -----------------------------------------------------------------------
#' # Label everything with a 4-fold change that is also significant (p < 0.05)
#' res4 <- plot_volcano(
#'   x         = dat,
#'   y         = "fold_change",
#'   z         = "pvalue",
#'   features  = "feature",
#'   annotate2 = c(0.25, 4),
#'   up        = 1.5,
#'   down      = 0.5
#' )
#' print(res4$plots[[1]])
#' }
#'
#' @author John Lennon L. Calorio
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text
#' @importFrom ggplot2 scale_colour_manual labs theme_classic theme_minimal theme_bw theme_gray
#' @importFrom ggplot2 theme_light theme_dark theme element_text element_blank scale_y_reverse guide_legend
#' @importFrom ggrepel geom_text_repel
#' @seealso \code{\link{run_foldchange}}, \code{\link{run_diff}}, \code{\link{get_volcanodata}}, \code{\link{run_DIpreprocess}}
#' @export
plot_volcano <- function(
    x,
    y,
    z,
    group            = NULL,
    features         = NULL,
    up               = 1.5,
    down             = 0.5,
    pval             = 0.05,
    annotate         = NULL,
    annotate2        = NULL,
    use_ggrepel      = TRUE,
    up_color         = c("red",  "purple", "yellow"),
    down_color       = c("blue", "green",  "cyan"),
    ns_color         = "grey60",
    neglog10         = TRUE,
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

  # ---------------------------------------------------------------------------
  # 0.  Accept run_foldchange / run_diff list as `x`
  # ---------------------------------------------------------------------------
  # Users may supply a named list of the form:
  #   list(fc = <run_foldchange result>, diff = <run_diff summary-table result>)
  # The two objects are matched on feature name and a flat data frame is
  # reconstructed so that the rest of plot_volcano works unchanged.
  #
  # Supported combinations
  # ----------------------
  # A) run_foldchange only  -> log2 FC available; p-values are NOT present,
  #    so `z` must be left as NULL and the horizontal threshold line is
  #    suppressed (p-values set to NA / 1).  A warning is issued.
  #
  # B) run_diff (multi-outcome summary_table) only -> p-values available;
  #    fold changes are NOT present.  A warning is issued and FC columns
  #    are set to 1 (log2 FC = 0) so the plot still renders.
  #
  # C) Both  -> joined on feature name; each pairwise comparison from
  #    run_foldchange becomes one facet group (or a single plot when only
  #    one comparison exists), with matching p-values from run_diff.
  #
  # Auto-detected column names (can still be overridden via y / z):
  #   y  <- "log2fc"   (log2 fold change, created internally)
  #   z  <- "pvalue"   (raw p-value, created internally)
  #   features <- "feature"
  #   group    <- comparison label column (only when >1 comparison)
  # ---------------------------------------------------------------------------
  if (is.list(x) &&
      (inherits(x, "run_foldchange") ||
       any(vapply(x, function(el) inherits(el, c("run_foldchange", "run_diff")),
                  logical(1L))))) {

    # Delegate to get_volcanodata to handle merging and matching.
    # This ensures consistency between the plot and standalone data extraction.
    v_args <- list(filter = FALSE)
    if (!missing(up))   v_args$up   <- up
    if (!missing(down)) v_args$down <- down
    if (!missing(pval)) v_args$pval <- pval
    
    # Call get_volcanodata with the objects in x
    v_data <- do.call(get_volcanodata, c(if (is.list(x) && !inherits(x, "run_foldchange")) x else list(x), v_args))
    
    # Inherit thresholds back for vertical/horizontal plot lines
    fc_obj <- NULL
    diff_obj <- NULL
    objs <- if (is.list(x) && !inherits(x, "run_foldchange")) x else list(x)
    for (o in objs) {
      if (inherits(o, "run_foldchange")) fc_obj <- o
      if (inherits(o, "run_diff") || (is.list(o) && !is.null(o$summary_table))) diff_obj <- o
    }

    if (missing(up))   up   <- fc_obj$params$up   %||% 1.5
    if (missing(down)) down <- fc_obj$params$down %||% 0.5
    if (missing(pval)) {
      pval <- diff_obj$params$p_value %||% diff_obj$params$pval %||% 
              diff_obj$params$adj_pval %||% 0.05
    }

    # Resolve column mapping
    if (missing(y) || is.null(y)) {
      # Favor raw fold_change column from run_foldchange output
      y <- if ("fold_change" %in% names(v_data)) "fold_change" else "log2_fc"
    }
    if (missing(z) || is.null(z)) {
      z <- "p_value_used"
    }
    if (missing(features) || is.null(features)) {
      features <- "feature"
    }
    if ((missing(group) || is.null(group)) && "comparison" %in% names(v_data)) {
      if (length(unique(v_data$comparison)) > 1L) {
        group <- "comparison"
      }
    }

    x <- v_data
  }   # end preprocessing block

  # ---------------------------------------------------------------------------
  # 1.  Coerce x
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x))
    stop("`x` must be a data.frame, tibble, or named matrix.", call. = FALSE)
  if (is.matrix(x)) {
    if (is.null(colnames(x)))
      stop("`x` is a matrix but has no column names.", call. = FALSE)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (nrow(x) == 0L) stop("`x` has zero rows.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2.  Validate y, z
  # ---------------------------------------------------------------------------
  .chk_col <- function(nm, arg) {
    if (!is.character(nm) || length(nm) != 1L || is.na(nm))
      stop("`", arg, "` must be a single non-NA character string.", call. = FALSE)
    if (!nm %in% colnames(x))
      stop("`", arg, "` = \"", nm, "\" is not a column in `x`.", call. = FALSE)
    if (!is.numeric(x[[nm]]))
      stop("Column `", arg, "` (\"", nm, "\") must be numeric.", call. = FALSE)
  }
  .chk_col(y, "y")
  .chk_col(z, "z")

  # ---------------------------------------------------------------------------
  # 3.  Validate group
  # ---------------------------------------------------------------------------
  grp_levels <- NULL
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L)
      stop("`group` must be a single character string.", call. = FALSE)
    if (!group %in% colnames(x))
      stop("`group` = \"", group, "\" is not a column in `x`.", call. = FALSE)
    grp_levels <- sort(unique(as.character(x[[group]])))
    if (length(grp_levels) < 2L)
      stop("`group` must have at least 2 distinct levels.", call. = FALSE)
  }
  n_grp <- if (is.null(grp_levels)) 1L else length(grp_levels)

  # ---------------------------------------------------------------------------
  # 4.  Validate / recycle up_color and down_color
  # ---------------------------------------------------------------------------
  if (!is.character(up_color)   || length(up_color)   == 0L)
    stop("`up_color` must be a non-empty character vector.", call. = FALSE)
  if (!is.character(down_color) || length(down_color) == 0L)
    stop("`down_color` must be a non-empty character vector.", call. = FALSE)

  .recycle_colors <- function(cols, n, param_nm) {
    if (length(cols) < n) {
      warning("`", param_nm, "` has ", length(cols), " colour(s) but `group` has ",
              n, " level(s). Colours will be recycled.", call. = FALSE)
      cols <- rep_len(cols, n)
    }
    cols[seq_len(n)]
  }

  up_cols   <- .recycle_colors(up_color,   n_grp, "up_color")
  down_cols <- .recycle_colors(down_color, n_grp, "down_color")

  # ---------------------------------------------------------------------------
  # 5.  Validate thresholds / numerics
  # ---------------------------------------------------------------------------
  # .chk_num1 <- function(v, nm) {
  #   if (!is.numeric(v) || length(v) != 1L || is.na(v))
  #     stop("`", nm, "` must be a single numeric value.", call. = FALSE)
  # }
  # .chk_num1(up,   "up")
  # .chk_num1(down, "down")
  if (!is.numeric(up) || length(up) != 1L || is.na(up) || up <= 1)
    stop("`up` must be a single numeric value greater than 1 (e.g., 1.5 for a 1.5-fold increase).",
         call. = FALSE)
  if (!is.numeric(down) || length(down) != 1L || is.na(down) || down <= 0 || down >= 1)
    stop("`down` must be a single numeric value in (0, 1) (e.g., 0.667 for a ~1.5-fold decrease).",
         call. = FALSE)

  if (!is.numeric(pval) || length(pval) != 1L || is.na(pval) ||
      pval <= 0 || pval > 1)
    stop("`pval` must be a single numeric value in (0, 1].", call. = FALSE)

  .chk_pos1 <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0)
      stop("`", nm, "` must be a single positive numeric value.", call. = FALSE)
  }
  .chk_pos1(size,             "size")
  .chk_pos1(global_font_size, "global_font_size")

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha > 1)
    stop("`alpha` must be a single numeric value in (0, 1].", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 5b.  Convert raw fold change thresholds to log2 scale for internal use
  # ---------------------------------------------------------------------------
  up_log2   <- log2(up)
  down_log2 <- log2(down)

  # ---------------------------------------------------------------------------
  # 6.  Validate logicals
  # ---------------------------------------------------------------------------
  .chk_bool <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop("`", nm, "` must be TRUE or FALSE.", call. = FALSE)
  }
  .chk_bool(use_ggrepel, "use_ggrepel")
  .chk_bool(neglog10,    "neglog10")

  # ---------------------------------------------------------------------------
  # 7.  Validate legend_nrow
  # ---------------------------------------------------------------------------
  if (!is.null(legend_nrow)) {
    if (!is.numeric(legend_nrow) || length(legend_nrow) != 1L ||
        is.na(legend_nrow) || legend_nrow < 1L || legend_nrow %% 1 != 0)
      stop("`legend_nrow` must be a single positive integer.", call. = FALSE)
    legend_nrow <- as.integer(legend_nrow)
  }

  # ---------------------------------------------------------------------------
  # 8.  Validate legend_title
  # ---------------------------------------------------------------------------
  if (!is.null(legend_title)) {
    if (!is.character(legend_title) || length(legend_title) != 1L)
      stop("`legend_title` must be a single character string or NULL.", call. = FALSE)
  }
  # Empty string sentinel used in scale so the title space is truly removed
  resolved_legend_title <- if (is.null(legend_title)) "" else legend_title

  # ---------------------------------------------------------------------------
  # 9.  Validate features and annotate
  # ---------------------------------------------------------------------------
  # 9a. Handle annotate2 (threshold-based selection)
  if (!is.null(annotate2)) {
    if (is.null(features)) {
      stop("`features` must be specified (a column name of `x`) when `annotate2` is supplied.",
           call. = FALSE)
    }
    if (!is.numeric(annotate2) || length(annotate2) != 2L) {
      stop("`annotate2` must be a numeric vector of length 2 (e.g. c(0.5, 2.0) or c(0.6, NA)).",
           call. = FALSE)
    }

    # Constructive logic checks
    lower_val <- annotate2[1]
    upper_val <- annotate2[2]

    if (!is.na(lower_val) && !is.na(upper_val) && lower_val >= upper_val) {
      warning("`annotate2[1]` (down) is greater than or equal to `annotate2[2]` (up). ",
              "This might lead to unintended results.", call. = FALSE)
    }

    # Identify feature IDs based on FC thresholds
    raw_fc_vals <- x[[y]]
    p_vals      <- x[[z]]
    extra_ids <- character(0)

    if (!is.na(lower_val)) {
      idx_low <- !is.na(raw_fc_vals) & raw_fc_vals <= lower_val & !is.na(p_vals) & p_vals < pval
      extra_ids <- c(extra_ids, as.character(x[[features]][idx_low]))
    }
    if (!is.na(upper_val)) {
      idx_up <- !is.na(raw_fc_vals) & raw_fc_vals >= upper_val & !is.na(p_vals) & p_vals < pval
      extra_ids <- c(extra_ids, as.character(x[[features]][idx_up]))
    }

    if (length(extra_ids) > 0L) {
      new_unique_ids <- unique(extra_ids)
      if (is.null(annotate)) {
        annotate <- new_unique_ids
      } else {
        existing_ids <- unname(annotate)
        ids_to_add <- setdiff(new_unique_ids, existing_ids)
        if (length(ids_to_add) > 0L) {
          if (!is.null(names(annotate))) {
            annotate <- c(annotate, stats::setNames(ids_to_add, ids_to_add))
          } else {
            annotate <- c(annotate, ids_to_add)
          }
        }
      }
    }
  }

  # 9b. Primary annotation processing
  ann_lookup <- NULL
  ann_ids    <- NULL

  if (!is.null(annotate)) {
    if (is.null(features))
      stop("`features` must be specified (a column name of `x`) when `annotate` is supplied.",
           call. = FALSE)
    if (!is.character(features) || length(features) != 1L)
      stop("`features` must be a single character string.", call. = FALSE)
    if (!features %in% colnames(x))
      stop("`features` = \"", features, "\" is not a column in `x`.", call. = FALSE)

    feat_vals <- as.character(x[[features]])

    if (!is.character(annotate) || length(annotate) == 0L)
      stop("`annotate` must be a non-empty character vector.", call. = FALSE)

    ann_ids     <- unname(annotate)
    missing_ids <- setdiff(ann_ids, feat_vals)
    if (length(missing_ids) > 0L)
      stop(
        "The following `annotate` identifiers were not found in the `features` column (\"",
        features, "\"): ", paste(missing_ids, collapse = ", "), ".", call. = FALSE
      )

    # Named vector -> display names; plain vector -> IDs as labels
    if (!is.null(names(annotate)) && any(nzchar(names(annotate)))) {
      ann_lookup <- stats::setNames(names(annotate), ann_ids)
    } else {
      ann_lookup <- stats::setNames(ann_ids, ann_ids)
    }

    if (use_ggrepel && !requireNamespace("ggrepel", quietly = TRUE))
      stop("Package 'ggrepel' is required when `use_ggrepel = TRUE`. ",
           "Install with: install.packages(\"ggrepel\")", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 10.  Validate theme and line colours
  # ---------------------------------------------------------------------------
  valid_themes <- c("nature", "minimal", "bw", "classic", "plos",
                    "gray",   "light",   "dark")
  if (!is.character(theme) || length(theme) != 1L || !theme %in% valid_themes)
    stop("`theme` must be one of: ", paste(valid_themes, collapse = ", "), ".",
         call. = FALSE)

  .chk_color1 <- function(v, nm, n_allowed = NULL) {
    if (is.null(v)) return(invisible(NULL))
    if (!is.character(v)) stop("`", nm, "` must be a character vector.", call. = FALSE)
    if (!is.null(n_allowed) && !length(v) %in% n_allowed)
      stop("`", nm, "` must have length ", paste(n_allowed, collapse = " or "), ".",
           call. = FALSE)
  }
  .chk_color1(ns_color, "ns_color", 1L)
  .chk_color1(h_color,  "h_color",  1L)
  .chk_color1(v_color,  "v_color",  1L:2L)

  h_col   <- if (!is.null(h_color)) h_color else ns_color
  v_col_d <- if (!is.null(v_color)) v_color[1L]              else ns_color
  v_col_u <- if (!is.null(v_color)) utils::tail(v_color, 1L) else ns_color

  # ---------------------------------------------------------------------------
  # 11.  Resolve font sizes
  # ---------------------------------------------------------------------------
  gfs <- global_font_size
  .rsz <- function(usr, prop) if (!is.null(usr)) usr else gfs * prop

  fs_title    <- .rsz(title_size,     1.20)
  fs_subtitle <- .rsz(subtitle_size,  0.90)
  fs_xlab     <- .rsz(xlab_size,      1.00)
  fs_ylab     <- .rsz(ylab_size,      1.00)
  fs_axis     <- .rsz(axis_text_size, 0.85)
  fs_legend   <- .rsz(legend_size,    0.85)
  fs_label    <- .rsz(label_size,     0.75)

  # ---------------------------------------------------------------------------
  # 12.  Handle zero p-values
  # ---------------------------------------------------------------------------
  pvals <- x[[z]]
  zero_p <- !is.na(pvals) & pvals == 0
  if (any(zero_p)) {
    warning(sum(zero_p),
            " p-value(s) equal to 0 replaced with .Machine$double.eps.",
            call. = FALSE)
    pvals[zero_p] <- .Machine$double.eps
    x[[z]] <- pvals
  }

  # ---------------------------------------------------------------------------
  # 13.  Log2-transform x[[y]] for plotting and classification
  # ---------------------------------------------------------------------------
  raw_fc <- x[[y]]

  if (any(!is.na(raw_fc) & raw_fc <= 0))
    stop("Column `y` (\"", y, "\") contains non-positive values. ",
         "Fold changes must be positive (e.g., 0.5 for a 2-fold decrease).",
         call. = FALSE)

  x$.x_plot <- log2(raw_fc)

  if (neglog10) {
    x$.y_plot <- -base::log10(x[[z]])
    h_line_y  <- -base::log10(pval)
    y_lbl     <- if (!is.null(ylab)) ylab else expression(-log[10](p))
  } else {
    x$.y_plot <- x[[z]]
    h_line_y  <- pval
    y_lbl     <- if (!is.null(ylab)) ylab else "p-value"
  }
  x_lbl <- if (!is.null(xlab)) xlab else expression(log[2]~"Fold Change")

  # ---------------------------------------------------------------------------
  # 14.  Classify features
  # ---------------------------------------------------------------------------
  x$.regulation <- ifelse(
    !is.na(x$.x_plot) & !is.na(x[[z]]) & x$.x_plot >= up_log2   & x[[z]] < pval, "Up",
    ifelse(
      !is.na(x$.x_plot) & !is.na(x[[z]]) & x$.x_plot <= down_log2 & x[[z]] < pval, "Down",
      "NS"
    )
  )

  # ---------------------------------------------------------------------------
  # 15.  Annotation label column
  # ---------------------------------------------------------------------------
  if (!is.null(annotate)) {
    feat_vals  <- as.character(x[[features]])
    x$.ann_lbl <- ifelse(
      feat_vals %in% ann_ids,
      ann_lookup[feat_vals],
      NA_character_
    )
  }

  # ---------------------------------------------------------------------------
  # 16.  Assign per-row colour key
  # ---------------------------------------------------------------------------
  if (is.null(grp_levels)) {
    # No group: key is simply "Up", "Down", "NS"
    x$.colour_key <- x$.regulation
  } else {
    grp_char      <- as.character(x[[group]])
    x$.colour_key <- ifelse(
      x$.regulation == "NS",
      "NS",
      paste0(grp_char, "_", x$.regulation)
    )
  }

  # ---------------------------------------------------------------------------
  # 17.  Theme builder
  # ---------------------------------------------------------------------------
  .build_theme <- function(theme_name) {
    base <- switch(
      theme_name,
      nature  = , plos = ggplot2::theme_classic(base_size = gfs),
      minimal = ggplot2::theme_minimal(base_size = gfs),
      bw      = ggplot2::theme_bw(base_size = gfs),
      classic = ggplot2::theme_classic(base_size = gfs),
      gray    = ggplot2::theme_gray(base_size = gfs),
      light   = ggplot2::theme_light(base_size = gfs),
      dark    = ggplot2::theme_dark(base_size = gfs)
    )
    # Legend always at bottom-centre, horizontal layout
    base + ggplot2::theme(
      plot.title       = ggplot2::element_text(size = fs_title,    face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = fs_subtitle, colour = "grey40"),
      axis.title.x     = ggplot2::element_text(size = fs_xlab),
      axis.title.y     = ggplot2::element_text(size = fs_ylab),
      axis.text        = ggplot2::element_text(size = fs_axis),
      legend.position  = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text      = ggplot2::element_text(size = fs_legend),
      # Suppress legend title spacing when no title requested
      legend.title     = if (is.null(legend_title)) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_text(size = fs_legend, face = "bold")
      }
    )
  }

  # ---------------------------------------------------------------------------
  # 18.  Colour scale builder
  # ---------------------------------------------------------------------------
  .build_colour_scale <- function() {
    guide_args <- list(
      override.aes = list(size = size * 1.5, alpha = 1),
      nrow         = legend_nrow,   # NULL = auto
      title        = resolved_legend_title
    )
    guide_args <- Filter(Negate(is.null), guide_args)
    # Always pass title even if "" so scale name is set
    guide_args$title <- resolved_legend_title

    if (is.null(grp_levels)) {
      # Single comparison - 3 levels
      vals <- c(Up = up_cols[1L], Down = down_cols[1L], NS = ns_color)
      ggplot2::scale_colour_manual(
        values = vals,
        breaks = c("Up", "Down", "NS"),
        labels = c("Up", "Down", "NS"),
        name   = resolved_legend_title,
        guide  = do.call(ggplot2::guide_legend, guide_args)
      )
    } else {
      # Multi-group
      up_vals   <- stats::setNames(up_cols,   paste0(grp_levels, "_Up"))
      down_vals <- stats::setNames(down_cols, paste0(grp_levels, "_Down"))
      ns_val    <- c(NS = ns_color)

      all_vals   <- c(up_vals, down_vals, ns_val)
      up_labels  <- paste0(grp_levels, " \u2191")   # \u2191 for up
      dn_labels  <- paste0(grp_levels, " \u2193")   # \u2193 for down
      all_labels <- c(up_labels, dn_labels, "NS")

      ggplot2::scale_colour_manual(
        values = all_vals,
        breaks = names(all_vals),
        labels = all_labels,
        name   = resolved_legend_title,
        guide  = do.call(ggplot2::guide_legend, guide_args)
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 19.  Plot builder
  # ---------------------------------------------------------------------------
  .build_plot <- function(df, title_str) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x      = .data[[".x_plot"]],
        y      = .data[[".y_plot"]],
        colour = .data[[".colour_key"]]
      )
    ) +
      ggplot2::geom_point(size = size, alpha = alpha, na.rm = TRUE) +
      .build_colour_scale() +
      ggplot2::geom_hline(yintercept = h_line_y, linetype = "dashed",
                          colour = h_col,   linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = up_log2,        linetype = "dashed",
                          colour = v_col_u, linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = down_log2,      linetype = "dashed",
                          colour = v_col_d, linewidth = 0.6) +
      ggplot2::labs(
        title    = title_str,
        subtitle = plot_subtitle,
        x        = x_lbl,
        y        = y_lbl
      ) +
      .build_theme(theme)

    if (!neglog10) p <- p + ggplot2::scale_y_reverse()

    # Annotation layer - ggrepel draws segment lines to dots by default
    if (!is.null(annotate)) {
      ann_df <- df[!is.na(df$.ann_lbl), , drop = FALSE]
      if (nrow(ann_df) > 0L) {
        extra      <- list(...)
        repel_fmls <- if (use_ggrepel) names(formals(ggrepel::geom_text_repel))
                      else             names(formals(ggplot2::geom_text))
        extra_filt <- extra[names(extra) %in% repel_fmls]

        if (use_ggrepel) {
          # segment.color = "black" ensures the connector line is always drawn;
          # min.segment.length = 0 forces a segment even for nearby points.
          repel_defaults <- list(
            data             = ann_df,
            mapping          = ggplot2::aes(label = .data[[".ann_lbl"]]),
            size             = fs_label,
            na.rm            = TRUE,
            max.overlaps     = Inf,
            show.legend      = FALSE,
            min.segment.length = 0,
            segment.color    = "black",
            segment.size     = 0.4
          )
          # User-supplied extras override defaults where keys clash
          repel_args <- utils::modifyList(repel_defaults, extra_filt)
          p <- p + do.call(ggrepel::geom_text_repel, repel_args)
        } else {
          txt_defaults <- list(
            data        = ann_df,
            mapping     = ggplot2::aes(label = .data[[".ann_lbl"]]),
            size        = fs_label,
            na.rm       = TRUE,
            show.legend = FALSE
          )
          txt_args <- utils::modifyList(txt_defaults, extra_filt)
          p <- p + do.call(ggplot2::geom_text, txt_args)
        }
      }
    }

    p
  }

  # ---------------------------------------------------------------------------
  # 20.  Produce plot
  # ---------------------------------------------------------------------------
  auto_title <- if (!is.null(plot_title)) plot_title else "Volcano Plot"
  plots_out  <- stats::setNames(list(.build_plot(x, auto_title)), auto_title)

  # ---------------------------------------------------------------------------
  # 21.  Clean up internal columns before returning classified data
  # ---------------------------------------------------------------------------
  drop_cols <- intersect(c(".x_plot", ".y_plot", ".ann_lbl", ".colour_key"), colnames(x))
  x_out     <- x[, setdiff(colnames(x), drop_cols), drop = FALSE]

  structure(
    list(
      plots      = plots_out,
      classified = x_out,
      params     = list(
        y            = y,
        z            = z,
        group        = group,
        features     = features,
        up           = up,
        down         = down,
        pval         = pval,
        neglog10     = neglog10,
        annotate     = annotate,
        annotate2    = annotate2,
        use_ggrepel  = use_ggrepel,
        up_color     = up_cols,
        down_color   = down_cols,
        ns_color     = ns_color,
        size         = size,
        alpha        = alpha,
        theme        = theme,
        legend_size  = fs_legend,
        legend_title = legend_title,
        legend_nrow  = legend_nrow
      )
    )
  )
}


# =============================================================================
#  S3 print method
# =============================================================================

#' @title Print Method for \code{plot_volcano} Objects
#' @description Compact console summary of a \code{plot_volcano} result.
#' @param x An object as ouput of \code{"plot_volcano"}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @author John Lennon L. Calorio
#' @export
print.plot_volcano <- function(x, ...) {
  clf  <- x$classified
  up_n <- sum(clf$.regulation == "Up",   na.rm = TRUE)
  dn_n <- sum(clf$.regulation == "Down", na.rm = TRUE)
  ns_n <- sum(clf$.regulation == "NS",   na.rm = TRUE)
  cat("-- plot_volcano results ------------------------------------\n")
  cat(sprintf("  Features   : %d  (Up: %d | Down: %d | NS: %d)\n",
              nrow(clf), up_n, dn_n, ns_n))
  cat(sprintf("  Thresholds : FC >= %.4f (up) | <= %.4f (down)  |  p < %.4f\n",
              x$params$up, x$params$down, x$params$pval))
  if (!is.null(x$params$group))
    cat(sprintf("  Group col  : %s\n", x$params$group))
  cat(sprintf("  Plots      : %s\n", paste(names(x$plots), collapse = ", ")))
  invisible(x)
}