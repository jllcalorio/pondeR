#' Create Box or Violin Plots with Statistical Comparisons
#'
#' @title Create Box or Violin Plots with Statistical Comparisons
#'
#' @description
#' Generates box plots or violin plots with statistical tests and pairwise
#' comparisons. Automatically selects appropriate statistical tests based on
#' assumptions, or allows manual specification. Statistical testing is
#' performed internally via \code{\link{.run_diff_single}}, the unexported
#' single-outcome engine that powers \code{\link{run_diff}}.
#'
#' @param x A data frame containing the variables to plot.
#' @param outcome Character string specifying the column name for the numeric
#'   outcome variable. Used as the default for \code{plot_vars} if it is not
#'   specified.
#' @param group Character string specifying the column name for the grouping
#'   variable.
#' @param plot_vars Character vector specifying column name(s) for numeric
#'   variables to plot. If \code{NULL} (default), uses the value of
#'   \code{outcome}.
#' @param plot_type Character string specifying plot type: \code{"boxplot"} or
#'   \code{"violin"}. Default is \code{"boxplot"}.
#' @param test_type Character string or \code{NULL}. Statistical test strategy:
#'   \code{"auto"}, \code{"parametric"}, or \code{"nonparametric"}. See
#'   \code{\link{run_diff}} for details. If \code{NULL} (default), uses
#'   \code{"auto"}.
#' @param posthoc Character string or \code{NULL}. Retained for API
#'   compatibility; post-hoc selection is handled automatically by the
#'   internal engine. Default is \code{NULL}.
#' @param p_adjust_method Character string. P-value adjustment method for
#'   post-hoc comparisons: \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"},
#'   \code{"none"}. Default is \code{"BH"}.
#' @param hide_ns Logical. If \code{TRUE}, hide non-significant p-values from
#'   plot brackets. Default is \code{TRUE}.
#' @param show_p_numeric Logical. If \code{TRUE}, show exact p-values instead
#'   of significance symbols. Default is \code{FALSE}.
#' @param paired Logical. If \code{TRUE}, use paired tests. Default is
#'   \code{FALSE}.
#' @param test_alpha Numeric. Significance level for determining significance
#'   on the plot. Default is 0.05.
#' @param colors Character vector of colors or \code{NULL} for automatic soft
#'   colors. Default is \code{NULL}.
#' @param plot_title Character string or \code{NULL}. Custom plot title.
#'   Default is \code{NULL}.
#' @param plot_subtitle Character string or \code{NULL}. Custom subtitle. If
#'   \code{NULL}, the test name is shown. Default is \code{NULL}.
#' @param xlab Character string or \code{NULL}. X-axis label. Default is
#'   \code{NULL}.
#' @param ylab Character string or \code{NULL}. Y-axis label. Default is
#'   \code{NULL}.
#' @param linewidth Numeric. The width of the lines for the box/violin outlines 
#'   and statistical brackets. Default is 3.
#' @param show_global_p Logical. If \code{TRUE}, annotates the plot with the
#'   omnibus p-value. Default is \code{FALSE}.
#' @param global_p_position Numeric or \code{NULL}. Y-axis position for the
#'   global p-value annotation. Default is \code{NULL} (automatic).
#' @param global_p_x Numeric or \code{NULL}. X-axis position for the global
#'   p-value annotation. Default is \code{NULL} (automatic).
#' @param global_font_size Numeric or \code{NULL}. Base font size for
#'   proportional scaling of all text elements. Default is \code{30}.
#' @param title_size Numeric or \code{NULL}. Title font size. Default is
#'   \code{NULL}.
#' @param subtitle_size Numeric or \code{NULL}. Subtitle font size. Default
#'   is \code{NULL}.
#' @param xlab_size Numeric or \code{NULL}. X-axis label font size. Default
#'   is \code{NULL}.
#' @param ylab_size Numeric or \code{NULL}. Y-axis label font size. Default
#'   is \code{NULL}.
#' @param pvalue_size Numeric or \code{NULL}. Font size for p-value
#'   annotations (ggplot2 \code{size} units). Default is \code{NULL}.
#' @param axis_text_size Numeric or \code{NULL}. Axis tick label font size.
#'   Default is \code{NULL}.
#' @param format_p_numeric Character string or numeric. Controls p-value
#'   formatting when \code{show_p_numeric = TRUE}. \code{"auto"} (default)
#'   uses 3 decimal places.
#' @param group_order Character vector or \code{NULL}. Group order on x-axis.
#'   Default is \code{NULL}.
#' @param label_points Character string or \code{NULL}. \code{"all"} labels
#'   all jitter points; \code{"outliers"} labels only outliers. Requires
#'   \code{sample_id}. Default is \code{NULL}.
#' @param sample_id Character string or \code{NULL}. Column name for sample
#'   identifiers used in point labeling. Default is \code{NULL}.
#' @param show_mean Logical. If \code{TRUE}, adds a red dot (mean indicator)
#'   to each group. Default is \code{FALSE}.
#' @param show_agg_val Logical. If \code{TRUE} (default), shows the numeric
#'   mean (parametric) or median (non-parametric) value as text on the plot.
#' @param horizontal Logical. If \code{TRUE}, flips to horizontal orientation.
#'   Default is \code{FALSE}.
#' @param bracket_spacing Numeric. Multiplier for vertical spacing between
#'   comparison brackets. Default is 0.08.
#' @param subgroup Character vector or \code{NULL}. Column name(s) for
#'   subgroup analysis. Default is \code{NULL}.
#' @param ... Additional arguments passed to \code{ggboxplot} or
#'   \code{ggviolin} from \pkg{ggpubr}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{plots}{Named list of all generated \code{ggplot} objects.}
#'     \item{significant_plots}{Named list of plots for variables with
#'       significant results only.}
#'     \item{statistics}{Named list of internal \code{.run_diff_single} result
#'       objects for each variable, providing access to raw test objects,
#'       descriptive statistics, and post-hoc tables.}
#'     \item{subgroup_analysis}{Nested list of results per subgroup variable
#'       and level, or \code{NULL} if no subgroup was specified.}
#'   }
#'
#' @import ggpubr
#' @import dplyr
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' # Example 1: Box plot with mean dot and mean value annotation
#' results <- plot_diff(
#'   x            = PlantGrowth,
#'   outcome      = "weight",
#'   group        = "group",
#'   show_mean    = TRUE,
#'   show_agg_val = TRUE
#' )
#' print(results$plots$weight)
#'
#' # Example 2: Global font size scaling and omnibus p-value annotation
#' results2 <- plot_diff(
#'   x                = PlantGrowth,
#'   outcome          = "weight",
#'   group            = "group",
#'   global_font_size = 16,
#'   show_global_p    = TRUE
#' )
#' print(results2$plots$weight)
#'
#' # Example 3: Paired comparison with numeric p-values
#' sleep_data    <- sleep
#' sleep_data$ID <- rep(1:10, 2)
#' results3 <- plot_diff(
#'   x                = sleep_data,
#'   outcome          = "extra",
#'   group            = "group",
#'   paired           = TRUE,
#'   show_p_numeric   = TRUE,
#'   format_p_numeric = 4
#' )
#' print(results3$plots$extra)
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @export
plot_diff <- function(x,
                      outcome,
                      group,
                      plot_vars         = NULL,
                      plot_type         = "boxplot",
                      test_type         = NULL,
                      posthoc           = NULL,
                      p_adjust_method   = "BH",
                      hide_ns           = TRUE,
                      show_p_numeric    = FALSE,
                      paired            = FALSE,
                      test_alpha        = 0.05,
                      colors            = NULL,
                      plot_title        = NULL,
                      plot_subtitle     = NULL,
                      xlab              = NULL,
                      ylab              = NULL,
                      linewidth         = 3,
                      show_global_p     = FALSE,
                      global_p_position = NULL,
                      global_p_x        = NULL,
                      global_font_size  = 30,
                      title_size        = NULL,
                      subtitle_size     = NULL,
                      xlab_size         = NULL,
                      ylab_size         = NULL,
                      pvalue_size       = NULL,
                      axis_text_size    = NULL,
                      format_p_numeric  = "auto",
                      group_order       = NULL,
                      label_points      = NULL,
                      sample_id         = NULL,
                      show_mean         = FALSE,
                      show_agg_val      = TRUE,
                      horizontal        = FALSE,
                      bracket_spacing   = 0.08,
                      subgroup          = NULL,
                      ...) {

  # ---------------------------------------------------------------------------
  # 0. Package checks
  # ---------------------------------------------------------------------------
  .check_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' not installed. Attempting to install...", pkg))
      tryCatch(
        install.packages(pkg, repos = "https://cran.r-project.org"),
        error = function(e) stop(sprintf("Failed to install '%s': %s", pkg, e$message))
      )
      if (!requireNamespace(pkg, quietly = TRUE))
        stop(sprintf("Package '%s' could not be loaded after installation.", pkg))
    }
  }
  lapply(c("ggpubr", "dplyr", "ggplot2", "ggrepel", "rstatix"), .check_pkg)

  # ---------------------------------------------------------------------------
  # 1. Font size resolution
  # ---------------------------------------------------------------------------
  .rs <- function(usr, gfs, prop, hard) {
    if (!is.null(usr)) return(usr)
    if (!is.null(gfs)) return(gfs * prop)
    hard
  }
  fs_title    <- .rs(title_size,     global_font_size, 1.20, 14)
  fs_subtitle <- .rs(subtitle_size,  global_font_size, 0.90, 11)
  fs_xlab     <- .rs(xlab_size,      global_font_size, 1.00, 12)
  fs_ylab     <- .rs(ylab_size,      global_font_size, 1.00, 12)
  fs_axis     <- .rs(axis_text_size, global_font_size, 0.85, 10)
  fs_pvalue   <- .rs(pvalue_size,    global_font_size, 0.35, 3.5)
  fs_agg      <- .rs(NULL,           global_font_size, 0.30, 3.0)

  # ---------------------------------------------------------------------------
  # 2. p-value formatter
  # ---------------------------------------------------------------------------
  .fmt_p <- function(p_vec) {
    decimals <- if (identical(format_p_numeric, "auto")) 3L else as.integer(format_p_numeric)
    fmt_str  <- paste0("%.", decimals, "f")
    sapply(p_vec, function(pv) {
      if (is.na(pv)) return(NA_character_)
      r <- round(pv, decimals)
      if (r < 10^(-decimals))
        return(paste0("<.", paste(rep("0", decimals - 1), collapse = ""), "1"))
      if (r > (1 - 10^(-decimals)))
        return(paste0(">.", paste(rep("9", decimals), collapse = "")))
      sub("^0\\.", ".", sprintf(fmt_str, r))
    })
  }

  # ---------------------------------------------------------------------------
  # 3. Small helpers
  # ---------------------------------------------------------------------------
  .signif_stars <- function(p) {
    sapply(p, function(pv) {
      if (is.na(pv)) return("ns")
      if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns"
    })
  }

  format_agg_val <- function(vals, label_prefix = "M") {
    vapply(vals, function(val) {
      if (is.na(val)) return(NA_character_)
      fmt <- if (abs(round(val, 2)) < 0.01) gsub("E", "e", sprintf("%.2e", val))
             else sprintf("%.2f", val)
      paste0(label_prefix, "=", fmt)
    }, character(1))
  }

  default_colors <- c("#7BAFD4", "#9AC5A8", "#E8B4A0", "#C4A5D8",
                      "#F4C48D", "#A8D5BA", "#D4A5C0", "#B8C4D6",
                      "#C9B8A0", "#A5C4C4", "#D8B8C4", "#B4D4A5")
  plot_colors <- if (!is.null(colors)) colors else default_colors

  # ---------------------------------------------------------------------------
  # 4. Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x))
    stop("'x' must be a data frame.")
  if (!all(c(outcome, group) %in% names(x)))
    stop("Outcome or group variable not found in x.")

  if (!is.null(subgroup)) {
    if (!is.character(subgroup))
      stop("'subgroup' must be a character vector of column names.")
    missing_sg <- setdiff(subgroup, names(x))
    if (length(missing_sg) > 0)
      stop(paste("Subgroup column(s) not found in data:", paste(missing_sg, collapse = ", ")))
    non_cat <- subgroup[sapply(subgroup, function(s) is.numeric(x[[s]]) && !is.factor(x[[s]]))]
    if (length(non_cat) > 0)
      stop(paste("Subgroup column(s) must be categorical:", paste(non_cat, collapse = ", ")))
  }

  if (is.null(plot_vars)) plot_vars <- outcome
  if (!all(plot_vars %in% names(x)))
    stop(paste("Variable(s) to plot not found in data:",
               paste(setdiff(plot_vars, names(x)), collapse = ", ")))

  # ---------------------------------------------------------------------------
  # 5. Internal stat engine (calls the unexported helper directly)
  # ---------------------------------------------------------------------------
  .get_stats <- function(data_in, var) {
    res <- NULL
    tryCatch({
      res <- .run_diff_single(
        x                     = data_in,
        outcome               = var,
        group                 = group,
        test_type             = if (is.null(test_type)) "auto" else test_type,
        p_adjust_method       = p_adjust_method,
        paired                = paired,
        test_alpha            = test_alpha,
        calculate_effect_size = TRUE,
        perform_posthoc       = TRUE,
        group_order           = group_order,
        verbose               = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Statistical test failed for '%s': %s", var, e$message), call. = FALSE)
    })
    res
  }

  # ---------------------------------------------------------------------------
  # 6. Main loop over plot variables
  # ---------------------------------------------------------------------------
  all_plots         <- list()
  significant_plots <- list()
  statistics        <- list()
  subgroup_analysis <- list()

  for (var in plot_vars) {
    if (!is.numeric(x[[var]])) {
      warning(paste("Variable", var, "is not numeric. Skipping."))
      next
    }

    stat_results <- .get_stats(x, var)
    if (is.null(stat_results)) next

    n_groups        <- stat_results$n_groups
    main_test_p_val <- stat_results$test_result$p.value
    posthoc_data    <- stat_results$posthoc_result
    is_parametric   <- isTRUE(stat_results$parametric)

    # # Prepare plotting data from raw_data stored in the rich list
    # plot_data           <- stat_results$raw_data
    # colnames(plot_data) <- c(var, group)

    # if (!is.null(sample_id)) {
    #   if (!sample_id %in% names(x))
    #     stop(paste("Sample ID variable", sample_id, "not found in data."))
    #   plot_data[[sample_id]] <- x[[sample_id]][stats::complete.cases(x[c(var, group)])]
    # }

    # y_max      <- max(plot_data[[var]], na.rm = TRUE)
    # y_min      <- min(plot_data[[var]], na.rm = TRUE)
    # y_range    <- y_max - y_min
    # y_max_plot <- y_max + y_range * 0.15

    # # -------------------------------------------------------------------------
    # # 6a. Bracket annotation table
    # # -------------------------------------------------------------------------
    # stat_result_plot <- data.frame()

    # if (n_groups > 2 && !is.null(posthoc_data) && nrow(posthoc_data) > 0) {
    #   req_cols <- c("group1", "group2", "p.adj", "p.adj.signif")
    #   if (!all(req_cols %in% names(posthoc_data))) {
    #     warning(paste("Missing expected columns in post-hoc results for", var))
    #   } else {
    #     stat_result_plot <- posthoc_data %>%
    #       dplyr::select(group1, group2, p.adj, p.adj.signif) %>%
    #       dplyr::mutate(
    #         y.position = y_max_plot + y_range * (0.05 + bracket_spacing * (dplyr::row_number() - 1))
    #       )
    #     if (hide_ns)
    #       stat_result_plot <- dplyr::filter(stat_result_plot, p.adj < test_alpha)
    #   }
    # } else if (n_groups == 2) {
    #   p_sig <- .signif_stars(main_test_p_val)
    #   if (!hide_ns || p_sig != "ns") {
    #     grp_names <- levels(plot_data[[group]])
    #     if (is.null(grp_names)) grp_names <- unique(as.character(plot_data[[group]]))
    #     stat_result_plot <- data.frame(
    #       group1       = grp_names[1],
    #       group2       = grp_names[2],
    #       p.adj        = main_test_p_val,
    #       p.adj.signif = p_sig,
    #       y.position   = y_max_plot + y_range * 0.05
    #     )
    #   }
    # }

    # if (nrow(stat_result_plot) > 0 && show_p_numeric)
    #   stat_result_plot$p.adj.fmt <- .fmt_p(stat_result_plot$p.adj)

    # is_significant <- if (nrow(stat_result_plot) > 0)
    #   any(stat_result_plot$p.adj < test_alpha, na.rm = TRUE)
    # else
    #   main_test_p_val < test_alpha

    # # -------------------------------------------------------------------------
    # # 6b. Labels
    # # -------------------------------------------------------------------------
    # title_text    <- if (!is.null(plot_title))    plot_title    else paste(var, "by", group)
    # subtitle_text <- if (!is.null(plot_subtitle)) plot_subtitle else stat_results$test_used
    # xlab_text     <- if (!is.null(xlab)) xlab else group
    # ylab_text     <- if (!is.null(ylab)) ylab else var
    # p_label_col   <- if (show_p_numeric) "p.adj.fmt" else "p.adj.signif"

    # # -------------------------------------------------------------------------
    # # 6c. Base plot
    # # -------------------------------------------------------------------------
    # base_theme <- ggplot2::theme(
    #   legend.position = "none",
    #   plot.title      = ggplot2::element_text(face = "bold", size = fs_title),
    #   plot.subtitle   = ggplot2::element_text(size = fs_subtitle),
    #   axis.title.x    = ggplot2::element_text(size = fs_xlab),
    #   axis.title.y    = ggplot2::element_text(size = fs_ylab),
    #   axis.text       = ggplot2::element_text(size = fs_axis)
    # )

    # # Backtick-protect var and group names so ggpubr/ggplot2 handle special characters
    # # (hyphens, spaces, parentheses) in column names safely.
    # safe_var   <- paste0("`", var, "`")
    # safe_group <- paste0("`", group, "`")

    # # -------------------------------------------------------------------------
    # # 6c. Base plot
    # # -------------------------------------------------------------------------
    # p <- if (plot_type == "boxplot") {
    #   ggpubr::ggboxplot(plot_data, x = safe_group, y = safe_var,
    #                     color = safe_group, palette = plot_colors,
    #                     add = "jitter", order = group_order, ...)
    # } else {
    #   ggpubr::ggviolin(plot_data, x = safe_group, y = safe_var,
    #                   color = safe_group, fill = safe_group, palette = plot_colors,
    #                   add = "boxplot", add.params = list(fill = "white"),
    #                   order = group_order, ...)
    # }

    # p <- p +
    #   ggplot2::labs(title = title_text, subtitle = subtitle_text,
    #                 x = xlab_text, y = ylab_text) +
    #   base_theme

    # # -------------------------------------------------------------------------
    # # 6d. Mean dot
    # # -------------------------------------------------------------------------
    # if (show_mean) {
    #   p <- p + ggplot2::stat_summary(
    #     mapping = ggplot2::aes(x = .data[[group]], y = .data[[var]]),
    #     fun     = mean, geom = "point",
    #     shape   = 21, size = 3, fill = "red", color = "black"
    #   )
    # }

    # # -------------------------------------------------------------------------
    # # 6e. Aggregate value text
    # # -------------------------------------------------------------------------
    # if (show_agg_val) {
    #   # Pre-compute aggregate labels into a summary data frame.
    #   # This avoids eval(parse()) and after_stat() entirely, both of which
    #   # break on column names containing spaces, hyphens, or parentheses.
    #   agg_fun    <- if (is_parametric) mean else median
    #   agg_prefix <- if (is_parametric) "M" else "Med"
    #   agg_color  <- if (is_parametric) "red" else "blue"

    #   agg_df <- do.call(rbind, lapply(unique(plot_data[[group]]), function(g) {
    #     vals <- plot_data[[var]][plot_data[[group]] == g]
    #     agg  <- agg_fun(vals, na.rm = TRUE)
    #     data.frame(
    #       grp_col = g,
    #       agg_val = agg,
    #       agg_lbl = format_agg_val(agg, label_prefix = agg_prefix),
    #       stringsAsFactors = FALSE
    #     )
    #   }))
    #   names(agg_df)[1] <- group

    #   p <- p + ggplot2::geom_text(
    #     data    = agg_df,
    #     mapping = ggplot2::aes(
    #       x     = .data[[group]],
    #       y     = agg_val,
    #       label = agg_lbl
    #     ),
    #     vjust  = -0.7,
    #     color  = agg_color,
    #     size   = fs_agg,
    #     inherit.aes = FALSE
    #   )
    # }

    # # -------------------------------------------------------------------------
    # # 6f. Point labels
    # # -------------------------------------------------------------------------
    # if (!is.null(label_points)) {
    #   if (label_points == "all") {
    #     p <- p + ggrepel::geom_text_repel(
    #       data    = plot_data,
    #       mapping = ggplot2::aes(x = .data[[group]], y = .data[[var]],
    #                              label = .data[[sample_id]]),
    #       size = fs_agg, max.overlaps = Inf
    #     )
    #   } else if (label_points == "outliers") {
    #     outliers <- rstatix::identify_outliers(plot_data, var)
    #     if (nrow(outliers) > 0) {
    #       outlier_df <- outliers %>% dplyr::left_join(plot_data, by = c(var, group))
    #       p <- p + ggrepel::geom_text_repel(
    #         data    = outlier_df,
    #         mapping = ggplot2::aes(x = .data[[group]], y = .data[[var]],
    #                                label = .data[[sample_id]]),
    #         size = fs_agg, color = "red", max.overlaps = Inf
    #       )
    #     }
    #   }
    # }

    # Prepare plotting data from raw_data stored in the rich list
    plot_data           <- stat_results$raw_data
    
    # -------------------------------------------------------------------------
    # CRITICAL FIX: Rename columns to safe internal names to bypass ggpubr's
    # formula parser crashing on special characters (spaces, colons, hyphens).
    # -------------------------------------------------------------------------
    safe_var   <- ".outcome_var"
    safe_group <- ".group_var"
    colnames(plot_data) <- c(safe_var, safe_group)

    if (!is.null(sample_id)) {
      if (!sample_id %in% names(x))
        stop(paste("Sample ID variable", sample_id, "not found in data."))
      plot_data[[sample_id]] <- x[[sample_id]][stats::complete.cases(x[c(var, group)])]
    }

    y_max      <- max(plot_data[[safe_var]], na.rm = TRUE)
    y_min      <- min(plot_data[[safe_var]], na.rm = TRUE)
    y_range    <- y_max - y_min
    y_max_plot <- y_max + y_range * 0.15

    # -------------------------------------------------------------------------
    # 6a. Bracket annotation table
    # -------------------------------------------------------------------------
    stat_result_plot <- data.frame()

    if (n_groups > 2 && !is.null(posthoc_data) && nrow(posthoc_data) > 0) {
      req_cols <- c("group1", "group2", "p.adj", "p.adj.signif")
      if (!all(req_cols %in% names(posthoc_data))) {
        warning(paste("Missing expected columns in post-hoc results for", var))
      } else {
        stat_result_plot <- posthoc_data %>%
          dplyr::select(group1, group2, p.adj, p.adj.signif) %>%
          dplyr::mutate(
            y.position = y_max_plot + y_range * (0.05 + bracket_spacing * (dplyr::row_number() - 1))
          )
        if (hide_ns)
          stat_result_plot <- dplyr::filter(stat_result_plot, p.adj < test_alpha)
      }
    } else if (n_groups == 2) {
      p_sig <- .signif_stars(main_test_p_val)
      if (!hide_ns || p_sig != "ns") {
        grp_names <- levels(plot_data[[safe_group]])
        if (is.null(grp_names)) grp_names <- unique(as.character(plot_data[[safe_group]]))
        stat_result_plot <- data.frame(
          group1       = grp_names[1],
          group2       = grp_names[2],
          p.adj        = main_test_p_val,
          p.adj.signif = p_sig,
          y.position   = y_max_plot + y_range * 0.05
        )
      }
    }

    if (nrow(stat_result_plot) > 0 && show_p_numeric)
      stat_result_plot$p.adj.fmt <- .fmt_p(stat_result_plot$p.adj)

    is_significant <- if (nrow(stat_result_plot) > 0)
      any(stat_result_plot$p.adj < test_alpha, na.rm = TRUE)
    else
      main_test_p_val < test_alpha

    # -------------------------------------------------------------------------
    # 6b. Labels
    # -------------------------------------------------------------------------
    title_text    <- if (!is.null(plot_title))    plot_title    else paste(var, "by", group)
    subtitle_text <- if (!is.null(plot_subtitle)) plot_subtitle else stat_results$test_used
    xlab_text     <- if (!is.null(xlab)) xlab else group
    ylab_text     <- if (!is.null(ylab)) ylab else var
    p_label_col   <- if (show_p_numeric) "p.adj.fmt" else "p.adj.signif"

    # -------------------------------------------------------------------------
    # 6c. Base plot
    # -------------------------------------------------------------------------
    base_theme <- ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(face = "bold", size = fs_title),
      plot.subtitle   = ggplot2::element_text(size = fs_subtitle),
      axis.title.x    = ggplot2::element_text(size = fs_xlab),
      axis.title.y    = ggplot2::element_text(size = fs_ylab),
      axis.text       = ggplot2::element_text(size = fs_axis)
    )

    p <- if (plot_type == "boxplot") {
      ggpubr::ggboxplot(plot_data, x = safe_group, y = safe_var,
                        color = safe_group, palette = plot_colors,
                        linewidth = linewidth,
                        add = "jitter", order = group_order, ...)
    } else {
      ggpubr::ggviolin(plot_data, x = safe_group, y = safe_var,
                      color = safe_group, fill = safe_group, palette = plot_colors,
                      add = "boxplot", add.params = list(fill = "white"),
                      linewidth = linewidth,
                      add = "boxplot", add.params = list(fill = "white", size = linewidth),
                      order = group_order, ...)
    }

    p <- p +
      ggplot2::labs(title = title_text, subtitle = subtitle_text,
                    x = xlab_text, y = ylab_text) +
      base_theme

    # -------------------------------------------------------------------------
    # 6d. Mean dot
    # -------------------------------------------------------------------------
    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        mapping = ggplot2::aes(x = .data[[safe_group]], y = .data[[safe_var]]),
        # fun     = mean, geom = "point",
        fun     = mean, geom = "point", stroke = linewidth * 0.5,
        shape   = 21, size = 3, fill = "red", color = "black"
      )
    }

    # -------------------------------------------------------------------------
    # 6e. Aggregate value text
    # -------------------------------------------------------------------------
    if (show_agg_val) {
      agg_fun    <- if (is_parametric) mean else median
      agg_prefix <- if (is_parametric) "M" else "Med"
      agg_color  <- if (is_parametric) "red" else "blue"

      agg_df <- do.call(rbind, lapply(unique(plot_data[[safe_group]]), function(g) {
        vals <- plot_data[[safe_var]][plot_data[[safe_group]] == g]
        agg  <- agg_fun(vals, na.rm = TRUE)
        data.frame(
          grp_col = g,
          agg_val = agg,
          agg_lbl = format_agg_val(agg, label_prefix = agg_prefix),
          stringsAsFactors = FALSE
        )
      }))
      names(agg_df)[1] <- safe_group

      p <- p + ggplot2::geom_text(
        data    = agg_df,
        mapping = ggplot2::aes(
          x     = .data[[safe_group]],
          y     = agg_val,
          label = agg_lbl
        ),
        vjust  = -0.7,
        color  = agg_color,
        size   = fs_agg,
        inherit.aes = FALSE
      )
    }

    # -------------------------------------------------------------------------
    # 6f. Point labels
    # -------------------------------------------------------------------------
    if (!is.null(label_points)) {
      if (label_points == "all") {
        p <- p + ggrepel::geom_text_repel(
          data    = plot_data,
          mapping = ggplot2::aes(x = .data[[safe_group]], y = .data[[safe_var]],
                                 label = .data[[sample_id]]),
          size = fs_agg, max.overlaps = Inf
        )
      } else if (label_points == "outliers") {
        outliers <- rstatix::identify_outliers(plot_data, safe_var)
        if (nrow(outliers) > 0) {
          outlier_df <- outliers %>% dplyr::left_join(plot_data, by = c(safe_var, safe_group))
          p <- p + ggrepel::geom_text_repel(
            data    = outlier_df,
            mapping = ggplot2::aes(x = .data[[safe_group]], y = .data[[safe_var]],
                                   label = .data[[sample_id]]),
            size = fs_agg, color = "red", max.overlaps = Inf
          )
        }
      }
    }

    # -------------------------------------------------------------------------
    # 6g. P-value brackets
    # -------------------------------------------------------------------------
    if (nrow(stat_result_plot) > 0) {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_result_plot, label = p_label_col,
        hide.ns = FALSE, tip.length = 0.01, size = fs_pvalue,
        bracket.size = linewidth * 0.3
      )
    }

    # -------------------------------------------------------------------------
    # 6h. Global p-value annotation
    # -------------------------------------------------------------------------
    if (show_global_p) {
      global_y  <- if (!is.null(global_p_position)) global_p_position
                   else y_max_plot + y_range * 0.25
      global_x  <- if (!is.null(global_p_x)) global_p_x else n_groups * 0.95
      gp_label  <- paste0(stat_results$test_used, "\np = ", .fmt_p(main_test_p_val))
      p <- p + ggplot2::annotate("text", x = global_x, y = global_y,
                                 label = gp_label, hjust = 1, size = fs_pvalue)
    }

    if (horizontal) p <- p + ggplot2::coord_flip()

    all_plots[[var]]  <- p
    statistics[[var]] <- stat_results
    if (is_significant) significant_plots[[var]] <- p
  }

  # ---------------------------------------------------------------------------
  # 7. Subgroup analysis
  # ---------------------------------------------------------------------------
  if (!is.null(subgroup)) {
    for (sg_var in subgroup) {
      sg_levels <- if (is.factor(x[[sg_var]])) {
        levels(droplevels(as.factor(x[[sg_var]])))
      } else {
        sort(unique(na.omit(as.character(x[[sg_var]]))))
      }

      sg_var_results <- list()

      for (lvl in sg_levels) {
        x_sub <- x[as.character(x[[sg_var]]) == lvl, , drop = FALSE]

        if (nrow(x_sub) == 0) {
          warning(sprintf("Subgroup '%s' level '%s' has no observations. Skipping.",
                          sg_var, lvl))
          next
        }

        sg_all_plots         <- list()
        sg_significant_plots <- list()
        sg_statistics        <- list()

        for (var in plot_vars) {
          if (!is.numeric(x_sub[[var]])) next
          grp_sub <- droplevels(as.factor(x_sub[[group]]))
          if (nlevels(grp_sub) < 2) {
            warning(sprintf("Subgroup '%s' = '%s': '%s' has <2 group levels. Skipping.",
                            sg_var, lvl, var))
            next
          }

          sg_stat <- .get_stats(x_sub, var)
          if (is.null(sg_stat)) next

          sg_plot_result <- tryCatch({
            plot_diff(
              x                 = x_sub,
              outcome           = var,
              group             = group,
              plot_vars         = var,
              plot_type         = plot_type,
              test_type         = if (is.null(test_type)) "auto" else test_type,
              posthoc           = posthoc,
              p_adjust_method   = p_adjust_method,
              hide_ns           = hide_ns,
              show_p_numeric    = show_p_numeric,
              paired            = paired,
              test_alpha        = test_alpha,
              colors            = colors,
              plot_title        = if (!is.null(plot_title)) plot_title else
                                    sprintf("%s by %s\n[%s = %s]", var, group, sg_var, lvl),
              plot_subtitle     = plot_subtitle,
              xlab              = xlab,
              ylab              = ylab,
              linewidth         = linewidth,
              show_global_p     = show_global_p,
              global_p_position = global_p_position,
              global_p_x        = global_p_x,
              global_font_size  = global_font_size,
              title_size        = title_size,
              subtitle_size     = subtitle_size,
              xlab_size         = xlab_size,
              ylab_size         = ylab_size,
              pvalue_size       = pvalue_size,
              axis_text_size    = axis_text_size,
              format_p_numeric  = format_p_numeric,
              group_order       = group_order,
              label_points      = label_points,
              sample_id         = sample_id,
              show_mean         = show_mean,
              show_agg_val      = show_agg_val,
              horizontal        = horizontal,
              bracket_spacing   = bracket_spacing,
              subgroup          = NULL   # prevent infinite recursion
            )
          }, error = function(e) {
            warning(sprintf("plot_diff failed for subgroup '%s'='%s', var '%s': %s",
                            sg_var, lvl, var, e$message), call. = FALSE)
            NULL
          })

          sg_statistics[[var]] <- sg_stat
          if (!is.null(sg_plot_result)) {
            sg_all_plots[[var]] <- sg_plot_result$plots[[var]]
            if (!is.null(sg_plot_result$significant_plots[[var]]))
              sg_significant_plots[[var]] <- sg_plot_result$significant_plots[[var]]
          }
        }

        sg_var_results[[lvl]] <- list(
          plots             = sg_all_plots,
          significant_plots = sg_significant_plots,
          statistics        = sg_statistics
        )
      }

      subgroup_analysis[[sg_var]] <- sg_var_results
    }
  }

  list(
    plots             = all_plots,
    significant_plots = significant_plots,
    statistics        = statistics,
    subgroup_analysis = if (length(subgroup_analysis) > 0) subgroup_analysis else NULL
  )
}