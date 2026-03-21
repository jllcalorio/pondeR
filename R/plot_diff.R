#' Create Box or Violin Plots with Statistical Comparisons
#'
#' This function generates box plots or violin plots with statistical tests and
#' pairwise comparisons. It automatically selects appropriate statistical tests
#' based on assumptions, or allows manual specification. Statistical testing is
#' performed using the companion function \code{\link{run_diff}}.
#'
#' @param x A data frame containing the variables to plot.
#' @param outcome Character string specifying the column name for the numeric outcome variable.
#'   Used as the default for \code{plot_vars} if it is not specified.
#' @param group Character string specifying the column name for the grouping variable.
#' @param plot_vars Character vector specifying column name(s) for numeric variables to plot.
#'   If \code{NULL} (default), uses the value of \code{outcome}.
#' @param plot_type Character string specifying plot type: \code{"boxplot"} or \code{"violin"}.
#'   Default is \code{"boxplot"}.
#' @param test_type Character string or NULL. Statistical test strategy: \code{"auto"},
#'   \code{"parametric"}, or \code{"nonparametric"}. See \code{\link{run_diff}} for details.
#'   If NULL (default), uses \code{"auto"}.
#' @param posthoc Character string or NULL. Post-hoc test for multiple groups:
#'   \code{"auto"}, \code{"none"}, \code{"tukey"}, \code{"games-howell"}, \code{"dunn"},
#'   \code{"t_test"}, or \code{"wilcox_test"}. If NULL (default), uses \code{"auto"}.
#' @param p_adjust_method Character string. P-value adjustment method for multiple comparisons.
#'   Options: \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}. Default is \code{"BH"}.
#' @param hide_ns Logical. If TRUE, hide non-significant p-values from plot brackets.
#'   Default is TRUE.
#' @param show_p_numeric Logical. If TRUE, show exact p-values as numbers instead of
#'   significance symbols (\code{*}, \code{**}, \code{***}). Default is FALSE.
#' @param paired Logical. If TRUE, use paired tests. Default is FALSE.
#' @param test_alpha Numeric. Significance level for determining significance on the plot.
#'   Default is 0.05.
#' @param colors Character vector of colors (names or HEX codes) or NULL for automatic
#'   soft colors. Default is NULL.
#' @param plot_title Character string or NULL. Custom title for the plot. If NULL, generates
#'   an automatic title. Default is NULL.
#' @param plot_subtitle Character string or NULL. Custom subtitle for the plot. If NULL,
#'   the name of the statistical test used is shown automatically. Default is NULL.
#' @param xlab Character string or NULL. X-axis label. If NULL, uses the group variable name.
#'   Default is NULL.
#' @param ylab Character string or NULL. Y-axis label. If NULL, uses the plot variable name.
#'   Default is NULL.
#' @param show_global_p Logical. If TRUE, annotates the plot with the global (omnibus)
#'   p-value from the main statistical test (e.g., Kruskal-Wallis, ANOVA). Shown for all
#'   group counts when TRUE. Default is FALSE.
#' @param global_p_position Numeric or NULL. Y-axis position (data coordinates) for the
#'   global p-value annotation. If NULL (default), positioned automatically near the top
#'   of the plot area.
#' @param global_p_x Numeric or NULL. X-axis position for the global p-value annotation.
#'   If NULL (default), positioned at 95% of the x-axis range (right side).
#' @param global_font_size Numeric or NULL. Base font size used to derive all text sizes
#'   proportionally. When provided, \code{title_size}, \code{subtitle_size},
#'   \code{xlab_size}, \code{ylab_size}, \code{pvalue_size}, \code{axis_text_size}, and
#'   annotation text sizes (mean/median values, global p) are all scaled relative to this
#'   value. Individual size parameters override \code{global_font_size} only for their
#'   specific element when set to a non-NULL numeric. Default is NULL (each parameter uses
#'   its own explicit default).
#' @param title_size Numeric or NULL. Font size for the plot title. If NULL and
#'   \code{global_font_size} is set, derived proportionally (1.2x). Default is NULL.
#' @param subtitle_size Numeric or NULL. Font size for the plot subtitle. If NULL and
#'   \code{global_font_size} is set, derived proportionally (0.9x). Default is NULL.
#' @param xlab_size Numeric or NULL. Font size for the x-axis label. If NULL and
#'   \code{global_font_size} is set, derived proportionally (1.0x). Default is NULL.
#' @param ylab_size Numeric or NULL. Font size for the y-axis label. If NULL and
#'   \code{global_font_size} is set, derived proportionally (1.0x). Default is NULL.
#' @param pvalue_size Numeric or NULL. Font size (in ggplot2 \code{size} units, i.e. mm)
#'   for p-value annotations on brackets and the global p-value. If NULL and
#'   \code{global_font_size} is set, derived proportionally (0.35x). Default is NULL.
#' @param axis_text_size Numeric or NULL. Font size for axis tick labels (both x and y axes).
#'   If NULL and \code{global_font_size} is set, derived proportionally (0.85x).
#'   Default is NULL.
#' @param format_p_numeric Character string or numeric. Controls how p-values are
#'   formatted when \code{show_p_numeric = TRUE}. \code{"auto"} (default) rounds to 3
#'   decimal places, showing \code{"<.001"} if below that threshold and \code{">.999"} if
#'   above. A numeric value specifies the exact number of decimal places to display.
#' @param group_order Character vector or NULL. Specifies the order of groups on the x-axis.
#'   If NULL (default), uses factor levels or alphabetical order.
#' @param label_points Character string or NULL. Point labeling: \code{"all"} (label all
#'   jitter points), \code{"outliers"} (label only outliers), or NULL (no labels).
#'   Requires \code{sample_id}. Default is NULL.
#' @param sample_id Character string or NULL. Column name containing sample identifiers
#'   for point labeling. Required when \code{label_points} is not NULL.
#' @param show_mean Logical. If TRUE, adds a red dot (mean indicator) to each group.
#'   Default is FALSE.
#' @param show_agg_val Logical. If TRUE (default), shows the numeric mean (parametric) or
#'   median (non-parametric) value as text on the plot.
#' @param horizontal Logical. If TRUE, flips the plot to horizontal orientation.
#'   Default is FALSE.
#' @param bracket_spacing Numeric. Multiplier controlling the vertical spacing between
#'   comparison brackets. Default is 0.08.
#' @param subgroup Character vector or NULL. Column name(s) in \code{x} for subgroup
#'   analysis. Each must be categorical (factor or character). When specified, the full
#'   analysis is repeated for each level of each subgroup column. Default is NULL.
#' @param ... Additional arguments passed to \code{ggboxplot} or \code{ggviolin} from
#'   \pkg{ggpubr}.
#'
#' @return A list containing:
#'   \item{plots}{List of all generated plots.}
#'   \item{significant_plots}{List of plots with significant results only.}
#'   \item{statistics}{List of \code{\link{run_diff}} results for each variable.}
#'   \item{subgroup_analysis}{Nested list of results per subgroup variable and level,
#'     or NULL if no subgroup was specified.}
#'
#' @import ggpubr
#' @import dplyr
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' # Example 1: Show mean dot AND mean value
#' results_anova <- plot_diff(
#'   x            = PlantGrowth,
#'   outcome      = "weight",
#'   group        = "group",
#'   plot_type    = "boxplot",
#'   show_mean    = TRUE,
#'   show_agg_val = TRUE
#' )
#' print(results_anova$plots$weight)
#'
#' # Example 2: Global font size controls all text proportionally
#' results_anova2 <- plot_diff(
#'   x               = PlantGrowth,
#'   outcome         = "weight",
#'   group           = "group",
#'   global_font_size = 16,
#'   show_global_p   = TRUE  # show the Kruskal-Wallis / ANOVA p-value on the plot
#' )
#' print(results_anova2$plots$weight)
#'
#' # Example 3: Paired t-test with custom subtitle, formatted p-values
#' sleep_data     <- sleep
#' sleep_data$ID  <- rep(1:10, 2)
#' results_paired <- plot_diff(
#'   x              = sleep_data,
#'   outcome        = "extra",
#'   group          = "group",
#'   paired         = TRUE,
#'   plot_subtitle  = "Paired comparison: Drug 1 vs Drug 2",
#'   show_p_numeric = TRUE,
#'   format_p_numeric = 4  # 4 decimal places
#' )
#' print(results_paired$plots$extra)
#' }
#'
#' @export
plot_diff <- function(x,
                      outcome,
                      group,
                      plot_vars        = NULL,
                      plot_type        = "boxplot",
                      test_type        = NULL,
                      posthoc          = NULL,
                      p_adjust_method  = "BH",
                      hide_ns          = TRUE,
                      show_p_numeric   = FALSE,
                      paired           = FALSE,
                      test_alpha       = 0.05,
                      colors           = NULL,
                      plot_title       = NULL,
                      plot_subtitle    = NULL,
                      xlab             = NULL,
                      ylab             = NULL,
                      show_global_p    = FALSE,
                      global_p_position = NULL,
                      global_p_x       = NULL,
                      global_font_size = NULL,
                      title_size       = NULL,
                      subtitle_size    = NULL,
                      xlab_size        = NULL,
                      ylab_size        = NULL,
                      pvalue_size      = NULL,
                      axis_text_size   = NULL,
                      format_p_numeric = "auto",
                      group_order      = NULL,
                      label_points     = NULL,
                      sample_id        = NULL,
                      show_mean        = FALSE,
                      show_agg_val     = TRUE,
                      horizontal       = FALSE,
                      bracket_spacing  = 0.08,
                      subgroup         = NULL,
                      ...) {

  # ---------------------------------------------------------------------------
  # 0. Package checks
  # ---------------------------------------------------------------------------
  check_and_install_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Package '", pkg, "' is not installed. Attempting to install..."))
      tryCatch(
        install.packages(pkg, repos = "https://cran.r-project.org"),
        error = function(e) stop(paste0("Failed to install '", pkg, "': ", e$message))
      )
      if (!requireNamespace(pkg, quietly = TRUE))
        stop(paste0("Package '", pkg, "' could not be loaded after installation."))
      message(paste0("Package '", pkg, "' installed successfully."))
    }
  }
  check_and_install_package("ggpubr")
  check_and_install_package("dplyr")
  check_and_install_package("ggplot2")
  check_and_install_package("ggrepel")
  check_and_install_package("rstatix")

  # ---------------------------------------------------------------------------
  # 1. Resolve font sizes
  #    Proportions relative to global_font_size:
  #      title      -> 1.20 x
  #      subtitle   -> 0.90 x
  #      xlab/ylab  -> 1.00 x
  #      axis_text  -> 0.85 x
  #      pvalue     -> 0.35 x  (ggplot2 'size' units ≈ pt/2.845)
  #    Hard defaults (used when global_font_size is also NULL):
  #      title=14, subtitle=11, xlab/ylab=12, axis_text=10, pvalue=3.5
  # ---------------------------------------------------------------------------
  .resolve_size <- function(user_val, gfs, prop, hard_default) {
    if (!is.null(user_val)) return(user_val)           # explicit override wins
    if (!is.null(gfs))      return(gfs * prop)         # scale from global
    return(hard_default)                               # built-in default
  }

  fs_title    <- .resolve_size(title_size,      global_font_size, 1.20, 14)
  fs_subtitle <- .resolve_size(subtitle_size,   global_font_size, 0.90, 11)
  fs_xlab     <- .resolve_size(xlab_size,       global_font_size, 1.00, 12)
  fs_ylab     <- .resolve_size(ylab_size,       global_font_size, 1.00, 12)
  fs_axis     <- .resolve_size(axis_text_size,  global_font_size, 0.85, 10)
  fs_pvalue   <- .resolve_size(pvalue_size,     global_font_size, 0.35, 3.5)
  # agg val text (mean/median labels) — slightly smaller than pvalue
  fs_agg      <- .resolve_size(NULL,            global_font_size, 0.30, 3.0)

  # ---------------------------------------------------------------------------
  # 2. p-value formatter
  # ---------------------------------------------------------------------------
  .fmt_p <- function(p_vec) {
    # Determine number of decimal places
    if (identical(format_p_numeric, "auto")) {
      decimals <- 3L
    } else {
      decimals <- as.integer(format_p_numeric)
    }
    fmt_str <- paste0("%.", decimals, "f")

    sapply(p_vec, function(pv) {
      if (is.na(pv)) return(NA_character_)
      rounded <- round(pv, decimals)
      if (rounded < 10^(-decimals)) {
        return(paste0("<.", paste(rep("0", decimals - 1), collapse = ""), "1"))
      }
      if (rounded > (1 - 10^(-decimals))) {
        return(paste0(">.", paste(rep("9", decimals), collapse = "")))
      }
      # Drop leading zero for APA style (e.g. ".045" not "0.045")
      val_str <- sprintf(fmt_str, rounded)
      sub("^0\\.", ".", val_str)
    })
  }

  # ---------------------------------------------------------------------------
  # 3. Helpers
  # ---------------------------------------------------------------------------
  get_signif_stars <- function(p) {
    sapply(p, function(pv) {
      if (is.na(pv))    return("ns")
      if (pv < 0.001)   return("***")
      if (pv < 0.01)    return("**")
      if (pv < 0.05)    return("*")
      return("ns")
    })
  }

  format_agg_val <- function(vals, label_prefix = "M") {
    vapply(vals, function(val) {
      if (is.na(val)) return(NA_character_)
      if (abs(round(val, 2)) < 0.01) {
        formatted <- gsub("E", "e", sprintf("%.2e", val))
      } else {
        formatted <- sprintf("%.2f", val)
      }
      paste0(label_prefix, "=", formatted)
    }, character(1))
  }

  default_colors <- c("#7BAFD4", "#9AC5A8", "#E8B4A0", "#C4A5D8",
                      "#F4C48D", "#A8D5BA", "#D4A5C0", "#B8C4D6",
                      "#C9B8A0", "#A5C4C4", "#D8B8C4", "#B4D4A5")
  plot_colors <- if (!is.null(colors)) colors else default_colors

  # ---------------------------------------------------------------------------
  # 4. Input validation
  # ---------------------------------------------------------------------------
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
  # 5. Main loop over plot variables
  # ---------------------------------------------------------------------------
  all_plots        <- list()
  significant_plots <- list()
  statistics       <- list()
  subgroup_analysis <- list()

  for (var in plot_vars) {
    if (!is.numeric(x[[var]])) {
      warning(paste("Variable", var, "is not numeric. Skipping."))
      next
    }

    stat_results <- NULL
    tryCatch({
      stat_results <- run_diff(
        x                     = x,
        outcome               = var,
        group                 = group,
        test_type             = if (is.null(test_type)) "auto" else test_type,
        p_adjust_method       = p_adjust_method,
        paired                = paired,
        test_alpha            = test_alpha,
        calculate_effect_size = TRUE,
        perform_posthoc       = TRUE,
        verbose               = FALSE,
        subgroup              = NULL
      )
    }, error = function(e) {
      warning(paste("Statistical test failed for", var, ":", e$message), call. = FALSE)
    }, warning = function(w) {
      warning(paste("Warning in run_diff for", var, ":", w$message), call. = FALSE)
    })

    if (is.null(stat_results)) next

    n_groups        <- nlevels(stat_results$raw_data$group)
    main_test_p_val <- stat_results$test_result$p.value
    posthoc_data    <- stat_results$posthoc_result
    is_parametric   <- isTRUE(stat_results$parametric)

    # Prepare plotting data
    plot_data           <- stat_results$raw_data
    colnames(plot_data) <- c(var, group)

    if (!is.null(sample_id)) {
      if (!sample_id %in% names(x))
        stop(paste("Sample ID variable", sample_id, "not found in data."))
      plot_data[[sample_id]] <- x[[sample_id]][complete.cases(x[c(var, group)])]
    }

    y_max       <- max(plot_data[[var]], na.rm = TRUE)
    y_min       <- min(plot_data[[var]], na.rm = TRUE)
    y_range     <- y_max - y_min
    y_max_plot  <- y_max + y_range * 0.15

    # -------------------------------------------------------------------------
    # 5a. Prepare post-hoc table for bracket annotations
    # -------------------------------------------------------------------------
    stat_result_plot <- data.frame()

    if (n_groups > 2 && !is.null(posthoc_data) && nrow(posthoc_data) > 0) {
      required_cols <- c("group1", "group2", "p.adj", "p.adj.signif")
      if (!all(required_cols %in% names(posthoc_data))) {
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
      p_signif <- get_signif_stars(main_test_p_val)
      if (!hide_ns || p_signif != "ns") {
        grp_names <- levels(plot_data[[group]])
        if (is.null(grp_names)) grp_names <- unique(as.character(plot_data[[group]]))
        stat_result_plot <- data.frame(
          group1      = grp_names[1],
          group2      = grp_names[2],
          p.adj       = main_test_p_val,
          p.adj.signif = p_signif,
          y.position  = y_max_plot + y_range * 0.05
        )
      }
    }

    # Add formatted p-value column for numeric display
    if (nrow(stat_result_plot) > 0 && show_p_numeric) {
      stat_result_plot$p.adj.fmt <- .fmt_p(stat_result_plot$p.adj)
    }

    is_significant <- if (nrow(stat_result_plot) > 0) {
      any(stat_result_plot$p.adj < test_alpha, na.rm = TRUE)
    } else {
      main_test_p_val < test_alpha
    }

    # -------------------------------------------------------------------------
    # 5b. Labels
    # -------------------------------------------------------------------------
    title_text    <- if (!is.null(plot_title))    plot_title    else paste(var, "by", group)
    subtitle_text <- if (!is.null(plot_subtitle)) plot_subtitle else stat_results$test_used
    # always use user-supplied xlab/ylab; fall back to var/group names
    xlab_text     <- if (!is.null(xlab)) xlab else group
    ylab_text     <- if (!is.null(ylab)) ylab else var

    # Column to use for bracket labels
    p_label_col <- if (show_p_numeric) "p.adj.fmt" else "p.adj.signif"

    # -------------------------------------------------------------------------
    # 5c. Base plot (labs() is authoritative for axis labels)
    # -------------------------------------------------------------------------
    base_theme <- ggplot2::theme(
      legend.position  = "none",
      plot.title       = ggplot2::element_text(face = "bold", size = fs_title),
      plot.subtitle    = ggplot2::element_text(size = fs_subtitle),
      axis.title.x     = ggplot2::element_text(size = fs_xlab),
      axis.title.y     = ggplot2::element_text(size = fs_ylab),
      axis.text        = ggplot2::element_text(size = fs_axis)
    )

    if (plot_type == "boxplot") {
      p <- ggpubr::ggboxplot(
        plot_data, x = group, y = var,
        color   = group, palette = plot_colors,
        add     = "jitter", order = group_order, ...
      )
    } else {
      p <- ggpubr::ggviolin(
        plot_data, x = group, y = var,
        color = group, fill = group, palette = plot_colors,
        add = "boxplot", add.params = list(fill = "white"),
        order = group_order, ...
      )
    }

    # Apply labs() AFTER ggpubr call so our labels override ggpubr defaults
    p <- p +
      ggplot2::labs(
        title    = title_text,
        subtitle = subtitle_text,
        x        = xlab_text,
        y        = ylab_text
      ) +
      base_theme

    # -------------------------------------------------------------------------
    # 5d. Mean dot
    # -------------------------------------------------------------------------
    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        fun   = mean, geom = "point",
        shape = 21, size = 3, fill = "red", color = "black"
      )
    }

    # -------------------------------------------------------------------------
    # 5e. Aggregate value text (mean / median)
    # -------------------------------------------------------------------------
    if (show_agg_val) {
      if (is_parametric) {
        p <- p + ggplot2::stat_summary(
          fun  = mean, geom = "text",
          aes(label = after_stat(format_agg_val(y, label_prefix = "M"))),
          vjust = -0.7, color = "red", size = fs_agg
        )
      } else {
        p <- p + ggplot2::stat_summary(
          fun  = median, geom = "text",
          aes(label = after_stat(format_agg_val(y, label_prefix = "Med"))),
          vjust = -0.7, color = "blue", size = fs_agg
        )
      }
    }

    # -------------------------------------------------------------------------
    # 5f. Point labels
    # -------------------------------------------------------------------------
    if (!is.null(label_points)) {
      if (label_points == "all") {
        p <- p + ggrepel::geom_text_repel(
          data = plot_data,
          aes_string(label = sample_id),
          size = fs_agg, max.overlaps = Inf
        )
      } else if (label_points == "outliers") {
        outliers <- rstatix::identify_outliers(plot_data, var)
        if (nrow(outliers) > 0) {
          outlier_df <- outliers %>% dplyr::left_join(plot_data, by = c(var, group))
          p <- p + ggrepel::geom_text_repel(
            data = outlier_df,
            aes_string(x = group, y = var, label = sample_id),
            size = fs_agg, color = "red", max.overlaps = Inf
          )
        }
      }
    }

    # -------------------------------------------------------------------------
    # 5g. P-value brackets (pairwise)
    # -------------------------------------------------------------------------
    if (nrow(stat_result_plot) > 0) {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_result_plot,
        label      = p_label_col,
        hide.ns    = FALSE,
        tip.length = 0.01,
        size       = fs_pvalue
      )
    }

    # -------------------------------------------------------------------------
    # 5h. Global p-value annotation (show_global_p)
    #     Previously this was always shown for n_groups > 2; now it is opt-in
    #     via show_global_p and works for any number of groups.
    # -------------------------------------------------------------------------
    if (show_global_p) {
      global_y <- if (!is.null(global_p_position)) {
        global_p_position
      } else {
        y_max_plot + y_range * 0.25
      }
      global_x <- if (!is.null(global_p_x)) {
        global_p_x
      } else {
        n_groups * 0.95
      }

      # Format using the same formatter as bracket p-values
      if (show_p_numeric || identical(format_p_numeric, "auto") || is.numeric(format_p_numeric)) {
        gp_label <- paste0("p = ", .fmt_p(main_test_p_val))
      } else {
        gp_label <- paste0("p = ", .fmt_p(main_test_p_val))
      }

      # Prepend the test name for clarity
      gp_label <- paste0(stat_results$test_used, "\n", gp_label)

      p <- p + ggplot2::annotate(
        "text",
        x     = global_x,
        y     = global_y,
        label = gp_label,
        hjust = 1,
        size  = fs_pvalue
      )
    }

    if (horizontal) p <- p + ggplot2::coord_flip()

    # -------------------------------------------------------------------------
    # 5i. Store
    # -------------------------------------------------------------------------
    all_plots[[var]]  <- p
    statistics[[var]] <- stat_results
    if (is_significant) significant_plots[[var]] <- p
  }

  # ---------------------------------------------------------------------------
  # 6. Subgroup analysis
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
          warning(paste0("Subgroup '", sg_var, "' level '", lvl, "' has no observations. Skipping."))
          next
        }

        sg_all_plots         <- list()
        sg_significant_plots <- list()
        sg_statistics        <- list()

        for (var in plot_vars) {
          if (!is.numeric(x_sub[[var]])) next
          grp_sub <- droplevels(as.factor(x_sub[[group]]))
          if (nlevels(grp_sub) < 2) {
            warning(paste0("Subgroup '", sg_var, "' = '", lvl, "': variable '", var,
                           "' has fewer than 2 group levels. Skipping."))
            next
          }

          sg_stat <- NULL
          tryCatch({
            sg_stat <- run_diff(
              x                     = x_sub,
              outcome               = var,
              group                 = group,
              test_type             = if (is.null(test_type)) "auto" else test_type,
              p_adjust_method       = p_adjust_method,
              paired                = paired,
              test_alpha            = test_alpha,
              calculate_effect_size = TRUE,
              perform_posthoc       = TRUE,
              group_order           = group_order,
              verbose               = FALSE,
              subgroup              = NULL
            )
          }, error = function(e) {
            warning(paste0("run_diff failed for subgroup '", sg_var, "' = '", lvl,
                           "', variable '", var, "': ", e$message), call. = FALSE)
          })

          if (is.null(sg_stat)) next

          sg_plot_result <- tryCatch({
            plot_diff(
              x                = x_sub,
              outcome          = var,
              group            = group,
              plot_vars        = var,
              plot_type        = plot_type,
              test_type        = if (is.null(test_type)) "auto" else test_type,
              posthoc          = posthoc,
              p_adjust_method  = p_adjust_method,
              hide_ns          = hide_ns,
              show_p_numeric   = show_p_numeric,
              paired           = paired,
              test_alpha       = test_alpha,
              colors           = colors,
              plot_title       = if (!is.null(plot_title)) plot_title else
                paste0(var, " by ", group, "\n[", sg_var, " = ", lvl, "]"),
              plot_subtitle    = plot_subtitle,
              xlab             = xlab,
              ylab             = ylab,
              show_global_p    = show_global_p,
              global_p_position = global_p_position,
              global_p_x       = global_p_x,
              global_font_size = global_font_size,
              title_size       = title_size,
              subtitle_size    = subtitle_size,
              xlab_size        = xlab_size,
              ylab_size        = ylab_size,
              pvalue_size      = pvalue_size,
              axis_text_size   = axis_text_size,
              format_p_numeric = format_p_numeric,
              group_order      = group_order,
              label_points     = label_points,
              sample_id        = sample_id,
              show_mean        = show_mean,
              show_agg_val     = show_agg_val,
              horizontal       = horizontal,
              bracket_spacing  = bracket_spacing,
              subgroup         = NULL  # prevent infinite recursion
            )
          }, error = function(e) {
            warning(paste0("plot_diff failed for subgroup '", sg_var, "' = '", lvl,
                           "', variable '", var, "': ", e$message), call. = FALSE)
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

  return(list(
    plots             = all_plots,
    significant_plots = significant_plots,
    statistics        = statistics,
    subgroup_analysis = if (length(subgroup_analysis) > 0) subgroup_analysis else NULL
  ))
}