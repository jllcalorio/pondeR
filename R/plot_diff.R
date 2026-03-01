#' Create Box or Violin Plots with Statistical Comparisons
#'
#' This function generates box plots or violin plots with statistical tests and
#' pairwise comparisons. It automatically selects appropriate statistical tests
#' based on assumptions, or allows manual specification. Statistical testing is
#' performed using the companion function \code{\link{run_diff}}.
#'
#' @param x A data frame containing the variables to plot.
#' @param outcome Character string specifying the column name for the numeric outcome variable. Used as the default for \code{plot_vars} if it is not specified.
#' @param group Character string specifying the column name for the grouping variable.
#' @param plot_vars Character vector specifying column name(s) for numeric variables to plot. If \code{NULL} (default), uses the value of \code{outcome}.
#' @param plot_type Character string specifying plot type: "boxplot" or "violin". Default is "boxplot".
#' @param test_type Character string or NULL. Statistical test strategy to use: "auto", "parametric", "nonparametric". See \code{\link{run_diff}} for details. If NULL (default), uses "auto".
#' @param posthoc Character string or NULL. Post-hoc test for multiple groups:
#'   "auto", "none", "tukey", "games-howell", "dunn", "t_test", or "wilcox_test". See \code{\link{run_diff}} for details. If NULL (default), uses "auto".
#' @param p_adjust_method Character string. P-value adjustment method for multiple comparisons.
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'   Default is "BH".
#' @param hide_ns Logical. If TRUE, hide non-significant p-values from plot brackets. Default is TRUE.
#' @param show_p_numeric Logical. If TRUE, show exact p-values as numbers instead of
#'   significance symbols (*, **, ***). Default is FALSE.
#' @param paired Logical. If TRUE, use paired tests. Default is FALSE.
#' @param test_alpha Numeric. Significance level for determining significance on the plot. Default is 0.05.
#' @param colors Character vector of colors (names or HEX codes) or NULL for automatic soft colors.
#'   Can be a single color or a vector of colors. Default is NULL.
#' @param plot_title Character string or NULL. Custom title for the plot. If NULL, generates
#'   an automatic title. Default is NULL.
#' @param xlab Character string or NULL. X-axis label. If NULL, uses group variable name. Default is NULL.
#' @param ylab Character string or NULL. Y-axis label. If NULL, uses plot variable name. Default is NULL.
#' @param global_p_position Numeric or NULL. Y-axis position for the global p-value as data coordinates.
#'   If NULL (default), automatically positions at the top-right of the plot area.
#' @param global_p_x Numeric or NULL. X-axis position for the global p-value. If NULL (default),
#'   positions at 95% of the x-axis range (right side). Only used when n_groups > 2.
#' @param title_size Numeric. Font size for the plot title. Default is 12.
#' @param subtitle_size Numeric. Font size for the plot subtitle. Default is 10.
#' @param xlab_size Numeric. Font size for the x-axis label. Default is 11.
#' @param ylab_size Numeric. Font size for the y-axis label. Default is 11.
#' @param axis_text_size Numeric. Font size for axis text. Default is 10.
#' @param group_order Character vector or NULL. Specifies the order of groups on the x-axis.
#'   If NULL (default), uses factor levels or alphabetical order.
#' @param label_points Character string or NULL. Specifies point labeling: "all" (label all jitter points),
#'   "outliers" (label only outliers), or NULL (default, no labels). Requires sample_id parameter.
#' @param sample_id Character string or NULL. Column name in data containing sample identifiers
#'   for point labeling. Required when label_points is not NULL.
#' @param show_mean Logical. If TRUE, add a red dot (mean indicator) to plots. Default is FALSE.
#' @param show_agg_val Logical. If TRUE (default), shows the numeric value of the mean
#'   (for parametric tests) or median (for non-parametric tests) on the plot.
#' @param horizontal Logical. If TRUE, flip the plot to horizontal orientation. Default is FALSE.
#' @param bracket_spacing Numeric. Multiplier for spacing between comparison brackets. Default is 0.08.
#' @param ... Additional arguments passed to ggboxplot or ggviolin from ggpubr.
#'
#' @return A list containing:
#'   \item{plots}{List of all generated plots.}
#'   \item{significant_plots}{List of plots with significant results only.}
#'   \item{statistics}{List of statistical test results for each variable (output from \code{\link{run_diff}}).}
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
#'   x = PlantGrowth,
#'   outcome = "weight",
#'   group = "group",
#'   plot_type = "boxplot",
#'   show_mean = TRUE,
#'   show_agg_val = TRUE
#' )
#' print(results_anova$plots$weight)
#'
#' # Example 2: Paired t-test, show only mean dot
#' sleep_data <- sleep
#' sleep_data$ID <- rep(1:10, 2) # Create a patient ID column
#'
#' results_paired <- plot_diff(
#'   x = sleep_data,
#'   outcome = "extra",
#'   group = "group",
#'   paired = TRUE,
#'   show_mean = TRUE,
#'   show_agg_val = FALSE
#' )
#' print(results_paired$plots$extra)
#' }
#'
#' @export
plot_diff <- function(x,
                      outcome,
                      group,
                      plot_vars = NULL,
                      plot_type = "boxplot",
                      test_type = NULL,
                      posthoc = NULL,
                      p_adjust_method = "BH",
                      hide_ns = TRUE,
                      show_p_numeric = FALSE,
                      paired = FALSE,
                      test_alpha = 0.05,
                      colors = NULL,
                      plot_title = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      global_p_position = NULL,
                      global_p_x = NULL,
                      title_size = 12,
                      subtitle_size = 10,
                      xlab_size = 11,
                      ylab_size = 11,
                      axis_text_size = 10,
                      group_order = NULL,
                      label_points = NULL,
                      sample_id = NULL,
                      show_mean = FALSE,
                      show_agg_val = TRUE,
                      horizontal = FALSE,
                      bracket_spacing = 0.08,
                      ...) {

  # --- Check for and Install Required Packages ---
  check_and_install_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package '", pkg, "' is not installed. Attempting to install it...", sep = ""))
      tryCatch({
        install.packages(pkg, repos = "https://cran.r-project.org")
      }, error = function(e) {
        stop(paste("Failed to install package '", pkg, "'. Please install it manually. Error: ", e$message, sep = ""))
      })
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Package '", pkg, "' could not be loaded after installation.", sep = ""))
      }
      message(paste("Package '", pkg, "' installed and loaded successfully.", sep = ""))
    }
  }

  check_and_install_package("ggpubr")
  check_and_install_package("dplyr")
  check_and_install_package("ggplot2")
  check_and_install_package("ggrepel")
  check_and_install_package("rstatix")

  # --- Helper function for significance stars ---
  get_signif_stars <- function(p) {
    sapply(p, function(p_val) {
      if (is.na(p_val)) return("ns")
      if (p_val < 0.001) return("***")
      if (p_val < 0.01) return("**")
      if (p_val < 0.05) return("*")
      return("ns")
    })
  }

  default_colors <- c("#7BAFD4", "#9AC5A8", "#E8B4A0", "#C4A5D8",
                      "#F4C48D", "#A8D5BA", "#D4A5C0", "#B8C4D6",
                      "#C9B8A0", "#A5C4C4", "#D8B8C4", "#B4D4A5")
  plot_colors <- if (!is.null(colors)) colors else default_colors

  # --- Input Validation ---
  if (!all(c(outcome, group) %in% names(x))) stop(paste("Outcome or group variable not found in x"))

  # --- Set plot_vars if not provided ---
  if (is.null(plot_vars)) {
    plot_vars <- outcome
  }
  if (!all(plot_vars %in% names(x))) stop(paste("Variable(s) to plot not found in data:", paste(setdiff(plot_vars, names(x)), collapse=", ")))

  # --- Initialize Output Lists ---
  all_plots <- list()
  significant_plots <- list()
  statistics <- list()

  # --- Process Each Plot Variable ---
  for (var in plot_vars) {
    if (!is.numeric(x[[var]])) { warning(paste("Variable", var, "is not numeric. Skipping.")); next }

    stat_results <- NULL

    # --- FIX: Improve error handling for run_diff ---
    tryCatch({
      stat_results <- run_diff(
        x = x,
        outcome = var,
        group = group,
        test_type = if (is.null(test_type)) "auto" else test_type,
        p_adjust_method = p_adjust_method,
        paired = paired,
        test_alpha = test_alpha,
        calculate_effect_size = TRUE,
        perform_posthoc = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      warning(paste("Statistical test failed for", var, ":", e$message), call. = FALSE)
    }, warning = function(w) {
      warning(paste("Warning in run_diff for", var, ":", w$message), call. = FALSE)
    })

    if (is.null(stat_results)) { next }

    # --- Reuse results from run_diff ---
    n_groups <- nlevels(stat_results$raw_data$group)
    main_test_p_val <- stat_results$test_result$p.value
    posthoc_data <- stat_results$posthoc_result

    # Prepare plotting data
    plot_data <- stat_results$raw_data
    colnames(plot_data) <- c(var, group) # Rename for ggpubr

    if (!is.null(sample_id)) {
      if (!sample_id %in% names(x)) stop(paste("Sample ID variable", sample_id, "not found in data."))
      original_ids <- x[[sample_id]][complete.cases(x[c(var, group)])]
      plot_data[[sample_id]] <- original_ids
    }

    # Use data range for brackets
    y_max <- max(plot_data[[var]], na.rm = TRUE)
    y_min <- min(plot_data[[var]], na.rm = TRUE)
    y_range <- y_max - y_min
    y_max_plot <- y_max + y_range * 0.15  # Base buffer

    # --- Prepare Post-Hoc Data for Plotting ---
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
        if (hide_ns) {
          stat_result_plot <- dplyr::filter(stat_result_plot, p.adj < test_alpha)
        }
      }
    }
    else if (n_groups == 2) {
      p_signif <- get_signif_stars(main_test_p_val)
      if (!hide_ns || (hide_ns && p_signif != "ns")) {
        group_names <- levels(plot_data[[group]])
        if (is.null(group_names)) group_names <- unique(as.character(plot_data[[group]]))

        stat_result_plot <- data.frame(
          group1 = group_names[1],
          group2 = group_names[2],
          p.adj = main_test_p_val,
          p.adj.signif = p_signif,
          y.position = y_max_plot + y_range * 0.05
        )
      }
    }

    is_significant <- if (nrow(stat_result_plot) > 0) {
      any(stat_result_plot$p.adj < test_alpha, na.rm = TRUE)
    } else {
      main_test_p_val < test_alpha
    }

    # --- Create Base Plot ---
    title_text <- if (!is.null(plot_title)) plot_title else paste(var, "by", group)
    xlab_text <- if (!is.null(xlab)) xlab else group
    ylab_text <- if (!is.null(ylab)) ylab else var
    p_label_col <- if (show_p_numeric) "p.adj" else "p.adj.signif"
    subtitle_text <- stat_results$test_used

    if (plot_type == "boxplot") {
      p <- ggpubr::ggboxplot(plot_data, x = group, y = var,
                             color = group, palette = plot_colors,
                             add = "jitter", order = group_order, ...) +
        ggplot2::labs(title = title_text, subtitle = subtitle_text,
                      x = xlab_text, y = ylab_text) +
        ggplot2::theme(legend.position = "none",
                       plot.title = ggplot2::element_text(face = "bold", size = title_size),
                       plot.subtitle = ggplot2::element_text(size = subtitle_size),
                       axis.title.x = ggplot2::element_text(size = xlab_size),
                       axis.title.y = ggplot2::element_text(size = ylab_size),
                       axis.text = ggplot2::element_text(size = axis_text_size))
    } else {
      p <- ggpubr::ggviolin(plot_data, x = group, y = var,
                            color = group, fill = group, palette = plot_colors,
                            add = "boxplot", add.params = list(fill = "white"), order = group_order, ...) +
        ggplot2::labs(title = title_text, subtitle = subtitle_text,
                      x = xlab_text, y = ylab_text) +
        ggplot2::theme(legend.position = "none",
                       plot.title = ggplot2::element_text(face = "bold", size = title_size),
                       plot.subtitle = ggplot2::element_text(size = subtitle_size),
                       axis.title.x = ggplot2::element_text(size = xlab_size),
                       axis.title.y = ggplot2::element_text(size = ylab_size),
                       axis.text = ggplot2::element_text(size = axis_text_size))
    }

    # --- Helper function to format mean/median values (vectorized) ---
    format_agg_val <- function(x, label_prefix = "M") {
      vapply(x, function(val) {
        if (is.na(val)) return(NA_character_)
        if (abs(round(val, 2)) < 0.01) { #} && val != 0) {
          formatted <- sprintf("%.2e", val)
          formatted <- gsub("E", "e", formatted)
        } else {
          formatted <- sprintf("%.2f", val)
        }
        paste0(label_prefix, "=", formatted)
      }, character(1))
    }

    # --- *** ADDED SECTION for show_mean *** ---
    # Add mean dot (if requested)
    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        fun = mean,
        geom = "point",
        shape = 21, # Circle shape
        size = 3,
        fill = "red",
        color = "black"
      )
    }

    # --- Add Mean/Median Numeric Value ---
    if (show_agg_val) {
      is_parametric <- grepl("ANOVA|t-test|Parametric|Tukey|Games-Howell", subtitle_text, ignore.case = TRUE)

      if (is_parametric) {
        # Add mean value (with scientific notation if very small)
        p <- p + ggplot2::stat_summary(
          fun = mean,
          geom = "text",
          aes(label = after_stat(format_agg_val(y, label_prefix = "M"))),
          vjust = -0.7,
          color = "red",
          size = 3
        )
      } else {
        # Add median value (with scientific notation if very small)
        p <- p + ggplot2::stat_summary(
          fun = median,
          geom = "text",
          aes(label = after_stat(format_agg_val(y, label_prefix = "Med"))),
          vjust = -0.7,
          color = "blue",
          size = 3
        )
      }
    }

    # --- Add Point Labels (All or Outliers) ---
    if (!is.null(label_points)) {
      if (label_points == "all") {
        p <- p + ggrepel::geom_text_repel(data = plot_data, aes_string(label = sample_id), size = 3, max.overlaps = Inf)
      } else if (label_points == "outliers") {

        outliers <- rstatix::identify_outliers(plot_data, var)

        if (nrow(outliers) > 0) {
          outlier_df <- outliers %>%
            dplyr::left_join(plot_data, by = c(var, group))

          p <- p + ggrepel::geom_text_repel(
            data = outlier_df,
            aes_string(x = group, y = var, label = sample_id),
            size = 3,
            color = "red",
            max.overlaps = Inf
          )
        }
      }
    }

    # --- Add P-value Brackets ---
    if (nrow(stat_result_plot) > 0) {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_result_plot,
        label = p_label_col,
        hide.ns = FALSE,
        tip.length = 0.01
      )
    }

    # --- Add Global P-value (only for n_groups > 2) ---
    if (n_groups > 2) {
      global_y <- if (is.null(global_p_position)) y_max_plot + y_range * 0.25 else global_p_position
      global_x <- if (is.null(global_p_x)) n_groups * 0.95 else global_p_x
      p_label <- if (main_test_p_val < 0.001) paste0("p < 0.001") else paste0("p = ", format.pval(main_test_p_val, digits = 3))
      p <- p + ggplot2::annotate("text", x = global_x, y = global_y, label = p_label, hjust = 1, size = 3.5)
    }

    if (horizontal) p <- p + ggplot2::coord_flip()

    # --- Store Results ---
    all_plots[[var]] <- p
    statistics[[var]] <- stat_results
    if (is_significant) significant_plots[[var]] <- p
  }

  return(list(plots = all_plots, significant_plots = significant_plots, statistics = statistics))
}
