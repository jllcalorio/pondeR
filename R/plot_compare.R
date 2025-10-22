#' Create Box or Violin Plots with Statistical Comparisons
#'
#' This function generates box plots or violin plots with statistical tests and
#' pairwise comparisons. It automatically selects appropriate statistical tests
#' based on sample size and number of groups, or allows manual specification.
#' Statistical testing is performed using the companion function `compare_stats()`.
#'
#' @param data A data frame containing the variables to plot
#' @param group_var Character string specifying the column name for grouping variable
#' @param plot_vars Character vector specifying the column name(s) for numeric variables to plot
#' @param plot_type Character string specifying plot type: "boxplot" or "violin". Default is "boxplot"
#' @param test_method Character string or NULL. Statistical test to use: "t.test", "wilcox.test",
#'   "anova", "kruskal.test", "paired.t.test", "paired.wilcox.test". If NULL (default),
#'   automatically selects based on number of groups and sample size
#' @param posthoc_test Character string or NULL. Post-hoc test for multiple groups:
#'   "tukey", "dunn", "wilcox", "t.test". If NULL (default), automatically selects
#'   based on the main test used
#' @param p_adjust_method Character string. P-value adjustment method for multiple comparisons.
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'   Default is "BH"
#' @param hide_ns Logical. If TRUE, hide non-significant p-values. Default is TRUE
#' @param show_p_numeric Logical. If TRUE, show exact p-values as numbers instead of
#'   significance symbols (*, **, ***). Default is FALSE
#' @param paired Logical. If TRUE, use paired tests. Default is FALSE
#' @param alpha Numeric. Significance level for determining significance. Default is 0.05
#' @param colors Character vector of colors (names or HEX codes) or NULL for automatic soft colors.
#'   Can be a single color or a vector of colors. Default is NULL
#' @param plot_title Character string or NULL. Custom title for the plot. If NULL, generates
#'   automatic title. Default is NULL
#' @param xlab Character string or NULL. X-axis label. If NULL, uses group variable name. Default is NULL
#' @param ylab Character string or NULL. Y-axis label. If NULL, uses plot variable name. Default is NULL
#' @param global_p_position Numeric or NULL. Y-axis position for global p-value as data coordinates.
#'   If NULL (default), automatically positions at top-right of plot area.
#' @param global_p_x Numeric or NULL. X-axis position for global p-value. If NULL (default),
#'   positions at 95% of x-axis range (right side). Only used when n_groups > 2.
#' @param title_size Numeric. Font size for plot title. Default is 12
#' @param subtitle_size Numeric. Font size for plot subtitle. Default is 10
#' @param xlab_size Numeric. Font size for x-axis label. Default is 11
#' @param ylab_size Numeric. Font size for y-axis label. Default is 11
#' @param axis_text_size Numeric. Font size for axis text. Default is 10
#' @param group_order Character vector or NULL. Specifies the order of groups on the x-axis.
#'   If NULL (default), uses the factor levels or alphabetical order.
#' @param label_points Character string or NULL. Specifies point labeling: "all" (label all jitter points),
#'   "outliers" (label only outliers), or NULL (default, no labels). Requires sample_id parameter.
#' @param sample_id Character string or NULL. Column name in data containing sample identifiers
#'   for point labeling. Required when label_points is not NULL.
#' @param show_mean Logical. If TRUE, add mean indicator to plots. Default is FALSE
#' @param horizontal Logical. If TRUE, flip plot to horizontal orientation. Default is FALSE
#' @param bracket_spacing Numeric. Multiplier for spacing between comparison brackets. Default is 0.08
#' @param ... Additional arguments passed to `ggboxplot` or `ggviolin`
#'
#' @return A list containing:
#'   \item{plots}{List of all generated plots}
#'   \item{significant_plots}{List of plots with significant results only}
#'   \item{statistics}{List of statistical test results for each variable (output from compare_stats)}
#'
#' @import ggpubr
#' @import dplyr
#' @import rstatix
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats IQR quantile
#'
#' @examples
#' \dontrun{
#' # Example data
#' set.seed(123)
#' df <- data.frame(
#'   SampleID = paste0("Sample", 1:60),
#'   Group = rep(c("A", "B", "C"), each = 20),
#'   ShannonDiversity = c(rnorm(20, 3, 0.5), rnorm(20, 3.5, 0.5), rnorm(20, 4, 0.5)),
#'   SimpsonsIndex = c(rnorm(20, 0.7, 0.1), rnorm(20, 0.75, 0.1), rnorm(20, 0.8, 0.1))
#' )
#'
#' # Basic usage
#' results <- plot_compare(
#'   data = df,
#'   group_var = "Group",
#'   plot_vars = c("ShannonDiversity", "SimpsonsIndex")
#' )
#'
#' # With custom group order and labeled outliers
#' results <- plot_compare(
#'   data = df,
#'   group_var = "Group",
#'   plot_vars = c("ShannonDiversity", "SimpsonsIndex"),
#'   plot_type = "violin",
#'   group_order = c("C", "B", "A"),
#'   label_points = "outliers",
#'   sample_id = "SampleID",
#'   show_mean = TRUE,
#'   horizontal = TRUE
#' )
#'
#' # View all plots
#' results$plots
#'
#' # View only significant plots
#' results$significant_plots
#'
#' # Access statistical results
#' results$statistics$ShannonDiversity$pairwise_results
#' }
#'
#' @export
plot_compare <- function(data,
                         group_var,
                         plot_vars,
                         plot_type = "boxplot",
                         test_method = NULL,
                         posthoc_test = NULL,
                         p_adjust_method = "BH",
                         hide_ns = TRUE,
                         show_p_numeric = FALSE,
                         paired = FALSE,
                         alpha = 0.05,
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
                         horizontal = FALSE,
                         bracket_spacing = 0.08,
                         ...) {

  # Load required packages
  require(ggpubr)
  require(dplyr)
  require(rstatix)
  require(ggplot2)
  require(ggrepel)

  # Define soft, eye-friendly color palette
  default_colors <- c("#7BAFD4", "#9AC5A8", "#E8B4A0", "#C4A5D8",
                      "#F4C48D", "#A8D5BA", "#D4A5C0", "#B8C4D6",
                      "#C9B8A0", "#A5C4C4", "#D8B8C4", "#B4D4A5")

  # Use custom colors if provided, otherwise use default
  plot_colors <- if (!is.null(colors)) colors else default_colors

  # Validate inputs
  if (!group_var %in% names(data)) {
    stop(paste("Group variable", group_var, "not found in data"))
  }

  if (!all(plot_vars %in% names(data))) {
    missing_vars <- plot_vars[!plot_vars %in% names(data)]
    stop(paste("Variable(s) not found in data:", paste(missing_vars, collapse = ", ")))
  }

  if (!plot_type %in% c("boxplot", "violin")) {
    stop("plot_type must be either 'boxplot' or 'violin'")
  }

  # Validate label_points parameter
  if (!is.null(label_points) && !label_points %in% c("all", "outliers")) {
    stop("label_points must be either 'all', 'outliers', or NULL")
  }

  # Validate sample_id parameter
  if (!is.null(label_points) && is.null(sample_id)) {
    stop("sample_id parameter is required when label_points is not NULL")
  }

  if (!is.null(sample_id) && !sample_id %in% names(data)) {
    stop(paste("Sample ID variable", sample_id, "not found in data"))
  }

  # Convert group variable to factor if not already
  data[[group_var]] <- as.factor(data[[group_var]])

  # Apply group order if specified
  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(data[[group_var]]))) {
      stop("group_order contains levels not present in the data")
    }
    data[[group_var]] <- factor(data[[group_var]], levels = group_order)
  }

  # Get number of groups
  n_groups <- length(unique(data[[group_var]]))

  # Initialize output lists
  all_plots <- list()
  significant_plots <- list()
  statistics <- list()

  # Process each plot variable
  for (var in plot_vars) {

    # Check if variable is numeric
    if (!is.numeric(data[[var]])) {
      warning(paste("Variable", var, "is not numeric. Skipping."))
      next
    }

    # Perform statistical comparison using compare_stats
    tryCatch({
      stat_results <- compare_stats(
        data = data,
        group_var = group_var,
        test_var = var,
        test_method = test_method,
        posthoc_test = posthoc_test,
        p_adjust_method = p_adjust_method,
        paired = paired,
        alpha = alpha,
        group_order = group_order
      )
    }, error = function(e) {
      warning(paste("Statistical test failed for", var, ":", e$message))
      return(NULL)
    })

    if (is.null(stat_results)) {
      next
    }

    # Extract statistical results
    selected_test <- stat_results$test_used
    selected_posthoc <- stat_results$posthoc_used
    overall_p <- stat_results$overall_p
    stat_result <- stat_results$pairwise_results
    n_groups <- stat_results$n_groups
    sample_sizes <- stat_results$sample_sizes

    # Prepare plotting data
    plot_data <- data %>%
      dplyr::select(dplyr::all_of(c(group_var, var, if (!is.null(sample_id)) sample_id else NULL))) %>%
      dplyr::filter(!is.na(.data[[var]]))

    # Apply group order if specified
    if (!is.null(group_order)) {
      plot_data[[group_var]] <- factor(plot_data[[group_var]], levels = group_order)
    }

    # Calculate y-axis range for bracket positioning
    y_max <- max(plot_data[[var]], na.rm = TRUE)
    y_min <- min(plot_data[[var]], na.rm = TRUE)
    y_range <- y_max - y_min

    # Calculate y positions for brackets with proper spacing
    stat_result <- stat_result %>%
      dplyr::arrange(desc(group1), desc(group2)) %>%
      dplyr::mutate(
        y.position = y_max + y_range * (0.05 + bracket_spacing * (dplyr::row_number() - 1))
      )

    # Filter non-significant comparisons if requested
    if (hide_ns) {
      # For 2 groups, use 'p' column; for >2 groups, use 'p.adj'
      p_col <- if (n_groups == 2) "p" else "p.adj"
      stat_result_plot <- stat_result %>% filter(.data[[p_col]] < alpha)
    } else {
      stat_result_plot <- stat_result
    }

    # Determine if results are significant
    p_col_sig <- if (n_groups == 2) "p" else "p.adj"
    is_significant <- any(stat_result[[p_col_sig]] < alpha, na.rm = TRUE)

    # Apply custom labels
    title_text <- if (!is.null(plot_title)) plot_title else paste(var, "by", group_var)
    xlab_text <- if (!is.null(xlab)) xlab else group_var
    ylab_text <- if (!is.null(ylab)) ylab else var

    # Determine label to use based on show_p_numeric parameter
    p_label_col <- if (show_p_numeric) {
      if (n_groups == 2) "p.format" else "p.adj.format"
    } else {
      "p.adj.signif"
    }

    # Identify outliers if needed
    outliers_data <- NULL
    if (!is.null(label_points) && label_points == "outliers") {
      outliers_data <- plot_data %>%
        group_by(.data[[group_var]]) %>%
        mutate(
          Q1 = quantile(.data[[var]], 0.25, na.rm = TRUE),
          Q3 = quantile(.data[[var]], 0.75, na.rm = TRUE),
          IQR_val = IQR(.data[[var]], na.rm = TRUE),
          lower_bound = Q1 - 1.5 * IQR_val,
          upper_bound = Q3 + 1.5 * IQR_val,
          is_outlier = .data[[var]] < lower_bound | .data[[var]] > upper_bound
        ) %>%
        ungroup() %>%
        filter(is_outlier)
    }

    # Create base plot
    if (plot_type == "boxplot") {
      p <- ggboxplot(plot_data, x = group_var, y = var,
                     color = group_var, palette = plot_colors,
                     add = "jitter", ...) +
        labs(title = title_text,
             subtitle = paste("Test:", selected_test,
                              ifelse(n_groups > 2, paste("| Post-hoc:", selected_posthoc), "")),
             x = xlab_text,
             y = ylab_text) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", size = title_size),
              plot.subtitle = element_text(size = subtitle_size),
              axis.title.x = element_text(size = xlab_size),
              axis.title.y = element_text(size = ylab_size),
              axis.text = element_text(size = axis_text_size))
    } else {
      p <- ggviolin(plot_data, x = group_var, y = var,
                    color = group_var, fill = group_var, palette = plot_colors,
                    add = "boxplot", add.params = list(fill = "white"), ...) +
        labs(title = title_text,
             subtitle = paste("Test:", selected_test,
                              ifelse(n_groups > 2, paste("| Post-hoc:", selected_posthoc), "")),
             x = xlab_text,
             y = ylab_text) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", size = title_size),
              plot.subtitle = element_text(size = subtitle_size),
              axis.title.x = element_text(size = xlab_size),
              axis.title.y = element_text(size = ylab_size),
              axis.text = element_text(size = axis_text_size))
    }

    # Add mean indicator if requested
    if (show_mean) {
      mean_data <- plot_data %>%
        group_by(.data[[group_var]]) %>%
        summarise(mean_val = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

      p <- p + stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                            fill = "red", color = "black")
    }

    # Add point labels if requested
    if (!is.null(label_points)) {
      if (label_points == "all") {
        # Label all points
        p <- p + geom_text_repel(
          data = plot_data,
          aes(label = .data[[sample_id]]),
          size = 3,
          max.overlaps = Inf
        )
      } else if (label_points == "outliers" && !is.null(outliers_data) && nrow(outliers_data) > 0) {
        # Label only outliers
        p <- p + geom_text_repel(
          data = outliers_data,
          aes(label = .data[[sample_id]]),
          size = 3,
          max.overlaps = Inf,
          color = "red"
        )
      }
    }

    # Add pairwise comparison p-values using the exact calculated values
    if (nrow(stat_result_plot) > 0) {
      p <- p + stat_pvalue_manual(
        stat_result_plot,
        label = p_label_col,
        hide.ns = hide_ns,
        tip.length = 0.01
      )
    }

    # Add overall p-value for multiple groups
    if (n_groups > 2) {
      # Determine position for global p-value
      if (is.null(global_p_position)) {
        # Auto-calculate: top of plot with some padding
        global_y <- y_max + y_range * 0.35
      } else {
        global_y <- global_p_position
      }

      if (is.null(global_p_x)) {
        # Position at right side (95% of x-axis)
        global_x <- n_groups * 0.95
      } else {
        global_x <- global_p_x
      }

      # Format p-value text
      p_label <- ifelse(overall_p < 0.001,
                        paste0(selected_test, ", p < 0.001"),
                        paste0(selected_test, ", p = ", format.pval(overall_p, digits = 3)))

      p <- p + annotate("text",
                        x = global_x,
                        y = global_y,
                        label = p_label,
                        hjust = 1,
                        size = 3.5)
    }

    # Flip coordinates if horizontal orientation requested
    if (horizontal) {
      p <- p + coord_flip()
    }

    # Store results
    all_plots[[var]] <- p
    statistics[[var]] <- list(
      test = selected_test,
      posthoc = if (n_groups > 2) selected_posthoc else NA,
      overall_p = if (n_groups > 2) overall_p else NA,
      results = stat_result,
      n_groups = n_groups,
      sample_sizes = sample_sizes
    )

    # Only add to significant_plots if there are significant results
    if (is_significant) {
      significant_plots[[var]] <- p
    }
  }

  # Return results
  return(list(
    plots = all_plots,
    significant_plots = significant_plots,
    statistics = statistics
  ))
}
