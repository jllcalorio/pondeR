#' Perform Statistical Comparisons Between Groups
#'
#' This function performs statistical tests and pairwise comparisons between groups.
#' It automatically selects appropriate statistical tests based on sample size and
#' number of groups, or allows manual specification. This function is designed to be
#' used standalone or called by plotting functions.
#'
#' @param data A data frame containing the variables to test
#' @param group_var Character string specifying the column name for grouping variable
#' @param test_var Character string specifying the column name for numeric variable to test
#' @param test_method Character string or NULL. Statistical test to use: "t.test", "wilcox.test",
#'   "anova", "kruskal.test", "paired.t.test", "paired.wilcox.test". If NULL (default),
#'   automatically selects based on number of groups and sample size
#' @param posthoc_test Character string or NULL. Post-hoc test for multiple groups:
#'   "tukey", "dunn", "wilcox", "t.test". If NULL (default), automatically selects
#'   based on the main test used
#' @param p_adjust_method Character string. P-value adjustment method for multiple comparisons.
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'   Default is "BH"
#' @param paired Logical. If TRUE, use paired tests. Default is FALSE
#' @param alpha Numeric. Significance level for determining significance. Default is 0.05
#' @param group_order Character vector or NULL. Specifies the order of groups.
#'   If NULL (default), uses the factor levels or alphabetical order.
#'
#' @return A list containing:
#'   \item{test_used}{Character string of the statistical test performed}
#'   \item{posthoc_used}{Character string of the post-hoc test used (NA for 2 groups)}
#'   \item{overall_p}{Overall p-value from global test (NA for 2 groups)}
#'   \item{pairwise_results}{Data frame with pairwise comparison results including p-values and significance labels}
#'   \item{n_groups}{Number of groups compared}
#'   \item{sample_sizes}{Data frame with sample sizes per group}
#'   \item{summary_stats}{Data frame with summary statistics (mean, median, sd, etc.) per group}
#'
#' @import dplyr
#' @import rstatix
#' @importFrom stats kruskal.test aov wilcox.test t.test sd
#'
#' @examples
#' \dontrun{
#' # Example data
#' set.seed(123)
#' df <- data.frame(
#'   Group = rep(c("A", "B", "C"), each = 20),
#'   Value = c(rnorm(20, 3, 0.5), rnorm(20, 3.5, 0.5), rnorm(20, 4, 0.5))
#' )
#'
#' # Perform statistical comparison
#' stats <- compare_stats(
#'   data = df,
#'   group_var = "Group",
#'   test_var = "Value"
#' )
#'
#' # View results
#' stats$pairwise_results
#' stats$overall_p
#' }
#'
#' @export
compare_stats <- function(data,
                          group_var,
                          test_var,
                          test_method = NULL,
                          posthoc_test = NULL,
                          p_adjust_method = "BH",
                          paired = FALSE,
                          alpha = 0.05,
                          group_order = NULL) {

  # Validate inputs
  if (!group_var %in% names(data)) {
    stop(paste("Group variable", group_var, "not found in data"))
  }

  if (!test_var %in% names(data)) {
    stop(paste("Test variable", test_var, "not found in data"))
  }

  if (!is.numeric(data[[test_var]])) {
    stop(paste("Test variable", test_var, "must be numeric"))
  }

  # Convert group variable to factor
  data[[group_var]] <- as.factor(data[[group_var]])

  # Apply group order if specified
  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(data[[group_var]]))) {
      stop("group_order contains levels not present in the data")
    }
    data[[group_var]] <- factor(data[[group_var]], levels = group_order)
  }

  # Remove missing values
  test_data <- data %>%
    dplyr::select(dplyr::all_of(c(group_var, test_var))) %>%
    dplyr::filter(!is.na(.data[[test_var]]))

  # Get number of groups
  n_groups <- length(unique(test_data[[group_var]]))

  if (n_groups < 2) {
    stop("Need at least 2 groups for comparison")
  }

  # Calculate sample sizes
  sample_sizes <- test_data %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  min_sample_size <- min(sample_sizes$n)

  # Calculate summary statistics
  summary_stats <- test_data %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data[[test_var]], na.rm = TRUE),
      median = median(.data[[test_var]], na.rm = TRUE),
      sd = sd(.data[[test_var]], na.rm = TRUE),
      min = min(.data[[test_var]], na.rm = TRUE),
      max = max(.data[[test_var]], na.rm = TRUE),
      .groups = "drop"
    )

  # Automatically select test method if not specified
  if (is.null(test_method)) {
    if (n_groups == 2) {
      if (paired) {
        selected_test <- if (min_sample_size >= 30) "paired.t.test" else "paired.wilcox.test"
      } else {
        selected_test <- if (min_sample_size >= 30) "t.test" else "wilcox.test"
      }
    } else {
      if (paired) {
        selected_test <- if (min_sample_size >= 30) "anova" else "kruskal.test"
        warning("Paired tests with >2 groups require repeated measures ANOVA. Using standard tests.")
      } else {
        selected_test <- if (min_sample_size >= 30) "anova" else "kruskal.test"
      }
    }
  } else {
    selected_test <- test_method
  }

  # Automatically select post-hoc test if not specified and more than 2 groups
  if (n_groups > 2 && is.null(posthoc_test)) {
    if (selected_test %in% c("anova", "paired.t.test")) {
      selected_posthoc <- "tukey"
    } else if (selected_test %in% c("kruskal.test", "paired.wilcox.test")) {
      selected_posthoc <- "dunn"
    } else {
      selected_posthoc <- "tukey"
    }
  } else {
    selected_posthoc <- posthoc_test
  }

  # Perform statistical test
  overall_p <- NA

  if (n_groups == 2) {
    # Two group comparison
    if (selected_test %in% c("wilcox.test", "paired.wilcox.test")) {
      stat_result <- test_data %>%
        rstatix::wilcox_test(
          as.formula(paste0(test_var, "~", group_var)),
          p.adjust.method = p_adjust_method,
          paired = paired
        )
    } else {
      stat_result <- test_data %>%
        rstatix::t_test(
          as.formula(paste0(test_var, "~", group_var)),
          p.adjust.method = p_adjust_method,
          paired = paired
        )
    }

    # Add significance labels
    stat_result <- stat_result %>%
      dplyr::mutate(
        p.adj = p,
        p.adj.signif = dplyr::case_when(
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        p.format = dplyr::case_when(
          p < 0.001 ~ "p < 0.001",
          TRUE ~ paste0("p = ", format(round(p, 3), nsmall = 3))
        ),
        p.adj.format = dplyr::case_when(
          p < 0.001 ~ "p < 0.001",
          TRUE ~ paste0("p = ", format(round(p, 3), nsmall = 3))
        ),
        significant = p < alpha
      )

  } else {
    # Multiple groups: perform overall test then post-hoc
    if (selected_test == "anova") {
      overall_test <- aov(as.formula(paste0(test_var, "~", group_var)), data = test_data)
      overall_p <- summary(overall_test)[[1]][1, "Pr(>F)"]
    } else if (selected_test == "kruskal.test") {
      overall_test <- kruskal.test(as.formula(paste0(test_var, "~", group_var)), data = test_data)
      overall_p <- overall_test$p.value
    }

    # Perform post-hoc test
    if (selected_posthoc == "tukey") {
      stat_result <- test_data %>%
        rstatix::tukey_hsd(as.formula(paste0(test_var, "~", group_var)))
    } else if (selected_posthoc == "dunn") {
      stat_result <- test_data %>%
        rstatix::dunn_test(as.formula(paste0(test_var, "~", group_var)),
                           p.adjust.method = p_adjust_method)
    } else if (selected_posthoc == "wilcox") {
      stat_result <- test_data %>%
        rstatix::wilcox_test(
          as.formula(paste0(test_var, "~", group_var)),
          p.adjust.method = p_adjust_method
        )
    } else {
      stat_result <- test_data %>%
        rstatix::t_test(
          as.formula(paste0(test_var, "~", group_var)),
          p.adjust.method = p_adjust_method
        )
    }

    # Add significance labels
    stat_result <- stat_result %>%
      dplyr::mutate(
        p.adj.signif = dplyr::case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01 ~ "**",
          p.adj < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        p.format = dplyr::case_when(
          p < 0.001 ~ "p < 0.001",
          TRUE ~ paste0("p = ", format(round(p, 3), nsmall = 3))
        ),
        p.adj.format = dplyr::case_when(
          p.adj < 0.001 ~ "p < 0.001",
          TRUE ~ paste0("p = ", format(round(p.adj, 3), nsmall = 3))
        ),
        significant = p.adj < alpha
      )
  }

  # Add summary statistics to results, showing only mean or median depending on test method
  if (selected_test %in% c("t.test", "paired.t.test", "anova")) {
    # Parametric → means
    stat_result <- stat_result %>%
      dplyr::left_join(
        summary_stats %>% dplyr::select(group1 = !!rlang::sym(group_var), mean1 = mean),
        by = "group1"
      ) %>%
      dplyr::left_join(
        summary_stats %>% dplyr::select(group2 = !!rlang::sym(group_var), mean2 = mean),
        by = "group2"
      ) %>%
      dplyr::mutate(
        interpretation = dplyr::case_when(
          p.adj.signif == "ns" ~ "No significant difference between groups",
          mean1 < mean2 ~ paste0(group1, " has lower mean value than ", group2),
          mean1 > mean2 ~ paste0(group1, " has higher mean value than ", group2),
          TRUE ~ ""
        )
      )
  } else {
    # Nonparametric → medians
    stat_result <- stat_result %>%
      dplyr::left_join(
        summary_stats %>% dplyr::select(group1 = !!rlang::sym(group_var), median1 = median),
        by = "group1"
      ) %>%
      dplyr::left_join(
        summary_stats %>% dplyr::select(group2 = !!rlang::sym(group_var), median2 = median),
        by = "group2"
      ) %>%
      dplyr::mutate(
        interpretation = dplyr::case_when(
          p.adj.signif == "ns" ~ "No significant difference between groups",
          median1 < median2 ~ paste0(group1, " has lower median value than ", group2),
          median1 > median2 ~ paste0(group1, " has higher median value than ", group2),
          TRUE ~ ""
        )
      )
  }

  # Return results
  return(list(
    test_used = selected_test,
    posthoc_used = if (n_groups > 2) selected_posthoc else NA,
    overall_p = overall_p,
    pairwise_results = stat_result,
    n_groups = n_groups,
    sample_sizes = sample_sizes,
    summary_stats = summary_stats
  ))

}
