#' Plot PCA Scores Plot
#'
#' @description
#' Creates a scores plot showing sample positions in principal component space.
#' Supports discrete group coloring, continuous variable coloring, confidence
#' ellipses, and outlier detection.
#'
#' @param pca_result List. Output from `run_pca()`.
#' @param pc Character vector of length 2. Principal components to plot, e.g., c("PC1", "PC2").
#'   Can also be numeric, e.g., c(1, 2). Default: c(1, 2).
#' @param color_by Character. Name of column in `metadata` to use for coloring points.
#'   If the column contains categorical data, points are colored by discrete groups with optional
#'   confidence ellipses. If the column contains numeric data, points are colored using a
#'   continuous gradient scale. Required parameter with no default.
#' @param points_from Character. Name of column in `metadata` to use for labeling points.
#'   This determines which values are used for sample labels in the plot and for outlier
#'   labeling when enabled. If the specified column is not found in metadata, sequential
#'   labels will be used instead. Required parameter with no default.
#' @param title Character. Main plot title. Default: auto-generated based on PCs.
#' @param subtitle Character. Plot subtitle. Default: NULL (no subtitle).
#' @param position Character. Position of title and subtitle. Options: "left", "center", "right".
#'   Default: "center".
#' @param arrange_levels Character vector. Custom order for group levels in legend (only applies
#'   to categorical `color_by`). Default: NULL (alphabetical order).
#' @param ellipse Logical. Draw confidence ellipses around groups. Only applies when `color_by`
#'   is categorical. Requires at least 3 samples per group. Default: TRUE.
#' @param ellipse_type Character. Type of ellipse: "t" (t-distribution, default) or
#'   "norm" (normal distribution). Default: "t".
#' @param ellipse_level Numeric. Confidence level for ellipses (0-1). Default: 0.95 (95% confidence).
#' @param legend Character or NULL. Custom legend title. Default: NULL (uses `color_by` column name).
#' @param show_outliers Logical. Label samples outside the confidence ellipse defined by
#'   `ellipse_level`. Only works when ellipses are drawn (categorical coloring). Default: FALSE.
#' @param label_outliers Character. Which groups to check for outliers. Options: "all"
#'   (default) or specific group names from the `color_by` variable. Can be a vector.
#' @param point_size Numeric. Size of points. Default: 3.
#' @param point_alpha Numeric. Transparency of points (0-1). Default: 0.8.
#' @param theme_base_size Numeric. Base font size for theme. Default: 11.
#' @param verbose Logical. Print messages. Default: TRUE.
#'
#' @return A ggplot2 object.
#'
#' @details
#' **Scores Plot Interpretation:**
#' 
#' - **Clustering**: Samples that cluster together have similar metabolic profiles
#' - **Separation**: Distance between groups indicates metabolic differences
#' - **Outliers**: Samples far from their group may represent:
#'   - Biological outliers (interesting phenotypes)
#'   - Technical issues (sample prep, analysis)
#'   - Mislabeled samples
#' 
#' **Coloring Options:**
#' 
#' The `color_by` parameter determines how points are colored:
#' 
#' 1. **Categorical column**: Colors by discrete groups (e.g., Treatment, Genotype)
#'    - Points are both colored and shaped by this variable
#'    - Confidence ellipses can be drawn with `ellipse = TRUE`
#'    - Outliers can be detected and labeled
#'    
#' 2. **Numeric column**: Gradient coloring by continuous variable (e.g., Time, Dose, Age)
#'    - Uses viridis color scale
#'    - Points are shaped by the same variable (binned into quartiles)
#'    - No ellipses or outlier detection available
#' 
#' **Point Labeling:**
#' 
#' The `points_from` parameter controls how samples are labeled:
#' 
#' - If the specified column exists in metadata, those values are used as labels
#' - If the column doesn't exist, sequential labels (Sample_1, Sample_2, ...) are used
#' - Labels are shown when `show_outliers = TRUE` for outlier samples
#' - Common choices: "Sample", "SampleID", "SubjectID"
#' 
#' **Outlier Detection:**
#' 
#' Outliers are identified using Mahalanobis distance from the group centroid.
#' A sample is flagged as an outlier if it falls outside the confidence
#' ellipse defined by `ellipse_level`. The default (0.95) identifies samples
#' outside the 95% confidence region. This exactly mirrors `ggplot2::stat_ellipse()` behavior.
#' 
#' **When Outlier Detection Works:**
#' - Only with categorical `color_by`
#' - When `ellipse = TRUE`
#' - When groups have at least 3 samples
#' 
#' **Sample Exclusion:**
#' 
#' To exclude samples (e.g., QC samples) from the plot, exclude them during PCA
#' computation using `run_pca(exclude = ...)` rather than filtering in the plot.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). Multivariate Analysis.
#' London: Academic Press.
#' 
#' Brereton, R.G., & Lloyd, G.R. (2014). Partial least squares discriminant analysis:
#' taking the magic away. Journal of Chemometrics, 28(4), 213-225. \doi{10.1002/cem.2609}
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse scale_color_manual scale_shape_manual labs theme_minimal theme element_text
#' @importFrom viridis scale_color_viridis
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats mahalanobis cov qt qnorm
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 50, mean = 100, sd = 20),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' rownames(x) <- paste0("Sample", 1:100)
#' 
#' metadata <- data.frame(
#'   Sample = paste0("Sample", 1:100),
#'   Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
#'   Batch = rep(1:4, each = 25),
#'   Time = rep(c(0, 6, 12, 24), 25)
#' )
#' 
#' # Run PCA (excluding QC samples)
#' scaled_result <- run_scale(x, method = "auto")
#' pca_result <- run_pca(scaled_result$data, metadata, 
#'                       group = "Group",
#'                       exclude = "QC")
#' 
#' # Basic scores plot (colored by Group, labeled by Sample)
#' plot_score(pca_result, 
#'            color_by = "Group",
#'            points_from = "Sample")
#' 
#' # PC2 vs PC3 with custom title
#' plot_score(pca_result, 
#'            pc = c(2, 3),
#'            color_by = "Group",
#'            points_from = "Sample",
#'            title = "PCA: PC2 vs PC3",
#'            position = "left")
#' 
#' # Color by batch
#' plot_score(pca_result,
#'            color_by = "Batch",
#'            points_from = "Sample")
#' 
#' # Continuous coloring by time with custom legend
#' plot_score(pca_result,
#'            color_by = "Time",
#'            points_from = "Sample",
#'            legend = "Collection Time (h)")
#' 
#' # Show outliers in specific groups
#' plot_score(pca_result,
#'            color_by = "Group",
#'            points_from = "Sample",
#'            show_outliers = TRUE,
#'            label_outliers = c("Control", "Treatment"),
#'            ellipse_level = 0.95)
#' 
#' # Stricter outlier detection (99% confidence)
#' plot_score(pca_result,
#'            color_by = "Group",
#'            points_from = "Sample",
#'            show_outliers = TRUE,
#'            ellipse_level = 0.99)
#' }
plot_score <- function(
    pca_result,
    pc = c(1, 2),
    color_by,
    points_from,
    title = NULL,
    subtitle = NULL,
    position = "center",
    arrange_levels = NULL,
    ellipse = TRUE,
    ellipse_type = "t",
    ellipse_level = 0.95,
    legend = NULL,
    show_outliers = FALSE,
    label_outliers = "all",
    point_size = 3,
    point_alpha = 0.8,
    theme_base_size = 11,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  if (!inherits(pca_result, "run_pca")) {
    stop("'pca_result' must be output from run_pca() function.\n",
         "Solution: Run run_pca() first on your data before plotting.")
  }
  
  # Check required parameters
  if (missing(color_by)) {
    stop("'color_by' is required. Please specify a column name from metadata.\n",
         sprintf("Available columns: %s\n", paste(colnames(pca_result$metadata), collapse = ", ")),
         "Solution: Add color_by='ColumnName' to your function call.")
  }
  
  if (missing(points_from)) {
    stop("'points_from' is required. Please specify a column name from metadata for point labels.\n",
         sprintf("Available columns: %s\n", paste(colnames(pca_result$metadata), collapse = ", ")),
         "Solution: Add points_from='ColumnName' (e.g., points_from='Sample') to your function call.")
  }
  
  if (!is.character(color_by) || length(color_by) != 1) {
    stop("'color_by' must be a single character string.\n",
         "Solution: Use color_by='ColumnName' with a valid column from metadata.")
  }
  
  if (!is.character(points_from) || length(points_from) != 1) {
    stop("'points_from' must be a single character string.\n",
         "Solution: Use points_from='ColumnName' with a valid column from metadata.")
  }
  
  position <- tolower(position)
  ellipse_type <- tolower(ellipse_type)
  
  if (!position %in% c("left", "center", "right")) {
    stop("'position' must be 'left', 'center', or 'right'.\n",
         "Solution: Use position='center' (default), position='left', or position='right'.")
  }
  if (!ellipse_type %in% c("t", "norm")) {
    stop("'ellipse_type' must be 't' or 'norm'.\n",
         "Solution: Use ellipse_type='t' (default, t-distribution) or ellipse_type='norm' (normal).")
  }
  if (!is.numeric(ellipse_level) || length(ellipse_level) != 1 || 
      ellipse_level <= 0 || ellipse_level >= 1) {
    stop("'ellipse_level' must be a numeric value between 0 and 1.\n",
         "Solution: Use ellipse_level=0.95 (default, 95% confidence) or another value like 0.99.")
  }
  if (!is.numeric(theme_base_size) || length(theme_base_size) != 1 || theme_base_size <= 0) {
    stop("'theme_base_size' must be a positive numeric value.\n",
         "Solution: Use theme_base_size=11 (default) or another positive number.")
  }
  
  # Parse PC specification
  if (is.numeric(pc)) {
    pc_x <- pc[1]
    pc_y <- pc[2]
  } else if (is.character(pc)) {
    pc_x <- as.numeric(gsub("PC", "", pc[1]))
    pc_y <- as.numeric(gsub("PC", "", pc[2]))
  } else {
    stop("'pc' must be either numeric (e.g., c(1,2)) or character (e.g., c('PC1','PC2')).\n",
         "Solution: Use pc=c(1,2) or pc=c('PC1','PC2').")
  }
  
  if (pc_x < 1 || pc_x > pca_result$n_pcs || pc_y < 1 || pc_y > pca_result$n_pcs) {
    stop(sprintf("Requested PCs out of range. Available: PC1 to PC%d.\n", pca_result$n_pcs),
         sprintf("Solution: Choose PC values between 1 and %d.", pca_result$n_pcs))
  }
  if (pc_x == pc_y) {
    stop("Cannot plot the same PC against itself.\n",
         "Solution: Choose different PC values, e.g., pc=c(1,2) instead of pc=c(1,1).")
  }
  
  # Check color_by column
  if (!color_by %in% colnames(pca_result$metadata)) {
    stop(sprintf("Column '%s' specified in color_by not found in metadata.\n", color_by),
         sprintf("Available columns: %s\n", paste(colnames(pca_result$metadata), collapse = ", ")),
         "Solution: Use a valid column name from metadata.")
  }
  
  # Check for ggrepel when show_outliers is TRUE
  if (show_outliers && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for show_outliers = TRUE.\n",
         "Solution: Install it with: install.packages('ggrepel')")
  }
  
  # ============================================================================
  # DATA PREPARATION
  # ============================================================================
  
  metadata <- pca_result$metadata
  scores <- pca_result$scores
  
  # Check if color column is numeric
  is_numeric_color <- is.numeric(metadata[[color_by]])
  
  if (is_numeric_color) {
    msg(sprintf("Using continuous coloring by '%s'", color_by))
  } else {
    msg(sprintf("Using discrete coloring by '%s'", color_by))
  }
  
  # Create plot data
  pc_x_col <- paste0("PC", pc_x)
  pc_y_col <- paste0("PC", pc_y)
  
  # Handle point labels from points_from
  if (points_from %in% colnames(metadata)) {
    sample_labels <- as.character(metadata[[points_from]])
    msg(sprintf("Using labels from column '%s'", points_from))
  } else {
    warning(sprintf("Column '%s' specified in points_from not found in metadata.\n", points_from),
            sprintf("Available columns: %s\n", paste(colnames(metadata), collapse = ", ")),
            "Solution: Using sequential labels instead (Sample_1, Sample_2, ...).")
    sample_labels <- paste0("Sample_", seq_len(nrow(scores)))
  }
  
  plot_data <- data.frame(
    Label = sample_labels,
    ColorGroup = metadata[[color_by]],
    stringsAsFactors = FALSE
  )
  
  plot_data[[pc_x_col]] <- scores[, pc_x]
  plot_data[[pc_y_col]] <- scores[, pc_y]
  
  # Handle group arrangement (only for categorical)
  if (!is_numeric_color) {
    if (!is.null(arrange_levels)) {
      all_groups <- unique(plot_data$ColorGroup)
      missing_grps <- arrange_levels[!arrange_levels %in% all_groups]
      if (length(missing_grps) > 0) {
        warning(sprintf("The following groups in arrange_levels were not found: %s\n",
                        paste(missing_grps, collapse = ", ")),
                sprintf("Available groups: %s\n", paste(all_groups, collapse = ", ")),
                "Solution: Check spelling or remove non-existent groups from arrange_levels.")
      }
      valid_levels <- arrange_levels[arrange_levels %in% all_groups]
      remaining <- all_groups[!all_groups %in% valid_levels]
      plot_data$ColorGroup <- factor(plot_data$ColorGroup, 
                                     levels = c(valid_levels, sort(remaining)))
    } else {
      plot_data$ColorGroup <- factor(plot_data$ColorGroup)
    }
  }
  
  # ============================================================================
  # OUTLIER COMPUTATION
  # ============================================================================
  
  if (show_outliers) {
    if (is_numeric_color) {
      warning("Outlier detection not available with continuous coloring (color_by is numeric).\n",
              "Solution: Use a categorical variable for color_by if you need outlier detection.")
      show_outliers <- FALSE
    } else {
      msg(sprintf("Computing outliers (ellipse level: %.1f%%)...", ellipse_level * 100))
      
      plot_data <- .compute_outliers_discrete(plot_data, pc_x_col, pc_y_col,
                                              ellipse_type, ellipse_level)
      
      n_outliers <- sum(plot_data$is_outlier, na.rm = TRUE)
      msg(sprintf("Outliers detected: %d (%.1f%%)", 
                  n_outliers, 100 * n_outliers / nrow(plot_data)))
      
      # Resolve which groups to label
      resolved_outlier_groups <- .resolve_outlier_groups_simple(
        label_outliers, plot_data$ColorGroup
      )
      
      if (length(resolved_outlier_groups) == 0) {
        warning("No valid groups resolved for outlier labeling. Skipping outlier labels.\n",
                sprintf("Available groups: %s\n", paste(levels(plot_data$ColorGroup), collapse = ", ")),
                "Solution: Check that label_outliers contains valid group names or use 'all'.")
        show_outliers <- FALSE
      }
    }
  }
  
  # ============================================================================
  # CREATE PLOT
  # ============================================================================
  
  # Create axis labels with variance explained
  x_label <- sprintf("PC%d (%.1f%%)", pc_x, pca_result$variance_explained[pc_x])
  y_label <- sprintf("PC%d (%.1f%%)", pc_y, pca_result$variance_explained[pc_y])
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- sprintf("PCA Scores Plot (PC%d vs PC%d)", pc_x, pc_y)
  }
  
  # Create base plot
  if (is_numeric_color) {
    # Continuous coloring - use same variable for shape (binned)
    plot_data$ShapeGroup <- cut(plot_data$ColorGroup,
                                breaks = quantile(plot_data$ColorGroup, 
                                                 probs = c(0, 0.25, 0.5, 0.75, 1),
                                                 na.rm = TRUE),
                                labels = c("Q1", "Q2", "Q3", "Q4"),
                                include.lowest = TRUE)
    
    p <- ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = .data[[pc_x_col]],
                                    y = .data[[pc_y_col]],
                                    color = ColorGroup,
                                    shape = ShapeGroup))
  } else {
    # Discrete coloring - use same variable for shape
    p <- ggplot2::ggplot(plot_data, 
                        ggplot2::aes(x = .data[[pc_x_col]], 
                                    y = .data[[pc_y_col]],
                                    color = ColorGroup, 
                                    shape = ColorGroup))
  }
  
  # Add points
  p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha)
  
  # Add ellipses (only for discrete coloring)
  if (ellipse && !is_numeric_color) {
    # Check which groups have enough samples for ellipse
    grp_counts <- table(plot_data$ColorGroup)
    grps_for_ellipse <- names(grp_counts)[grp_counts >= 3]
    
    if (length(grps_for_ellipse) > 0) {
      df_ell <- plot_data[plot_data$ColorGroup %in% grps_for_ellipse, ]
      p <- p + ggplot2::stat_ellipse(
        data = df_ell,
        ggplot2::aes(group = ColorGroup, color = ColorGroup, fill = ColorGroup),
        type = ellipse_type,
        level = ellipse_level,
        geom = "polygon",
        alpha = 0.2,
        show.legend = FALSE
      )
    } else {
      warning("No groups with >= 3 samples. Ellipses cannot be drawn.\n",
              "Solution: Ensure groups have at least 3 samples each, or set ellipse=FALSE.")
    }
  }
  
  # Add color scale
  if (is_numeric_color) {
    legend_title <- if (!is.null(legend)) legend else color_by
    p <- p + viridis::scale_color_viridis(
      option = "plasma",
      name = legend_title
    )
  }
  
  # Add outlier labels
  if (show_outliers && !is_numeric_color) {
    df_out <- plot_data[
      !is.na(plot_data$is_outlier) & 
        plot_data$is_outlier &
        as.character(plot_data$ColorGroup) %in% resolved_outlier_groups,
    ]
    
    if (nrow(df_out) > 0) {
      p <- p + ggrepel::geom_label_repel(
        data = df_out,
        ggplot2::aes(label = Label, color = ColorGroup),
        fill = "white",
        size = theme_base_size / 3.5,  # Scale with base size
        label.padding = ggplot2::unit(0.15, "lines"),
        box.padding = ggplot2::unit(0.5, "lines"),
        point.padding = ggplot2::unit(0.3, "lines"),
        segment.color = "grey40",
        min.segment.length = 0,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    }
  }
  
  # Determine title justification
  hjust_val <- switch(position,
                     "left" = 0,
                     "center" = 0.5,
                     "right" = 1,
                     0.5)
  
  # Apply theme
  p <- p +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = hjust_val, size = theme_base_size + 3, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = hjust_val, size = theme_base_size + 1),
      axis.title = ggplot2::element_text(size = theme_base_size + 1),
      axis.text = ggplot2::element_text(size = theme_base_size - 1),
      legend.title = ggplot2::element_text(size = theme_base_size),
      legend.text = ggplot2::element_text(size = theme_base_size - 1)
    )
  
  msg("Scores plot created successfully")
  return(p)
}


# ==============================================================================
# HELPER FUNCTIONS FOR OUTLIER DETECTION
# ==============================================================================

.compute_outliers_discrete <- function(df, pc1_col, pc2_col, ellipse_type, level) {
  mahal_d <- rep(NA_real_, nrow(df))
  threshold <- rep(NA_real_, nrow(df))
  is_out <- rep(NA, nrow(df))
  
  for (grp in levels(df$ColorGroup)) {
    idx <- which(df$ColorGroup == grp)
    n <- length(idx)
    if (n < 3) next
    
    pts <- as.matrix(df[idx, c(pc1_col, pc2_col), drop = FALSE])
    cov_mat <- tryCatch(stats::cov(pts), error = function(e) NULL)
    if (is.null(cov_mat) || det(cov_mat) <= .Machine$double.eps) next
    
    # Compute threshold based on ellipse type
    if (ellipse_type == "t") {
      crit <- stats::qt(level/2 + 0.5, df = n - 1)^2
    } else {
      crit <- stats::qnorm(level/2 + 0.5)^2
    }
    
    d2 <- tryCatch(
      stats::mahalanobis(pts, colMeans(pts), cov_mat),
      error = function(e) rep(NA_real_, n)
    )
    
    mahal_d[idx] <- d2
    threshold[idx] <- crit
    is_out[idx] <- !is.na(d2) & (d2 > crit)
  }
  
  df$mahal_dist <- mahal_d
  df$ellipse_threshold <- threshold
  df$is_outlier <- is_out
  df
}

.resolve_outlier_groups_simple <- function(label_outliers, all_groups) {
  if ("all" %in% label_outliers) {
    return(levels(all_groups))
  }
  
  valid_groups <- label_outliers[label_outliers %in% levels(all_groups)]
  
  if (length(valid_groups) < length(label_outliers)) {
    invalid <- label_outliers[!label_outliers %in% levels(all_groups)]
    warning(sprintf("The following groups in label_outliers were not found: %s\n",
                    paste(invalid, collapse = ", ")),
            sprintf("Available groups: %s", paste(levels(all_groups), collapse = ", ")))
  }
  
  unique(valid_groups)
}