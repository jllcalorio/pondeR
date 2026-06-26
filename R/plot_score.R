#' @title Scores Plot for Multivariate Ordination Results
#'
#' @description
#' Produces a publication-ready scores plot from the output of multivariate
#' ordination and dimensionality-reduction methods, including Principal
#' Component Analysis (\code{run_pca}), Partial Least Squares and PLS-DA
#' (\code{run_pls}), and Principal Coordinate Analysis (\code{run_pcoa}).
#' Supports discrete group coloring with confidence ellipses, continuous
#' gradient coloring, outlier detection and labeling, and a range of
#' publication themes.
#'
#' @param res List. Output from \code{run_pca()}, \code{run_pls()}, or
#'   \code{run_pcoa()}. The appropriate plot method is selected automatically
#'   via S3 dispatch based on the class of \code{res}.
#' @param pc Character vector of length 2. Principal components to plot, e.g., `c("PC1", "PC2")`.
#'   Can also be numeric, e.g., `c(1, 2)`. Default: `c(1, 2)`.
#' @param color_by Character. Name of column in `metadata` to use for coloring points.
#'   If the column contains categorical data, points are colored by discrete groups with optional
#'   confidence ellipses. If the column contains numeric data, points are colored using a
#'   continuous gradient scale. Required parameter with no default.
#' @param points_from Character. Name of column in `metadata` to use for labeling points.
#'   This determines which values are used for sample labels in the plot and for outlier
#'   labeling when enabled. If the specified column is not found in metadata, sequential
#'   labels will be used instead. Required parameter with no default.
#' @param title Character. Main plot title. Default: auto-generated based on PCs.
#' @param subtitle Character. Plot subtitle. Default: `NULL` (no subtitle).
#' @param caption Character. Plot caption (bottom-right). Default: `NULL` (no caption).
#' @param position Character. Horizontal alignment of title, subtitle, and caption.
#'   Options: `"left"`, `"center"`, `"right"`. Default: `"center"`.
#' @param arrange_levels Character vector. Custom order for group levels in legend (only applies
#'   to categorical `color_by`). Default: `NULL` (alphabetical order).
#' @param ellipse Logical. Draw confidence ellipses around groups. Only applies when `color_by`
#'   is categorical. Requires at least 3 samples per group. Default: `TRUE`.
#' @param ellipse_type Character. Type of ellipse: `"t"` (t-distribution, default) or
#'   `"norm"` (normal distribution). Default: `"t"`.
#' @param ellipse_level Numeric. Confidence level for ellipses (0-1). Default: `0.95`.
#' @param legend Character or `NULL`. Custom legend title. Default: `NULL` (uses `color_by`
#'   column name).
#' @param legend_position Character. Position of the legend. Options: `"bottom"`, `"top"`,
#'   `"left"`, `"right"`, `"none"`. Default: `"bottom"`.
#' @param show_outliers Logical. Label samples outside the confidence ellipse defined by
#'   `ellipse_level`. Only works when ellipses are drawn (categorical coloring). Default: `FALSE`.
#' @param label_outliers Character. Which groups to check for outliers. Use `"all"` (default)
#'   or a character vector of specific group names from the `color_by` variable.
#' @param colors Character vector. Custom color palette for discrete groups. If `NULL`
#'   (default), the Okabe-Ito colorblind-friendly palette is used. Ignored when `color_by`
#'   is numeric (continuous scale is always used for numeric variables).
#' @param point_size Numeric. Size of points. Default: `10`.
#' @param point_alpha Numeric. Transparency of points (0-1). Default: `0.8`.
#' @param theme Character. ggplot2 theme to apply. Options: `"nature"`, `"minimal"`,
#'   `"classic"`, `"bw"`, `"light"`, `"dark"`. Default: `"nature"` (a clean, publication-ready
#'   style based on `theme_bw()`).
#' @param base_size Numeric. Base font size for the theme (pts). Default: `30`.
#' @param font_family Character. Font family for all text elements. Default: `"sans"`.
#' @param axis_title_size Numeric. Font size for axis titles. Default: `base_size + 1`.
#' @param axis_text_size Numeric. Font size for axis tick labels. Default: `base_size - 1`.
#' @param plot_title_size Numeric. Font size for the plot title. Default: `base_size + 3`.
#' @param legend_title_size Numeric. Font size for the legend title. Default: `base_size`.
#' @param legend_text_size Numeric. Font size for legend text. Default: `base_size - 1`.
#' @param zoom Numeric. Scaling factor applied uniformly to all text and point size elements.
#'   Values > 1 enlarge, values < 1 shrink. Default: `1`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `ggplot2` object.
#'
#' @details
#' **Scores Plot Interpretation:**
#'
#' - **Clustering**: Samples that cluster together have similar metabolic profiles.
#' - **Separation**: Distance between groups indicates metabolic differences.
#' - **Outliers**: Samples far from their group may represent biological outliers,
#'   technical issues, or mislabeled samples.
#'
#' **Coloring Options:**
#'
#' The `color_by` parameter determines how points are colored:
#'
#' 1. **Categorical column**: Colors by discrete groups (e.g., Treatment, Genotype).
#'    - Points are both colored and shaped by this variable.
#'    - Confidence ellipses can be drawn with `ellipse = TRUE`.
#'    - Outliers can be detected and labeled.
#'    - Colors default to the Okabe-Ito colorblind-safe palette; override with `colors`.
#'
#' 2. **Numeric column**: Gradient coloring by continuous variable (e.g., Time, Dose, Age).
#'    - Uses the viridis "plasma" color scale.
#'    - Points are shaped by the same variable (binned into quartiles).
#'    - No ellipses or outlier detection available.
#'
#' **Color Defaults:**
#'
#' By default, discrete groups use the Okabe-Ito palette, which is perceptually
#' distinct and safe for the most common forms of color vision deficiency (deuteranopia,
#' protanopia, tritanopia). The palette supports up to 8 groups. For more than 8 groups,
#' supply a custom `colors` vector.
#'
#' **Point Labeling:**
#'
#' The `points_from` parameter controls how samples are labeled:
#'
#' - If the specified column exists in metadata, those values are used as labels.
#' - If the column doesn't exist, sequential labels (`Sample_1`, `Sample_2`, ...) are used.
#' - Labels are shown when `show_outliers = TRUE` for outlier samples only.
#' - Common choices: `"Sample"`, `"SampleID"`, `"SubjectID"`.
#'
#' **Outlier Detection:**
#'
#' Outliers are identified using Mahalanobis distance from the group centroid.
#' A sample is flagged as an outlier if it falls outside the confidence region defined
#' by `ellipse_level`. The default (`0.95`) identifies samples outside the 95% confidence
#'   region. This exactly mirrors `ggplot2::stat_ellipse()` behavior.
#'
#' Outlier detection requires:
#' - Categorical `color_by`
#' - `ellipse = TRUE`
#' - At least 3 samples per group
#'
#' **Sample Exclusion:**
#'
#' To exclude samples (e.g., QC samples) from the plot, exclude them during PCA
#' computation using `run_pca(exclude = ...)` rather than filtering at the plot stage.
#'
#' **Themes:**
#'
#' The `"nature"` theme (default) is a clean, publication-ready style with a white
#' background, minimal gridlines, and no top/right panel border - suitable for
#' journal figures. Other options map directly to their ggplot2 equivalents.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). *Multivariate Analysis*.
#' London: Academic Press.
#'
#' Brereton, R.G., & Lloyd, G.R. (2014). Partial least squares discriminant analysis:
#' taking the magic away. *Journal of Chemometrics*, 28(4), 213-225.
#' \doi{10.1002/cem.2609}
#'
#' Okabe, M., & Ito, K. (2002). *Color Universal Design (CUD): How to make figures
#' and presentations that are friendly to colorblind people*.
#' \url{https://jfly.uni-koeln.de/color/}
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse scale_color_manual
#' @importFrom ggplot2 scale_shape_manual scale_fill_manual labs theme_bw theme_minimal theme_classic
#' @importFrom ggplot2 theme_light theme_dark theme element_text element_blank element_line element_rect unit
#' @importFrom viridis scale_color_viridis
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats mahalanobis cov qt qnorm
#'
#' @seealso \code{\link[dimsprepr]{run_DIpreprocess}}, \code{\link{run_pca}}, \code{\link{run_pls}}
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
#'   Group  = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
#'   Batch  = rep(1:4, each = 25),
#'   Time   = rep(c(0, 6, 12, 24), 25)
#' )
#'
#' # Run PCA (excluding QC samples)
#' scaled_result <- run_scale(x, method = "auto")
#' res    <- run_pca(scaled_result$data, metadata,
#'                          group   = "Group",
#'                          exclude = "QC")
#'
#' # Basic scores plot - nature theme, Okabe-Ito colors (defaults)
#' plot_score(res,
#'            color_by    = "Group",
#'            points_from = "Sample")
#'
#' # PC2 vs PC3 with left-aligned title
#' plot_score(res,
#'            pc          = c(2, 3),
#'            color_by    = "Group",
#'            points_from = "Sample",
#'            title       = "PCA: PC2 vs PC3",
#'            position    = "left")
#'
#' # Color by batch with custom palette
#' plot_score(res,
#'            color_by    = "Batch",
#'            points_from = "Sample",
#'            colors      = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))
#'
#' # Continuous coloring by time
#' plot_score(res,
#'            color_by    = "Time",
#'            points_from = "Sample",
#'            legend      = "Collection Time (h)")
#'
#' # Show outliers in specific groups, classic theme, larger zoom
#' plot_score(res,
#'            color_by       = "Group",
#'            points_from    = "Sample",
#'            show_outliers  = TRUE,
#'            label_outliers = c("Control", "Treatment"),
#'            ellipse_level  = 0.95,
#'            theme          = "classic",
#'            zoom           = 1.2)
#'
#' # Stricter outlier detection with serif font
#' plot_score(res,
#'            color_by      = "Group",
#'            points_from   = "Sample",
#'            show_outliers = TRUE,
#'            ellipse_level = 0.99,
#'            font_family   = "serif",
#'            base_size     = 13)
#' }
#' 
plot_score <- function(res, ...) UseMethod("plot_score")
#'
#' @rdname plot_score
#' @export
plot_score.run_pca <- function(
    res,
    pc               = c(1, 2),
    color_by,
    points_from,
    title            = NULL,
    subtitle         = NULL,
    caption          = NULL,
    position         = "center",
    arrange_levels   = NULL,
    ellipse          = TRUE,
    ellipse_type     = "t",
    ellipse_level    = 0.95,
    legend           = NULL,
    legend_position  = "bottom",
    show_outliers    = FALSE,
    label_outliers   = "all",
    colors           = NULL,
    point_size       = 10,
    point_alpha      = 0.8,
    theme            = "nature",
    base_size        = 30,
    font_family      = "sans",
    axis_title_size  = NULL,
    axis_text_size   = NULL,
    plot_title_size  = NULL,
    legend_title_size = NULL,
    legend_text_size = NULL,
    zoom             = 1,
    verbose          = TRUE,
    ...
) {

  msg <- function(...) if (verbose) message(...)

  # Okabe-Ito colorblind-safe palette (8 colors), ordered for maximum perceptual
  # contrast between adjacent groups: blue, vermillion, green, sky-blue, orange,
  # purple, yellow, black. Yellow is pushed to position 7 to avoid similarity
  # with orange when <= 6 groups are used.
  .okabe_ito <- c(
    "#0072B2", "#D55E00", "#009E73", "#56B4E9",
    "#E69F00", "#CC79A7", "#F0E442", "#000000"
  )

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  if (missing(color_by)) {
    stop(
      "'color_by' is required. Please specify a column name from metadata.\n",
      sprintf("Available columns: %s\n", paste(colnames(res$metadata), collapse = ", ")),
      "Solution: Add color_by = 'ColumnName' to your function call."
    )
  }

  if (missing(points_from)) {
    stop(
      "'points_from' is required. Please specify a column name from metadata for point labels.\n",
      sprintf("Available columns: %s\n", paste(colnames(res$metadata), collapse = ", ")),
      "Solution: Add points_from = 'ColumnName' (e.g., points_from = 'Sample') to your function call."
    )
  }

  if (!is.character(color_by) || length(color_by) != 1L) {
    stop(
      "'color_by' must be a single character string.\n",
      "Solution: Use color_by = 'ColumnName' with a valid column from metadata."
    )
  }

  if (!is.character(points_from) || length(points_from) != 1L) {
    stop(
      "'points_from' must be a single character string.\n",
      "Solution: Use points_from = 'ColumnName' with a valid column from metadata."
    )
  }

  position        <- tolower(position)
  ellipse_type    <- tolower(ellipse_type)
  theme_choice    <- tolower(theme)
  legend_position <- tolower(legend_position)

  if (!position %in% c("left", "center", "right")) {
    stop(
      "'position' must be 'left', 'center', or 'right'.\n",
      "Solution: Use position = 'center' (default), 'left', or 'right'."
    )
  }

  if (!ellipse_type %in% c("t", "norm")) {
    stop(
      "'ellipse_type' must be 't' or 'norm'.\n",
      "Solution: Use ellipse_type = 't' (default) or ellipse_type = 'norm'."
    )
  }

  if (!is.numeric(ellipse_level) || length(ellipse_level) != 1L ||
      ellipse_level <= 0 || ellipse_level >= 1) {
    stop(
      "'ellipse_level' must be a numeric value strictly between 0 and 1.\n",
      "Solution: Use ellipse_level = 0.95 (default) or another value such as 0.99."
    )
  }

  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  if (!theme_choice %in% valid_themes) {
    stop(
      sprintf("'theme' must be one of: %s.\n", paste(valid_themes, collapse = ", ")),
      "Solution: Use theme = 'nature' (default) or another listed option."
    )
  }

  valid_legend_pos <- c("bottom", "top", "left", "right", "none")
  if (!legend_position %in% valid_legend_pos) {
    stop(
      sprintf("'legend_position' must be one of: %s.\n", paste(valid_legend_pos, collapse = ", ")),
      "Solution: Use legend_position = 'bottom' (default) or another listed option."
    )
  }

  if (!is.numeric(zoom) || length(zoom) != 1L || zoom <= 0) {
    stop(
      "'zoom' must be a single positive numeric value.\n",
      "Solution: Use zoom = 1 (default), zoom = 1.5 to enlarge, or zoom = 0.8 to shrink."
    )
  }

  if (!is.numeric(base_size) || length(base_size) != 1L || base_size <= 0) {
    stop(
      "'base_size' must be a positive numeric value.\n",
      "Solution: Use base_size = 11 (default) or another positive number."
    )
  }

  if (!is.character(font_family) || length(font_family) != 1L) {
    stop(
      "'font_family' must be a single character string.\n",
      "Solution: Use font_family = 'sans' (default), 'serif', or 'mono'."
    )
  }

  # Resolve size defaults (allow per-element overrides, fall back to base_size offsets)
  axis_title_size   <- if (is.null(axis_title_size))   base_size + 1 else axis_title_size
  axis_text_size    <- if (is.null(axis_text_size))    base_size - 1 else axis_text_size
  plot_title_size   <- if (is.null(plot_title_size))   base_size + 3 else plot_title_size
  legend_title_size <- if (is.null(legend_title_size)) base_size     else legend_title_size
  legend_text_size  <- if (is.null(legend_text_size))  base_size - 1 else legend_text_size

  # Apply zoom uniformly
  .z <- function(x) x * zoom

  # Parse PC specification
  if (is.numeric(pc)) {
    pc_x <- as.integer(pc[1L])
    pc_y <- as.integer(pc[2L])
  } else if (is.character(pc)) {
    pc_x <- as.integer(gsub("PC", "", pc[1L], ignore.case = TRUE))
    pc_y <- as.integer(gsub("PC", "", pc[2L], ignore.case = TRUE))
  } else {
    stop(
      "'pc' must be numeric (e.g., c(1, 2)) or character (e.g., c('PC1', 'PC2')).\n",
      "Solution: Use pc = c(1, 2) or pc = c('PC1', 'PC2')."
    )
  }

  if (anyNA(c(pc_x, pc_y))) {
    stop(
      "Could not parse 'pc' values into valid integers.\n",
      "Solution: Use pc = c(1, 2) or pc = c('PC1', 'PC2')."
    )
  }

  if (pc_x < 1L || pc_x > res$n_pcs || pc_y < 1L || pc_y > res$n_pcs) {
    stop(
      sprintf("Requested PCs are out of range. Available: PC1 to PC%d.\n", res$n_pcs),
      sprintf("Solution: Choose pc values between 1 and %d.", res$n_pcs)
    )
  }

  if (pc_x == pc_y) {
    stop(
      "Cannot plot the same PC against itself.\n",
      "Solution: Choose different PC values, e.g., pc = c(1, 2) instead of pc = c(1, 1)."
    )
  }

  if (!color_by %in% colnames(res$metadata)) {
    stop(
      sprintf("Column '%s' specified in 'color_by' not found in metadata.\n", color_by),
      sprintf("Available columns: %s\n", paste(colnames(res$metadata), collapse = ", ")),
      "Solution: Use a valid column name from metadata."
    )
  }

  if (show_outliers && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package 'ggrepel' is required when show_outliers = TRUE.\n",
      "Solution: Install it with: install.packages('ggrepel')"
    )
  }

  # ============================================================================
  # DATA PREPARATION
  # ============================================================================

  metadata <- res$metadata
  scores   <- res$scores

  is_numeric_color <- is.numeric(metadata[[color_by]])

  if (is_numeric_color) {
    msg(sprintf("Using continuous coloring by '%s'", color_by))
  } else {
    msg(sprintf("Using discrete coloring by '%s'", color_by))
  }

  pc_x_col <- paste0("PC", pc_x)
  pc_y_col <- paste0("PC", pc_y)

  if (points_from %in% colnames(metadata)) {
    sample_labels <- as.character(metadata[[points_from]])
    msg(sprintf("Using labels from column '%s'", points_from))
  } else {
    warning(
      sprintf("Column '%s' specified in 'points_from' not found in metadata.\n", points_from),
      sprintf("Available columns: %s\n", paste(colnames(metadata), collapse = ", ")),
      "Solution: Using sequential labels instead (Sample_1, Sample_2, ...)."
    )
    sample_labels <- paste0("Sample_", seq_len(nrow(scores)))
  }

  plot_data <- data.frame(
    Label      = sample_labels,
    ColorGroup = metadata[[color_by]],
    stringsAsFactors = FALSE
  )
  plot_data[[pc_x_col]] <- scores[, pc_x]
  plot_data[[pc_y_col]] <- scores[, pc_y]

  # Arrange group levels (categorical only)
  if (!is_numeric_color) {
    if (!is.null(arrange_levels)) {
      all_groups   <- unique(plot_data$ColorGroup)
      missing_grps <- arrange_levels[!arrange_levels %in% all_groups]
      if (length(missing_grps) > 0L) {
        warning(
          sprintf("Groups in arrange_levels not found: %s\n", paste(missing_grps, collapse = ", ")),
          sprintf("Available groups: %s\n", paste(all_groups, collapse = ", ")),
          "Solution: Check spelling or remove non-existent groups from arrange_levels."
        )
      }
      valid_levels <- arrange_levels[arrange_levels %in% all_groups]
      remaining    <- sort(all_groups[!all_groups %in% valid_levels])
      plot_data$ColorGroup <- factor(plot_data$ColorGroup, levels = c(valid_levels, remaining))
    } else {
      plot_data$ColorGroup <- factor(plot_data$ColorGroup)
    }

    # Resolve color palette
    n_groups <- nlevels(plot_data$ColorGroup)
    if (is.null(colors)) {
      if (n_groups > length(.okabe_ito)) {
        warning(
          sprintf(
            "More groups (%d) than Okabe-Ito palette colors (%d). Colors will be recycled.\n",
            n_groups, length(.okabe_ito)
          ),
          "Solution: Supply a custom 'colors' vector with at least ", n_groups, " colors."
        )
      }
      pal <- rep_len(.okabe_ito, n_groups)
    } else {
      if (length(colors) < n_groups) {
        warning(
          sprintf(
            "'colors' has fewer entries (%d) than groups (%d). Colors will be recycled.\n",
            length(colors), n_groups
          ),
          "Solution: Supply at least ", n_groups, " colors."
        )
        pal <- rep_len(colors, n_groups)
      } else {
        pal <- colors[seq_len(n_groups)]
      }
    }
    names(pal) <- levels(plot_data$ColorGroup)
  }

  # ============================================================================
  # OUTLIER COMPUTATION
  # ============================================================================

  if (show_outliers) {
    if (is_numeric_color) {
      warning(
        "Outlier detection is not available with continuous coloring (color_by is numeric).\n",
        "Solution: Use a categorical variable for color_by if you need outlier detection."
      )
      show_outliers <- FALSE
    } else {
      msg(sprintf("Computing outliers (ellipse level: %.1f%%)...", ellipse_level * 100))

      plot_data  <- .compute_outliers_discrete(plot_data, pc_x_col, pc_y_col,
                                               ellipse_type, ellipse_level)
      n_outliers <- sum(plot_data$is_outlier, na.rm = TRUE)
      msg(sprintf("Outliers detected: %d (%.1f%%)",
                  n_outliers, 100 * n_outliers / nrow(plot_data)))

      resolved_outlier_groups <- .resolve_outlier_groups_simple(
        label_outliers, plot_data$ColorGroup
      )

      if (length(resolved_outlier_groups) == 0L) {
        warning(
          "No valid groups resolved for outlier labeling. Skipping outlier labels.\n",
          sprintf("Available groups: %s\n", paste(levels(plot_data$ColorGroup), collapse = ", ")),
          "Solution: Check that label_outliers contains valid group names or use 'all'."
        )
        show_outliers <- FALSE
      }
    }
  }

  # ============================================================================
  # THEME CONSTRUCTION
  # ============================================================================

  base_theme <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
      ggplot2::theme(
        panel.grid.major  = ggplot2::element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor  = ggplot2::element_blank(),
        panel.border      = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.6),
        strip.background  = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text        = ggplot2::element_text(face = "bold", family = font_family)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size = base_size, base_family = font_family),
    "light"   = ggplot2::theme_light(base_size = base_size, base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size = base_size, base_family = font_family)
  )

  hjust_val <- switch(position, "left" = 0, "center" = 0.5, "right" = 1, 0.5)

  custom_theme <- ggplot2::theme(
    panel.grid.major   = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    legend.position    = legend_position,
    plot.title         = ggplot2::element_text(
      hjust = hjust_val, size = .z(plot_title_size), face = "bold", family = font_family
    ),
    plot.subtitle      = ggplot2::element_text(
      hjust = hjust_val, size = .z(plot_title_size - 2), family = font_family
    ),
    plot.caption       = ggplot2::element_text(
      hjust = 1, size = .z(base_size - 2), color = "grey50", family = font_family
    ),
    axis.title         = ggplot2::element_text(size = .z(axis_title_size), family = font_family),
    axis.text          = ggplot2::element_text(size = .z(axis_text_size),  family = font_family),
    legend.title       = ggplot2::element_text(size = .z(legend_title_size), family = font_family),
    legend.text        = ggplot2::element_text(size = .z(legend_text_size),  family = font_family)
  )

  # ============================================================================
  # BUILD PLOT
  # ============================================================================

  x_label <- sprintf("PC%d (%.1f%%)", pc_x, res$variance_explained[pc_x])
  y_label <- sprintf("PC%d (%.1f%%)", pc_y, res$variance_explained[pc_y])

  if (is.null(title)) {
    title <- sprintf("PCA Scores Plot (PC%d vs PC%d)", pc_x, pc_y)
  }

  if (is_numeric_color) {
    # Continuous: bin into quartiles for shape aesthetic
    plot_data$ShapeGroup <- .bin_continuous_safe(plot_data$ColorGroup)

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x     = .data[[pc_x_col]],
        y     = .data[[pc_y_col]],
        color = ColorGroup,
        shape = ShapeGroup
      )
    )
  } else {
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x     = .data[[pc_x_col]],
        y     = .data[[pc_y_col]],
        color = ColorGroup,
        shape = ColorGroup
      )
    )
  }

  # Points
  p <- p + ggplot2::geom_point(size = .z(point_size), alpha = point_alpha)

  # Discrete color + shape + fill scales
  if (!is_numeric_color) {
    legend_title <- if (!is.null(legend)) legend else color_by
    p <- p +
      ggplot2::scale_color_manual(values = pal, name = legend_title) +
      ggplot2::scale_shape_manual(
        values = seq_len(nlevels(plot_data$ColorGroup)),
        name   = legend_title
      ) +
      ggplot2::scale_fill_manual(values = pal, name = legend_title, guide = "none")
  }

  # Confidence ellipses (categorical only)
  if (ellipse && !is_numeric_color) {
    grp_counts     <- table(plot_data$ColorGroup)
    grps_for_ellipse <- names(grp_counts)[grp_counts >= 3L]

    if (length(grps_for_ellipse) > 0L) {
      df_ell <- plot_data[plot_data$ColorGroup %in% grps_for_ellipse, ]
      p <- p + ggplot2::stat_ellipse(
        data         = df_ell,
        ggplot2::aes(group = ColorGroup, color = ColorGroup, fill = ColorGroup),
        type         = ellipse_type,
        level        = ellipse_level,
        geom         = "polygon",
        alpha        = 0.15,
        linewidth    = 0.6,
        show.legend  = FALSE
      )
    } else {
      warning(
        "No groups with >= 3 samples found. Ellipses cannot be drawn.\n",
        "Solution: Ensure groups have at least 3 samples each, or set ellipse = FALSE."
      )
    }
  }

  # Continuous color scale (viridis plasma)
  if (is_numeric_color) {
    legend_title <- if (!is.null(legend)) legend else color_by
    p <- p + viridis::scale_color_viridis(option = "plasma", name = legend_title)
  }

  # Outlier labels
  if (show_outliers && !is_numeric_color) {
    df_out <- plot_data[
      !is.na(plot_data$is_outlier) &
        plot_data$is_outlier &
        as.character(plot_data$ColorGroup) %in% resolved_outlier_groups,
    ]

    if (nrow(df_out) > 0L) {
      p <- p + ggrepel::geom_label_repel(
        data              = df_out,
        ggplot2::aes(label = Label, color = ColorGroup),
        fill              = "white",
        family            = font_family,
        size              = .z(base_size / 3.5),
        label.padding     = ggplot2::unit(0.15, "lines"),
        box.padding       = ggplot2::unit(0.5,  "lines"),
        point.padding     = ggplot2::unit(0.3,  "lines"),
        segment.color     = "grey40",
        min.segment.length = 0,
        max.overlaps      = Inf,
        show.legend       = FALSE
      )
    }
  }

  # Labels + theme
  p <- p +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = title,
      subtitle = subtitle,
      caption  = caption
    ) +
    base_theme +
    custom_theme

  msg("Scores plot created successfully.")
  return(p)
}

#' @rdname plot_score
#' @export
plot_score.run_pls <- function(
    res,
    pc               = c(1, 2),
    title            = NULL,
    subtitle         = NULL,
    ...
) {
  # -- Parse PC indices ----------------------------------------------------------
  if (is.numeric(pc)) {
    pc_x <- as.integer(pc[1L]); pc_y <- as.integer(pc[2L])
  } else {
    pc_x <- as.integer(gsub("PC", "", pc[1L], ignore.case = TRUE))
    pc_y <- as.integer(gsub("PC", "", pc[2L], ignore.case = TRUE))
  }

  # -- Method label --------------------------------------------------------------
  method_label <- switch(res$method_used,
                         "oplsda" = "OPLS-DA",
                         "plsda"  = "PLS-DA",
                         "splsda" = "sPLS-DA",
                         toupper(res$method_used))

  # ============================================================================
  # EARLY EXIT: no usable scores
  # Conditions that make the result unplottable:
  #   (a) scores is NULL or has 0 columns
  #   (b) n_pcs < 2  - can't form an X vs Y axes plot
  #   (c) all variance_explained are NA  - model produced no valid component
  # ============================================================================

  scores_unusable <- is.null(res$scores)  ||
    ncol(res$scores) < 2L               ||
    all(is.na(res$variance_explained))

  if (scores_unusable) {
    warning(
      sprintf(
        paste0(
          "[plot_score] %s model for this dataset produced no plottable components.\n",
          "  n_pcs = %d | variance_explained = %s\n",
          "  Possible causes: the predictive component was not significant ",
          "(low R\u00b2X/Q\u00b2), insufficient between-group separation, or ",
          "cross-validation eliminated all components.\n",
          "  Returning NULL invisibly."
        ),
        method_label,
        res$n_pcs,
        paste(round(res$variance_explained, 3L), collapse = ", ")
      ),
      call. = FALSE
    )
    return(invisible(NULL))
  }

  # -- Auto-title ----------------------------------------------------------------
  if (is.null(title)) {
    if (res$method_used == "oplsda") {
      x_lbl <- if (pc_x == 1L) "t[1]" else sprintf("to[%d]", pc_x - 1L)
      y_lbl <- if (pc_y == 1L) "t[1]" else sprintf("to[%d]", pc_y - 1L)
      title <- sprintf("%s Scores Plot (%s vs %s)", method_label, x_lbl, y_lbl)
    } else {
      title <- sprintf("%s Scores Plot (Component %d vs Component %d)",
                       method_label, pc_x, pc_y)
    }
  }

  # -- Auto-subtitle: model fit stats (ropls only) -------------------------------
  if (is.null(subtitle) && res$method_used %in% c("oplsda", "plsda")) {
    subtitle <- tryCatch({
      mdl    <- res$pls_object
      df_mod <- mdl@modelDF

      .fmt_p <- function(p) {
        if (is.null(p) || length(p) == 0L || is.na(p)) return("NA")
        if (p < 0.001) return("< 0.001")
        sprintf("%.3f", round(p, 3L))
      }
      .fmt_r <- function(x) {
        if (is.null(x) || length(x) == 0L || is.na(x)) return("NA")
        sprintf("%.2f", round(x, 2L))
      }

      r2x <- if ("R2X(cum)" %in% colnames(df_mod)) {
        df_mod[nrow(df_mod), "R2X(cum)"]
      } else {
        sum(df_mod[, "R2X"], na.rm = TRUE)
      }

      r2y <- if ("R2Y(cum)" %in% colnames(df_mod)) {
        df_mod[nrow(df_mod), "R2Y(cum)"]
      } else if ("R2Y" %in% colnames(df_mod)) {
        sum(df_mod[, "R2Y"], na.rm = TRUE)
      } else {
        NA_real_
      }

      q2 <- if ("Q2(cum)" %in% colnames(df_mod)) {
        df_mod[nrow(df_mod), "Q2(cum)"]
      } else if ("Q2" %in% colnames(df_mod)) {
        df_mod[nrow(df_mod), "Q2"]
      } else {
        NA_real_
      }

      sum_df <- mdl@summaryDF
      p_r2y  <- if (!is.null(sum_df) && "pR2Y" %in% colnames(sum_df)) sum_df[1L, "pR2Y"] else NA_real_
      p_q2   <- if (!is.null(sum_df) && "pQ2"  %in% colnames(sum_df)) sum_df[1L, "pQ2"]  else NA_real_
      perm_n <- res$parameters$permI

      sprintf(
        "R\u00b2X = %s, R\u00b2Y = %s, Q\u00b2 = %s | perm = %d; p(R\u00b2Y) = %s, p(Q\u00b2) = %s",
        .fmt_r(r2x), .fmt_r(r2y), .fmt_r(q2),
        perm_n,
        .fmt_p(p_r2y), .fmt_p(p_q2)
      )
    }, error = function(e) NULL)
  }

  # -- Delegate to run_pca plotting engine ---------------------------------------
  p <- plot_score.run_pca(res = res, pc = pc, title = title, subtitle = subtitle, ...)

  ve <- res$variance_explained

  # -- Re-label axes by method ---------------------------------------------------
  if (res$method_used == "oplsda") {
    .fmt_axis <- function(idx, val) {
      if (idx == 1L) sprintf("Predictive Score t[1] (%.1f%%)", val)
      else           sprintf("Orthogonal Score to[%d] (%.1f%%)", idx - 1L, val)
    }
    p <- p + ggplot2::labs(
      x = .fmt_axis(pc_x, ve[pc_x]),
      y = .fmt_axis(pc_y, ve[pc_y])
    )
  } else {
    p <- p + ggplot2::labs(
      x = sprintf("Component %d (%.1f%%)", pc_x, ve[pc_x]),
      y = sprintf("Component %d (%.1f%%)", pc_y, ve[pc_y])
    )
  }

  return(p)
}

# ==============================================================================
# HELPER FUNCTIONS FOR OUTLIER DETECTION
# ==============================================================================

.compute_outliers_discrete <- function(df, pc1_col, pc2_col, ellipse_type, level) {
  mahal_d   <- rep(NA_real_, nrow(df))
  threshold <- rep(NA_real_, nrow(df))
  is_out    <- rep(NA,       nrow(df))

  for (grp in levels(df$ColorGroup)) {
    idx <- which(df$ColorGroup == grp)
    n   <- length(idx)
    if (n < 3L) next

    pts     <- as.matrix(df[idx, c(pc1_col, pc2_col), drop = FALSE])
    cov_mat <- tryCatch(stats::cov(pts), error = function(e) NULL)
    if (is.null(cov_mat) || det(cov_mat) <= .Machine$double.eps) next

    crit <- if (ellipse_type == "t") {
      stats::qt(level / 2 + 0.5, df = n - 1L)^2
    } else {
      stats::qnorm(level / 2 + 0.5)^2
    }

    d2 <- tryCatch(
      stats::mahalanobis(pts, colMeans(pts), cov_mat),
      error = function(e) rep(NA_real_, n)
    )

    mahal_d[idx]   <- d2
    threshold[idx] <- crit
    is_out[idx]    <- !is.na(d2) & (d2 > crit)
  }

  df$mahal_dist      <- mahal_d
  df$ellipse_threshold <- threshold
  df$is_outlier      <- is_out
  df
}

.resolve_outlier_groups_simple <- function(label_outliers, all_groups) {
  if ("all" %in% label_outliers) return(levels(all_groups))

  valid_groups <- label_outliers[label_outliers %in% levels(all_groups)]

  if (length(valid_groups) < length(label_outliers)) {
    invalid <- label_outliers[!label_outliers %in% levels(all_groups)]
    warning(
      sprintf("Groups in label_outliers not found: %s\n", paste(invalid, collapse = ", ")),
      sprintf("Available groups: %s", paste(levels(all_groups), collapse = ", "))
    )
  }

  unique(valid_groups)
}

#' @rdname plot_score
#' @param metadata A \code{data.frame} with one row per observation (sample).
#'   Required for \code{plot_score.run_pcoa} because \code{run_pcoa()} does not
#'   bundle metadata inside its result (unlike \code{run_pca()}). Row order must
#'   match the row order of the data matrix originally passed to
#'   \code{run_pcoa()}.
#' @export
plot_score.run_pcoa <- function(
    res,
    pc               = c(1, 2),
    color_by,
    points_from,
    metadata         = NULL, # Reordered to fix S3 consistency
    title            = NULL,
    subtitle         = NULL,
    caption          = NULL,
    position         = "center",
    arrange_levels   = NULL,
    ellipse          = TRUE,
    ellipse_type     = "t",
    ellipse_level    = 0.95,
    legend           = NULL,
    legend_position  = "bottom",
    show_outliers    = FALSE,
    label_outliers   = "all",
    colors           = NULL,
    point_size       = 3,
    point_alpha      = 0.8,
    theme            = "nature",
    base_size        = 11,
    font_family      = "sans",
    axis_title_size  = NULL,
    axis_text_size   = NULL,
    plot_title_size  = NULL,
    legend_title_size = NULL,
    legend_text_size  = NULL,
    zoom             = 1,
    verbose          = TRUE,
    ...
) {

  # -- 1. Validate metadata ----------------------------------------------------
  if (missing(metadata) || is.null(metadata)) {
    stop(
      "'metadata' is required for plot_score.run_pcoa().\n",
      "  Provide a data.frame with one row per observation, in the same row ",
      "order as the data originally passed to run_pcoa().",
      call. = FALSE
    )
  }

  if (!is.data.frame(metadata) && !inherits(metadata, "tbl_df")) {
    stop(
      "'metadata' must be a data.frame or tibble.\n",
      "  Received an object of class: ", paste(class(metadata), collapse = ", "),
      call. = FALSE
    )
  }

  n_obs <- res$n_obs
  if (nrow(metadata) != n_obs) {
    stop(
      "Row count mismatch: 'metadata' has ", nrow(metadata), " row(s) but the ",
      "PCoA result has ", n_obs, " observation(s).\n",
      "  Ensure 'metadata' has exactly one row per observation, in the same ",
      "order as the original data passed to run_pcoa().",
      call. = FALSE
    )
  }

  # -- 2. Build a pseudo run_pca-compatible object -----------------------------
  n_axes  <- res$n_axes
  pc_nms  <- colnames(res$scores)   # already "PC1", "PC2", ...

  proxy <- list(
    scores            = res$scores,
    metadata          = metadata,
    variance_explained = res$variance_explained,
    n_pcs             = n_axes
  )
  class(proxy) <- c("run_pca", "list")   # borrow run_pca dispatch for the body

  # -- 3. Patch axis title builder so labels read "PCo" not "PC" ---------------
  if (is.null(title)) {
    title <- sprintf("PCoA Scores Plot (PCoA%d vs PCoA%d)", pc[1], pc[2])
  }

  # -- 3.5 Auto-generate PERMANOVA & Post-Hoc subtitle if present --------------
  if (is.null(subtitle)) {
    sub_parts <- character(0)

    # Helper: format a p-value to 3 dp or "<0.001"
    .fmt_p_sub <- function(p) {
      if (is.na(p))      return("NA")
      if (p < 0.001)     return("< 0.001")
      sprintf("%.3f", round(p, 3L))
    }

    # 1. Main PERMANOVA line (always shown when present)
    if (!is.null(res$permanova)) {
      r2_val <- res$permanova$R2
      p_val  <- res$permanova$p_value
      sub_parts <- c(sub_parts, sprintf(
        "PERMANOVA: R\u00b2 = %.2f; p = %s",
        r2_val, .fmt_p_sub(p_val)
      ))
    }

    # 1.1 PERMDISP line
    if (!is.null(res$permdisp)) {
      p_val_disp <- res$permdisp$tab[["Pr(>F)"]][1]
      sub_parts <- c(sub_parts, sprintf(
        "Dispersion test (PERMDISP): p = %s",
        .fmt_p_sub(p_val_disp)
      ))
    }

    # 2. Determine number of levels in the rhs grouping variable
    n_levels <- NA_integer_
    if (!is.null(res$permanova) && !missing(metadata)) {
      rhs_col  <- res$permanova$rhs
      if (rhs_col %in% colnames(metadata)) {
        n_levels <- length(unique(metadata[[rhs_col]]))
      }
    }

    # 3. Post-hoc block - only for 3+ groups
    if (!is.null(res$pairwise_adonis) && !is.na(n_levels) && n_levels >= 3L) {

      pa       <- as.data.frame(res$pairwise_adonis)
      p_col    <- intersect(c("p.adj", "p.adjusted", "padj", "p.value", "pval"), colnames(pa))[1]
      pair_col <- intersect(c("pairs", "group", "contrast"),                      colnames(pa))[1]

      if (!is.na(p_col) && !is.na(pair_col)) {
        sig_rows <- pa[which(pa[[p_col]] < 0.05), , drop = FALSE]

        if (nrow(sig_rows) > 0L) {
          pair_lines <- vapply(seq_len(nrow(sig_rows)), function(i) {
            sprintf("%s; p = %s",
                    sig_rows[[pair_col]][i],
                    .fmt_p_sub(sig_rows[[p_col]][i]))
          }, character(1L))

          if (length(pair_lines) > 6L) {
            # Too many to list inline - redirect to summary()
            posthoc_str <- sprintf(
              "Post-hoc Sig. pairs: %d (p < 0.05) - see summary() for full table",
              length(pair_lines)
            )
          } else {
            posthoc_str <- paste0(
              "Post-hoc Sig. pairs:\n  ",
              paste(pair_lines, collapse = "\n  ")
            )
          }
          sub_parts <- c(sub_parts, posthoc_str)

        } else {
          sub_parts <- c(sub_parts, "Post-hoc: No sig. pairs (p < 0.05)")
        }
      }
    }
    # n_levels == 2: post-hoc block silently omitted; PERMANOVA p already shown above

    if (length(sub_parts) > 0L) {
      subtitle <- paste(sub_parts, collapse = "\n")
    }
  }

  # -- 4. Delegate to run_pca method -------------------------------------------
  p <- plot_score.run_pca(
    res               = proxy,
    pc                = pc,
    color_by          = color_by,
    points_from       = points_from,
    title             = title,
    subtitle          = subtitle,    # Automatically populated here
    caption           = caption,
    position          = position,
    arrange_levels    = arrange_levels,
    ellipse           = ellipse,
    ellipse_type      = ellipse_type,
    ellipse_level     = ellipse_level,
    legend            = legend,
    legend_position   = legend_position,
    show_outliers     = show_outliers,
    label_outliers    = label_outliers,
    colors            = colors,
    point_size        = point_size,
    point_alpha       = point_alpha,
    theme             = theme,
    base_size         = base_size,
    font_family       = font_family,
    axis_title_size   = axis_title_size,
    axis_text_size    = axis_text_size,
    plot_title_size   = plot_title_size,
    legend_title_size = legend_title_size,
    legend_text_size  = legend_text_size,
    zoom              = zoom,
    verbose           = verbose
  )

  # -- 5. Re-label axes to "PCo" convention ------------------------------------
  ve <- res$variance_explained
  p  <- p + ggplot2::labs(
    x = sprintf("PCoA%d (%.1f%%)", pc[1], ve[pc[1]]),
    y = sprintf("PCoA%d (%.1f%%)", pc[2], ve[pc[2]])
  )

  return(p)
}

.bin_continuous_safe <- function(x) {
  # Try quantile-based cut first; if breaks are non-unique (duplicate quantiles,
  # common with integer or low-cardinality continuous variables), fall back to
  # rank-based ntile-style binning, which always produces 4 non-empty bins.
  breaks <- unique(
    quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  )
  lbs <- c("Q1", "Q2", "Q3", "Q4")

  if (length(breaks) == 5L) {
    # Happy path: all four quantile boundaries are distinct
    result <- cut(x, breaks = breaks, labels = lbs, include.lowest = TRUE)
  } else {
    # Fallback: rank-based equal-frequency binning - always produces 4 groups
    # regardless of ties in the original variable.
    n      <- length(x)
    rnk    <- rank(x, ties.method = "first", na.last = "keep")
    bin    <- ceiling(rnk / n * 4L)
    bin    <- pmin(bin, 4L)           # guard against floating-point ceiling > 4
    result <- factor(lbs[bin], levels = lbs)
  }
  result
}

# Also update the .default method error message to include run_pcoa:
#' @rdname plot_score
#' @export
plot_score.default <- function(res, ...) {
  stop(
    "'res' must be output from run_pca(), run_pls(), or run_pcoa().\n",
    sprintf("Got object of class: %s\n", paste(class(res), collapse = ", ")),
    "Solution: Pass the result of run_pca(), run_pls(), or run_pcoa() directly.",
    call. = FALSE
  )
}