#' Plot PCA Scree Plot
#'
#' @description
#' Creates a scree plot showing variance explained or eigenvalues per principal component.
#' Helps determine the number of meaningful PCs to retain.
#'
#' @param pca_result List. Output from `run_pca()`.
#' @param title Character. Main plot title. Default: "Scree Plot".
#' @param subtitle Character. Plot subtitle. Default: NULL (no subtitle).
#' @param position Character. Position of title and subtitle. Options: "left", "center", "right".
#'   Default: "center".
#' @param what Character. What to plot: "variance" (variance explained %) or "eigen" (eigenvalues).
#'   Default: "variance".
#' @param type Character. Plot type: "bar", "line", or "both". Default: "line".
#' @param max_pc Integer. Maximum number of PCs to display. Default: 15.
#' @param show_cumulative Logical. Show cumulative variance explained as a line on secondary axis.
#'   Only works when `what = "variance"`. Default: TRUE.
#' @param show_values Logical. Display values on plot using ggrepel to avoid overlap. Default: TRUE.
#' @param theme_base_size Numeric. Base font size for theme. Default: 11.
#' @param verbose Logical. Print messages. Default: TRUE.
#'
#' @return A ggplot2 object.
#'
#' @details
#' **Scree Plot Interpretation:**
#' 
#' The scree plot helps determine how many principal components to retain:
#' 
#' - **Elbow Method**: Look for an "elbow" where the curve flattens. PCs before
#'   the elbow capture meaningful variance; those after capture mostly noise.
#'   
#' - **Kaiser Criterion**: Retain PCs with eigenvalues > 1 (when `what = "eigen"`).
#'   These PCs explain more variance than a single original variable.
#'   
#' - **Cumulative Variance**: Retain enough PCs to explain 70-90% of total variance
#'   (shown when `show_cumulative = TRUE`).
#' 
#' **Plot Types:**
#' 
#' - **bar**: Traditional bar chart, good for discrete comparison
#' - **line**: Line plot, better for seeing trends
#' - **both**: Combines bars and line for comprehensive view
#' 
#' **Cumulative Variance Line:**
#' 
#' When `show_cumulative = TRUE`, a line showing cumulative variance is added
#' to the plot on a secondary y-axis at the top. This helps visualize how
#' many PCs are needed to explain a target percentage of variance (e.g., 80%).
#' 
#' **When to Use Each Metric:**
#' 
#' - **Variance explained**: More intuitive (percentages). Recommended for most users.
#' - **Eigenvalues**: Useful for Kaiser criterion and advanced analysis.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Cattell, R.B. (1966). The scree test for the number of factors.
#' Multivariate Behavioral Research, 1(2), 245-276. \doi{10.1207/s15327906mbr0102_10}
#' 
#' Kaiser, H.F. (1960). The application of electronic computers to factor analysis.
#' Educational and Psychological Measurement, 20(1), 141-151.
#' \doi{10.1177/001316446002000116}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_point scale_x_continuous scale_y_continuous sec_axis labs theme_minimal theme element_text element_blank element_line margin
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run PCA first
#' pca_result <- run_pca(scaled_data, metadata, scale_method = "auto")
#' 
#' # Basic scree plot
#' plot_scree(pca_result)
#' 
#' # Eigenvalue plot with custom styling
#' plot_scree(pca_result, 
#'           what = "eigen",
#'           type = "both",
#'           title = "PCA Eigenvalues",
#'           position = "left")
#' 
#' # Variance plot without cumulative line
#' plot_scree(pca_result,
#'           show_cumulative = FALSE,
#'           max_pc = 10)
#' }
plot_scree <- function(
    pca_result,
    title = "Scree Plot",
    subtitle = NULL,
    position = "center",
    what = "variance",
    type = "line",
    max_pc = 15,
    show_cumulative = TRUE,
    show_values = TRUE,
    theme_base_size = 11,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if (!inherits(pca_result, "run_pca")) {
    stop("'pca_result' must be output from run_pca() function.\n",
         "Solution: Run run_pca() first on your data before plotting.")
  }
  
  what <- tolower(what)
  type <- tolower(type)
  position <- tolower(position)
  
  if (!what %in% c("variance", "eigen")) {
    stop("'what' must be either 'variance' or 'eigen'.\n",
         "Solution: Set what='variance' (default) or what='eigen'.")
  }
  if (!type %in% c("bar", "line", "both")) {
    stop("'type' must be 'bar', 'line', or 'both'.\n",
         "Solution: Choose one of: type='line' (default), type='bar', or type='both'.")
  }
  if (!position %in% c("left", "center", "right")) {
    stop("'position' must be 'left', 'center', or 'right'.\n",
         "Solution: Use position='center' (default), position='left', or position='right'.")
  }
  if (!is.numeric(max_pc) || length(max_pc) != 1 || max_pc < 1) {
    stop("'max_pc' must be a positive integer.\n",
         "Solution: Set max_pc to a number between 1 and ", pca_result$n_pcs, 
         " (total PCs available).")
  }
  if (!is.character(title) || length(title) != 1) {
    stop("'title' must be a single character string.\n",
         "Solution: Provide a title like: title='My Scree Plot'")
  }
  if (!is.null(subtitle) && (!is.character(subtitle) || length(subtitle) != 1)) {
    stop("'subtitle' must be NULL or a single character string.\n",
         "Solution: Either omit subtitle or provide: subtitle='My subtitle'")
  }
  if (!is.numeric(theme_base_size) || length(theme_base_size) != 1 || theme_base_size <= 0) {
    stop("'theme_base_size' must be a positive numeric value.\n",
         "Solution: Use theme_base_size=11 (default) or another positive number.")
  }
  
  if (show_cumulative && what != "variance") {
    warning("'show_cumulative' only works with what='variance'. Setting show_cumulative to FALSE.\n",
            "Solution: Use what='variance' if you want to show cumulative variance.")
    show_cumulative <- FALSE
  }
  
  # Check for ggrepel
  if (show_values && !requireNamespace("ggrepel", quietly = TRUE)) {
    warning("Package 'ggrepel' not available. Using regular geom_text instead (labels may overlap).\n",
            "Solution: Install ggrepel with: install.packages('ggrepel')")
    use_ggrepel <- FALSE
  } else {
    use_ggrepel <- show_values
  }
  
  # Prepare data
  n_pcs <- pca_result$n_pcs
  max_show <- min(max_pc, n_pcs)
  
  if (max_pc > n_pcs) {
    warning(sprintf("Requested max_pc=%d but only %d PCs available. Showing %d PCs.\n",
                    max_pc, n_pcs, n_pcs),
            "Solution: Reduce max_pc to a value <= ", n_pcs)
  }
  
  if (what == "variance") {
    y_values <- pca_result$variance_explained[1:max_show]
    y_label <- "Variance Explained (%)"
  } else {
    y_values <- pca_result$eigenvalues[1:max_show]
    y_label <- "Eigenvalue"
  }
  
  scree_data <- data.frame(
    PC_num = 1:max_show,
    Value = y_values,
    CumValue = if (what == "variance") pca_result$cumulative_variance[1:max_show] else NA,
    stringsAsFactors = FALSE
  )
  
  msg(sprintf("Creating scree plot (%s, type: %s, max PC: %d)", what, type, max_show))
  
  # Base plot with numeric x-axis
  p <- ggplot2::ggplot(scree_data, ggplot2::aes(x = PC_num, y = Value))
  
  # Add plot elements based on type
  if (type == "bar") {
    p <- p + ggplot2::geom_col(fill = "steelblue", alpha = 0.7, width = 0.6)
  } else if (type == "line") {
    p <- p + 
      ggplot2::geom_line(color = "steelblue", linewidth = 1.2, group = 1) +
      ggplot2::geom_point(color = "steelblue", size = 3)
  } else {  # both
    p <- p +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.5, width = 0.6) +
      ggplot2::geom_line(color = "darkblue", linewidth = 1.2, group = 1) +
      ggplot2::geom_point(color = "darkblue", size = 3)
  }
  
  # Add value labels with ggrepel
  if (show_values) {
    if (what == "variance") {
      val_labels <- sprintf("%.1f", scree_data$Value)
    } else {
      val_labels <- sprintf("%.2f", scree_data$Value)
    }
    
    label_data <- data.frame(
      PC_num = scree_data$PC_num,
      Value = scree_data$Value,
      Label = val_labels
    )
    
    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = Label),
        size = theme_base_size / 3.5,  # Scale with base size
        color = "black",
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.size = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        direction = "y",
        nudge_y = max(scree_data$Value) * 0.05
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = Label),
        vjust = -0.5,
        size = theme_base_size / 3.5,  # Scale with base size
        color = "black"
      )
    }
  }
  
  # Determine title justification
  hjust_val <- switch(position,
                      "left" = 0,
                      "center" = 0.5,
                      "right" = 1,
                      0.5)
  
  # Add cumulative variance line on secondary axis (at top)
  if (show_cumulative && what == "variance") {
    cum_data <- data.frame(
      PC_num = scree_data$PC_num,
      CumValue = scree_data$CumValue
    )
    
    cum_labels <- sprintf("%.1f", cum_data$CumValue)
    
    # Add cumulative line
    p <- p + 
      ggplot2::geom_line(
        data = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue),
        color = "darkred",
        linewidth = 1,
        linetype = "dashed",
        group = 1
      ) +
      ggplot2::geom_point(
        data = cum_data,
        ggplot2::aes(x = PC_num, y = CumValue),
        color = "darkred",
        size = 2.5
      )
    
    # Add cumulative labels with ggrepel
    cum_label_data <- data.frame(
      PC_num = cum_data$PC_num,
      CumValue = cum_data$CumValue,
      Label = cum_labels
    )
    
    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = cum_label_data,
        ggplot2::aes(x = PC_num, y = CumValue, label = Label),
        size = theme_base_size / 4,  # Scale with base size, slightly smaller
        color = "darkred",
        fontface = "italic",
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.size = 0.2,
        min.segment.length = 0,
        max.overlaps = Inf,
        direction = "y",
        nudge_y = max(cum_data$CumValue) * 0.02
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = cum_label_data,
        ggplot2::aes(x = PC_num, y = CumValue, label = Label),
        vjust = -0.5,
        size = theme_base_size / 4,  # Scale with base size
        color = "darkred",
        fontface = "italic"
      )
    }
    
    # Add secondary axis label for cumulative variance
    p <- p + 
      ggplot2::scale_y_continuous(
        name = y_label,
        sec.axis = ggplot2::sec_axis(
          ~ .,
          name = "Cumulative Variance (%)"
        )
      )
  }
  
  # Add labels and theme
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = 1:max_show,
      labels = as.character(1:max_show)
    ) +
    ggplot2::labs(
      x = "Principal Component",
      y = if (!show_cumulative || what != "variance") y_label else NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = theme_base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = hjust_val, size = theme_base_size + 3, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = hjust_val, size = theme_base_size + 1),
      axis.title.x = ggplot2::element_text(size = theme_base_size + 1),
      axis.title.y = ggplot2::element_text(size = theme_base_size + 1),
      axis.title.y.right = if (show_cumulative && what == "variance") {
        ggplot2::element_text(size = theme_base_size + 1, color = "darkred")
      } else {
        ggplot2::element_blank()
      },
      axis.text.x = ggplot2::element_text(size = theme_base_size - 1),
      axis.text.y = if (show_values && !show_cumulative) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_text(size = theme_base_size - 1)
      },
      axis.text.y.right = if (show_cumulative && what == "variance") {
        ggplot2::element_text(size = theme_base_size - 1, color = "darkred")
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.y = if (show_values && !show_cumulative) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_line()
      },
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  msg("Scree plot created successfully")
  return(p)
}