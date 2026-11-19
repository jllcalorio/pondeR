#' @importFrom magrittr %>%
#' @importFrom rlang %||%
#' @importFrom utils globalVariables tail
#' @importFrom stats binom.test
NULL

# Global variables to avoid R CMD check notes (ggplot2/dplyr NSE column names)
utils::globalVariables(c(
  ".", "group1", "group2", "p.adj", "p.adj.signif", "y",
  "ColorGroup", "ShapeGroup", "Label", "PC_num",
  "Value", "CumValue", "var",
  "xend", "yend", "label", # plot_dend, plot_meanmap
  "agg_val", "agg_lbl", # plot_diff
  "Var1", "Var2", "estimate", # plot_heatmap
  "seg_col", "ymin", "ymax", "color", "variable", "x_pos", # plot_meanmap
  ".neg_log10_p", ".sig", ".label", # plot_mixedmodel
  "fpr", "tpr", "predictor", # run_auc
  "prop", "group", # plot_bar
  ".data", # standard for ggplot2/dplyr
  "value" # plot_meanmap
))
