#' @importFrom magrittr %>%
#' @importFrom utils globalVariables tail
#' @importFrom stats binom.test
NULL

#' Global variables to avoid R CMD check notes
#'
#' @noRd
#'
globalVariables(c(
  ".", "group1", "group2", "p.adj", "p.adj.signif", "y",
  "ColorGroup", "ShapeGroup", "Label", "PC_num", 
  "Value", "CumValue", "daxta", "var", "::", ":::",

  "injection", "value", "group", "batch", # plot_beforeafter
  "xend", "yend", "label", # plot_dend, plot_meanmap
  "agg_val", "agg_lbl", # plot_diff
  "Var1", "Var2", "estimate", # plot_heatmap
  "seg_col", "ymin", "ymax", "color", "variable", "x_pos", # plot_meanmap
  ".neg_log10_p", ".sig", ".label", # plot_mixedmodel
  "fpr", "tpr", "predictor", # run_auc
  "export_xl_workbook", # save_excel
  ".data" # standard for ggplot2/dplyr
))
