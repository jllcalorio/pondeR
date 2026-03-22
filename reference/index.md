# Package index

## Descriptive Statistics

Functions for generating publication-ready summary and descriptive
statistics tables.

- [`run_anthroindex()`](https://jllcalorio.github.io/pondeR/reference/run_anthroindex.md)
  : Compute WHO Anthropometric Z-Scores and Nutritional Status
- [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  : Performs Descriptive Statistics on a Dataframe

## Group Comparisons

Functions for comparing groups via parametric and non-parametric tests,
with automatic test selection and effect-size reporting.

- [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  : Automatic Statistical Comparison with Comprehensive Analysis
- [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md)
  : Create Box or Violin Plots with Statistical Comparisons
- [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  : Fold Change Analysis Across Groups
- [`run_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  : Automatic Statistical Comparison for Categorical Variables

## Regression & Classification

Logistic regression and AUC/AUROC analysis functions with
publication-ready outputs and optional Firth correction.

- [`run_logreg()`](https://jllcalorio.github.io/pondeR/reference/run_logreg.md)
  : Perform Binary Logistic Regression (Standard or Firth-Corrected)
  with Iterative Feature Combinations
- [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md)
  : Comprehensive Area Under the ROC Curve (AUROC) Analysis

## Dimensionality Reduction

PCA-based dimensionality reduction, scree plots, and score plots for
multivariate data exploration.

- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  : Perform Principal Component Analysis
- [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md)
  : Reduce Multiple Data Frames by Common or Unique Column/Row Names
- [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md)
  : Plot PCA Scree Plot
- [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
  : Plot PCA Scores Plot

## Visualisation

Standalone plotting functions for volcano plots and other figures.

- [`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
  : Volcano Plot for Fold Change and P-value Data

## Data Pre-processing

Metabolomics-focused data pre-processing utilities including imputation,
normalisation, filtering, scaling, and batch correction.

- [`run_mvimpute()`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md)
  : Impute Missing Values in Metabolomics Data
- [`run_normalize()`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md)
  : Normalize Metabolomics Data
- [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md)
  : Scale Metabolomics Data
- [`run_transform()`](https://jllcalorio.github.io/pondeR/reference/run_transform.md)
  : Transform Metabolomics Data
- [`run_filterRSD()`](https://jllcalorio.github.io/pondeR/reference/run_filterRSD.md)
  : Filter Features by Relative Standard Deviation in QC Samples
- [`run_filtermissing()`](https://jllcalorio.github.io/pondeR/reference/run_filtermissing.md)
  : Filter Features by Missing Value Threshold
- [`run_filtervariance()`](https://jllcalorio.github.io/pondeR/reference/run_filtervariance.md)
  : Filter Features by Low Variance
- [`run_driftBatchCorrect()`](https://jllcalorio.github.io/pondeR/reference/run_driftBatchCorrect.md)
  : Correct Signal Drift and Batch Effects Using QC Samples

## Utilities & Helpers

Utility functions for data management, relabelling, and S3 methods.

- [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md)
  : Recode and Factorize a Column

- [`print(`*`<plot_volcano>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.plot_volcano.md)
  :

  Print Method for `plot_volcano` Objects

- [`print(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc.md)
  : Print method for run_assoc objects

- [`print(`*`<run_auc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_auc.md)
  :

  Print Method for `run_auc` Objects

- [`print(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_diff.md)
  : Print method for run_diff objects

- [`print(`*`<run_foldchange>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_foldchange.md)
  :

  Print Method for `run_foldchange` Objects

- [`print(`*`<run_reduce>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_reduce.md)
  :

  Print Method for `run_reduce` Objects

- [`summary(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_assoc.md)
  : Summary method for run_assoc objects

- [`summary(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_diff.md)
  : Summary method for run_diff objects

## Other

Additional exported functions. New functions land here automatically
until formally categorised.

- [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md)
  : Create Box or Violin Plots with Statistical Comparisons

- [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
  : Plot PCA Scores Plot

- [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md)
  : Plot PCA Scree Plot

- [`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
  : Volcano Plot for Fold Change and P-value Data

- [`print(`*`<plot_volcano>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.plot_volcano.md)
  :

  Print Method for `plot_volcano` Objects

- [`print(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc.md)
  : Print method for run_assoc objects

- [`print(`*`<run_auc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_auc.md)
  :

  Print Method for `run_auc` Objects

- [`print(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_diff.md)
  : Print method for run_diff objects

- [`print(`*`<run_foldchange>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_foldchange.md)
  :

  Print Method for `run_foldchange` Objects

- [`print(`*`<run_reduce>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_reduce.md)
  :

  Print Method for `run_reduce` Objects

- [`run_anthroindex()`](https://jllcalorio.github.io/pondeR/reference/run_anthroindex.md)
  : Compute WHO Anthropometric Z-Scores and Nutritional Status

- [`run_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  : Automatic Statistical Comparison for Categorical Variables

- [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md)
  : Comprehensive Area Under the ROC Curve (AUROC) Analysis

- [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  : Automatic Statistical Comparison with Comprehensive Analysis

- [`run_driftBatchCorrect()`](https://jllcalorio.github.io/pondeR/reference/run_driftBatchCorrect.md)
  : Correct Signal Drift and Batch Effects Using QC Samples

- [`run_filterRSD()`](https://jllcalorio.github.io/pondeR/reference/run_filterRSD.md)
  : Filter Features by Relative Standard Deviation in QC Samples

- [`run_filtermissing()`](https://jllcalorio.github.io/pondeR/reference/run_filtermissing.md)
  : Filter Features by Missing Value Threshold

- [`run_filtervariance()`](https://jllcalorio.github.io/pondeR/reference/run_filtervariance.md)
  : Filter Features by Low Variance

- [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  : Fold Change Analysis Across Groups

- [`run_logreg()`](https://jllcalorio.github.io/pondeR/reference/run_logreg.md)
  : Perform Binary Logistic Regression (Standard or Firth-Corrected)
  with Iterative Feature Combinations

- [`run_mvimpute()`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md)
  : Impute Missing Values in Metabolomics Data

- [`run_normalize()`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md)
  : Normalize Metabolomics Data

- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  : Perform Principal Component Analysis

- [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md)
  : Reduce Multiple Data Frames by Common or Unique Column/Row Names

- [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md)
  : Recode and Factorize a Column

- [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md)
  : Scale Metabolomics Data

- [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  : Performs Descriptive Statistics on a Dataframe

- [`run_transform()`](https://jllcalorio.github.io/pondeR/reference/run_transform.md)
  : Transform Metabolomics Data

- [`summary(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_assoc.md)
  : Summary method for run_assoc objects

- [`summary(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_diff.md)
  : Summary method for run_diff objects
