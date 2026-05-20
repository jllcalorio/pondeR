# Package index

## Direct-Injection Metabolomics Preprocessing Pipeline

Function for running a sequential preprocessing pipeline for
metabolomics data by orchestrating the individual run\_\* functions of
pondeR.

- [`run_DIpreprocess()`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)
  : Direct-Injection Metabolomics Preprocessing Pipeline

## Summarizing Data & Checking Distributions

Functions for generating publication-ready summary tables,
anthropometric indices according to WHO standards, and testing for
normality or skewness.

- [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  : Performs Descriptive Statistics on a Dataframe
- [`run_anthroindex()`](https://jllcalorio.github.io/pondeR/reference/run_anthroindex.md)
  : Compute WHO Anthropometric Z-Scores and Nutritional Status

## Comparing Groups & Identifying Associations

Methods for detecting differences between cohorts, calculating fold
changes, and testing statistical associations.

- [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  : Automatic Statistical Comparison with Comprehensive Analysis
- [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md)
  : Create Box or Violin Plots with Statistical Comparisons
- [`run_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  : Automatic Statistical Comparison for Categorical Variables
- [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  : Fold Change Analysis Across Groups
- [`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
  : Volcano Plot for Fold Change and p-value Data
- [`run_correl()`](https://jllcalorio.github.io/pondeR/reference/run_correl.md)
  : Perform Correlation Analysis
- [`plot_heatmap()`](https://jllcalorio.github.io/pondeR/reference/plot_heatmap.md)
  : Plot a Correlation Heatmap
- [`plot_dend()`](https://jllcalorio.github.io/pondeR/reference/plot_dend.md)
  : Plot a Hierarchical Clustering Dendrogram

## Predictive Modeling & Classification

Regularized regression (Ridge, LASSO, and Elastic Net), Linear
Mixed-Effect Model, Logistic regression analysis and performance
evaluation (AUC/AUROC) with automated bias correction options.

- [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  : Regularized Regression with Cross-Validation (Ridge, Elastic Net,
  LASSO)
- [`run_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md)
  : Perform Linear Mixed-Effects Model Analysis
- [`plot_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/plot_mixedmodel.md)
  : Plot Results from a Linear Mixed-Effects Model Analysis
- [`run_logreg()`](https://jllcalorio.github.io/pondeR/reference/run_logreg.md)
  : Perform Binary Logistic Regression (Standard or Firth-Corrected)
  with Iterative Feature Combinations
- [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md)
  : Comprehensive Area Under the ROC Curve (AUROC) Analysis

## Multivariate Exploration and Discrimination

Dimensionality reduction and visualization tools for simplified
exploration of complex datasets.

- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  : Perform Principal Component Analysis
- [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  : Perform Partial Least Squares (PLS) Analysis
- [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md)
  : Principal Coordinate Analysis (PCoA) with Optional PERMANOVA and
  Post-hoc
- [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md)
  : Plot PCA Scree Plot
- [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
  : Scores Plot for Multivariate Ordination Results
- [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md)
  : Reduce Multiple Data Frames by Common or Unique Column/Row Names

## Data Cleaning & Quality Control

Filtering features based on missingness, variance, or RSD thresholds.

- [`run_filtermissing()`](https://jllcalorio.github.io/pondeR/reference/run_filtermissing.md)
  : Filter Features by Missing Value Threshold
- [`run_filterRSD()`](https://jllcalorio.github.io/pondeR/reference/run_filterRSD.md)
  : Filter Features by Relative Standard Deviation in QC Samples
- [`run_filtervariance()`](https://jllcalorio.github.io/pondeR/reference/run_filtervariance.md)
  : Filter Features by Low Variance

## Missing Value Imputation

Addressing missing data using various imputation algorithms.

- [`run_mvimpute()`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md)
  : Impute Missing Values in Metabolomics Data

## Data Normalization, Transformation, and Scaling

Transforming and scaling data to ensure comparability across samples.

- [`run_normalize()`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md)
  : Normalize Metabolomics Data
- [`run_transform()`](https://jllcalorio.github.io/pondeR/reference/run_transform.md)
  : Transform Metabolomics Data
- [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md)
  : Scale Metabolomics Data

## Signal Drift and Batch Effects Correction

Correcting for technical variation and instrument signal drift.

- [`run_driftBatchCorrect()`](https://jllcalorio.github.io/pondeR/reference/run_driftBatchCorrect.md)
  : Correct Signal Drift and Batch Effects Using QC Samples

## Plotting Functions of Two Data Frames

Plots two data frames at a time, can be before and after implementing a
method.

- [`plot_beforeafter()`](https://jllcalorio.github.io/pondeR/reference/plot_beforeafter.md)
  : Plot Pairwise Before-and-After Comparison for Selected Features
- [`plot_dist_beforeafter()`](https://jllcalorio.github.io/pondeR/reference/plot_dist_beforeafter.md)
  : Plot Distribution Comparison Before and After Data Transformation

## Other Plotting Functions

Other functions to plot metabolomics and metagenomics data.

- [`plot_meanmap()`](https://jllcalorio.github.io/pondeR/reference/plot_meanmap.md)
  : Mean Intensity Heatmap with Hierarchical Clustering Dendrogram
- [`plot_multi_y()`](https://jllcalorio.github.io/pondeR/reference/plot_multi_y.md)
  : Plot Multiple Variables on Primary and Secondary Y Axes

## Export Functions

Functions export data frames to Excel/CSV files and figure to PNGs, etc.

- [`save_excel()`](https://jllcalorio.github.io/pondeR/reference/save_excel.md)
  : Save Data Frames, Tibbles, or Matrices to an Excel File
- [`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md)
  : Save Plots to Image or Vector Files

## Utilities & Helpers

Utility functions for data management, relabelling, and S3 methods.

- [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md)
  : Recode and Factorize a Column

- [`run_move_names()`](https://jllcalorio.github.io/pondeR/reference/run_move_names.md)
  : Move Row or Column Names into the Data Area

- [`run_randomize()`](https://jllcalorio.github.io/pondeR/reference/run_randomize.md)
  : Randomize a Data Frame, Matrix, or Vector

- [`get_groupsizes()`](https://jllcalorio.github.io/pondeR/reference/get_groupsizes.md)
  :

  Inspect Group Sample Sizes That May Trigger CV Fold Reduction in
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)

- [`get_removed_samples()`](https://jllcalorio.github.io/pondeR/reference/get_removed_samples.md)
  :

  Get Samples Removed by
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  Due to Missing Covariate Data

- [`get_sigfeatures()`](https://jllcalorio.github.io/pondeR/reference/get_sigfeatures.md)
  : Get Features with Non-Zero Coefficients from a Regularized
  Regression Result

- [`print(`*`<plot_volcano>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.plot_volcano.md)
  :

  Print Method for `plot_volcano` Objects

- [`print(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc.md)
  : Print method for run_assoc objects

- [`print(`*`<run_assoc_error>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc_error.md)
  : Print method for run_assoc_error objects

- [`print(`*`<run_assoc_multi>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc_multi.md)
  : Print method for run_assoc_multi objects

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

- [`summary(`*`<run_assoc_multi>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_assoc_multi.md)
  : Summary method for run_assoc_multi objects

- [`summary(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_diff.md)
  : Summary method for run_diff objects

- [`.assumption_row()`](https://jllcalorio.github.io/pondeR/reference/dot-assumption_row.md)
  : Build the standard assumptions data frame row from a normality
  result

- [`.build_posthoc_df()`](https://jllcalorio.github.io/pondeR/reference/dot-build_posthoc_df.md)
  : Build a standardised post-hoc result data frame

- [`.build_summary_table()`](https://jllcalorio.github.io/pondeR/reference/dot-build_summary_table.md)
  : Build a summary_table data frame from a list of run_diff objects

- [`.check_normality()`](https://jllcalorio.github.io/pondeR/reference/dot-check_normality.md)
  : Check normality of a numeric vector

- [`.check_variance()`](https://jllcalorio.github.io/pondeR/reference/dot-check_variance.md)
  : Run variance-homogeneity test (Levene or Bartlett fallback)

- [`.compute_effect_size()`](https://jllcalorio.github.io/pondeR/reference/dot-compute_effect_size.md)
  : Compute the appropriate effect size for a classical run_diff result

- [`.format_df()`](https://jllcalorio.github.io/pondeR/reference/dot-format_df.md)
  : Format a degrees-of-freedom string consistently across test types

- [`.parallel_normality()`](https://jllcalorio.github.io/pondeR/reference/dot-parallel_normality.md)
  : Run normality checks in parallel (or sequentially)

- [`.run_diff_single()`](https://jllcalorio.github.io/pondeR/reference/dot-run_diff_single.md)
  : Internal single-outcome run_diff for plot_diff

- [`.run_independent_anova()`](https://jllcalorio.github.io/pondeR/reference/dot-run_independent_anova.md)
  : Execute an independent ANOVA-family test

- [`.run_independent_two_sample()`](https://jllcalorio.github.io/pondeR/reference/dot-run_independent_two_sample.md)
  : Execute an independent two-sample test

- [`.run_paired_two_sample()`](https://jllcalorio.github.io/pondeR/reference/dot-run_paired_two_sample.md)
  : Execute a paired two-sample test

- [`.run_rm_anova()`](https://jllcalorio.github.io/pondeR/reference/dot-run_rm_anova.md)
  : Execute RM or mixed ANOVA and return a full run_diff-compatible
  result list

## Other

Additional exported functions. New functions land here automatically
until formally categorised.

- [`get_groupsizes()`](https://jllcalorio.github.io/pondeR/reference/get_groupsizes.md)
  :

  Inspect Group Sample Sizes That May Trigger CV Fold Reduction in
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)

- [`get_removed_samples()`](https://jllcalorio.github.io/pondeR/reference/get_removed_samples.md)
  :

  Get Samples Removed by
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  Due to Missing Covariate Data

- [`get_sigfeatures()`](https://jllcalorio.github.io/pondeR/reference/get_sigfeatures.md)
  : Get Features with Non-Zero Coefficients from a Regularized
  Regression Result

- [`has_inf()`](https://jllcalorio.github.io/pondeR/reference/has_inf.md)
  : Check if an Object Contains Infinite Values

- [`has_na()`](https://jllcalorio.github.io/pondeR/reference/has_na.md)
  : Check if an Object Contains Missing Values

- [`has_names()`](https://jllcalorio.github.io/pondeR/reference/has_names.md)
  : Check if an Object has Names

- [`has_negative()`](https://jllcalorio.github.io/pondeR/reference/has_negative.md)
  : Check if an Object Contains Negative Values

- [`has_zero()`](https://jllcalorio.github.io/pondeR/reference/has_zero.md)
  : Check if an Object Contains Zero Values

- [`is_list()`](https://jllcalorio.github.io/pondeR/reference/is_list.md)
  : Check if an Object is a Plain List

- [`is_tabular()`](https://jllcalorio.github.io/pondeR/reference/is_tabular.md)
  : Check if an Object is a Tabular Data Structure

- [`plot_beforeafter()`](https://jllcalorio.github.io/pondeR/reference/plot_beforeafter.md)
  : Plot Pairwise Before-and-After Comparison for Selected Features

- [`plot_dend()`](https://jllcalorio.github.io/pondeR/reference/plot_dend.md)
  : Plot a Hierarchical Clustering Dendrogram

- [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md)
  : Create Box or Violin Plots with Statistical Comparisons

- [`plot_dist_beforeafter()`](https://jllcalorio.github.io/pondeR/reference/plot_dist_beforeafter.md)
  : Plot Distribution Comparison Before and After Data Transformation

- [`plot_heatmap()`](https://jllcalorio.github.io/pondeR/reference/plot_heatmap.md)
  : Plot a Correlation Heatmap

- [`plot_meanmap()`](https://jllcalorio.github.io/pondeR/reference/plot_meanmap.md)
  : Mean Intensity Heatmap with Hierarchical Clustering Dendrogram

- [`plot_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/plot_mixedmodel.md)
  : Plot Results from a Linear Mixed-Effects Model Analysis

- [`plot_multi_y()`](https://jllcalorio.github.io/pondeR/reference/plot_multi_y.md)
  : Plot Multiple Variables on Primary and Secondary Y Axes

- [`plot_relabund()`](https://jllcalorio.github.io/pondeR/reference/plot_relabund.md)
  : Plot Relative Abundance as Stacked Bar Chart

- [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
  : Scores Plot for Multivariate Ordination Results

- [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md)
  : Plot PCA Scree Plot

- [`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
  : Volcano Plot for Fold Change and p-value Data

- [`print(`*`<plot_volcano>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.plot_volcano.md)
  :

  Print Method for `plot_volcano` Objects

- [`print(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc.md)
  : Print method for run_assoc objects

- [`print(`*`<run_assoc_error>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc_error.md)
  : Print method for run_assoc_error objects

- [`print(`*`<run_assoc_multi>`*`)`](https://jllcalorio.github.io/pondeR/reference/print.run_assoc_multi.md)
  : Print method for run_assoc_multi objects

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

- [`run_DIpreprocess()`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)
  : Direct-Injection Metabolomics Preprocessing Pipeline

- [`run_anthroindex()`](https://jllcalorio.github.io/pondeR/reference/run_anthroindex.md)
  : Compute WHO Anthropometric Z-Scores and Nutritional Status

- [`run_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  : Automatic Statistical Comparison for Categorical Variables

- [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md)
  : Comprehensive Area Under the ROC Curve (AUROC) Analysis

- [`run_correl()`](https://jllcalorio.github.io/pondeR/reference/run_correl.md)
  : Perform Correlation Analysis

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

- [`run_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md)
  : Perform Linear Mixed-Effects Model Analysis

- [`run_move_names()`](https://jllcalorio.github.io/pondeR/reference/run_move_names.md)
  : Move Row or Column Names into the Data Area

- [`run_mvimpute()`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md)
  : Impute Missing Values in Metabolomics Data

- [`run_normalize()`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md)
  : Normalize Metabolomics Data

- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  : Perform Principal Component Analysis

- [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md)
  : Principal Coordinate Analysis (PCoA) with Optional PERMANOVA and
  Post-hoc

- [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  : Perform Partial Least Squares (PLS) Analysis

- [`run_randomize()`](https://jllcalorio.github.io/pondeR/reference/run_randomize.md)
  : Randomize a Data Frame, Matrix, or Vector

- [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md)
  : Reduce Multiple Data Frames by Common or Unique Column/Row Names

- [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  : Regularized Regression with Cross-Validation (Ridge, Elastic Net,
  LASSO)

- [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md)
  : Recode and Factorize a Column

- [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md)
  : Scale Metabolomics Data

- [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  : Performs Descriptive Statistics on a Dataframe

- [`run_transform()`](https://jllcalorio.github.io/pondeR/reference/run_transform.md)
  : Transform Metabolomics Data

- [`save_excel()`](https://jllcalorio.github.io/pondeR/reference/save_excel.md)
  : Save Data Frames, Tibbles, or Matrices to an Excel File

- [`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md)
  : Save Plots to Image or Vector Files

- [`summary(`*`<run_assoc>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_assoc.md)
  : Summary method for run_assoc objects

- [`summary(`*`<run_assoc_multi>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_assoc_multi.md)
  : Summary method for run_assoc_multi objects

- [`summary(`*`<run_diff>`*`)`](https://jllcalorio.github.io/pondeR/reference/summary.run_diff.md)
  : Summary method for run_diff objects
