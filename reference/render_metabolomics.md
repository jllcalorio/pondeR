# Render a standalone metabolomics analysis report

Runs the full direct-injection metabolomics workflow (preprocessing,
batch correction views, PCA, OPLS-DA / PLS-DA, fold change, t-test /
ANOVA, volcano plot, ROC) on one dataset and writes a single
self-contained HTML report (all tables and figures embedded) plus an
`.rds` of the results to a `results` folder next to the input file. The
workflow lives in the package, so projects never hold their own copy.

## Usage

``` r
render_metabolomics(
  features,
  metadata = NULL,
  features_sheet = NULL,
  metadata_sheet = NULL,
  metadata_cols = NULL,
  name = NULL,
  output_dir = NULL,
  sample_id_col = "Sample",
  group_col = "Group",
  batch_col = "Batch",
  injection_col = "InjectionSequence",
  qc_types = c("QC", "SQC", "EQC"),
  normalize_method = "median",
  norm_factor_col = NULL,
  normalize_qc_method = "none",
  control_group = NULL,
  ...,
  overwrite = FALSE,
  timestamp = FALSE,
  quiet = TRUE
)
```

## Arguments

- features:

  Path to the features file (or the combined file), or a data frame.

- metadata:

  Path to the metadata file, or a data frame. Omit when `metadata_sheet`
  is in the `features` workbook, or for a combined table.

- features_sheet, metadata_sheet:

  Sheet names (workbooks only).

- metadata_cols:

  Extra metadata column names for a combined table.

- name:

  Report name; default the features file name.

- output_dir:

  Folder that will contain `results/`. Default: the folder of the
  `features` file.

- sample_id_col, group_col, batch_col, injection_col:

  Metadata column names.

- qc_types:

  Group labels identifying QC samples.

- normalize_method:

  Normalization method (e.g. `"median"`, `"sum"`, `"pqn_global"`); see
  [`dimsprepr::run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html).

- norm_factor_col:

  Numeric metadata column with per-sample normalization values (e.g.
  specific gravity for urine). When given, it overrides
  `normalize_method`.

- normalize_qc_method:

  How QC samples are normalized when `norm_factor_col` is used (QC rows
  usually have no value): `"median"`, `"mean"` or `"none"` (default).
  Ignored without `norm_factor_col`.

- control_group:

  Reference group for ROC (two-group designs).

- ...:

  Per-function argument lists, see above.

- overwrite:

  If `FALSE` (default) an existing report is never replaced: a timestamp
  (`_YYYYmmdd-HHMMSS`) is appended to the new file names. `TRUE`
  replaces the existing report and results.

- timestamp:

  If `TRUE`, always append a timestamp to the file names.

- quiet:

  Set `FALSE` to see Quarto's progress output.

## Value

Invisibly, the path of the HTML report.

## Data layout

Samples are rows and features are columns in `features`; `metadata` has
one row per sample. Supply either

- two files (`features`, `metadata`), or one workbook with two sheets
  (`features_sheet`, `metadata_sheet`), or

- one table holding metadata and features together (give only
  `features`); all non-numeric columns plus `sample_id_col`,
  `group_col`, `batch_col`, `injection_col` and `metadata_cols` are
  treated as metadata, the rest as features.

Supported formats: `.csv`, `.tsv`/`.txt`, `.xlsx`/`.xls`. A data frame
can also be given instead of a path (then set `output_dir`). Features
may carry a sample ID column, or be in the same row order as `metadata`;
a features table with samples as columns is transposed automatically.

## Overriding function arguments

Settings of the underlying functions are passed through `...` as one
named list per function, so the same argument name in different
functions is never ambiguous:

    render_metabolomics("serum.xlsx", features_sheet = "features",
      metadata_sheet = "metadata",
      run_foldchange = list(arrange = c("Control", "Case"), up = 2),
      plot_volcano   = list(pval = 0.01),
      run_DIpreprocess = list(missing_threshold = 0.3))

Valid names: `run_DIpreprocess`, `plot_beforeafter`,
`plot_dist_beforeafter`,
[`run_pca`](https://jllcalorio.github.io/pondeR/reference/run_pca.md),
[`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md),
[`plot_score`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
(PCA and PLS plots), `plot_score.pca` and `plot_score.pls` (only that
plot; override
[`plot_score`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)),
[`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
[`get_volcanodata`](https://jllcalorio.github.io/pondeR/reference/get_volcanodata.md),
[`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md),
[`plot_meanmap`](https://jllcalorio.github.io/pondeR/reference/plot_meanmap.md),
[`run_auc`](https://jllcalorio.github.io/pondeR/reference/run_auc.md).
Argument names are checked against the function, and a supplied value
replaces the workflow's default (including `NULL`). Data inputs (`x`)
cannot be overridden.
