# Metabolomics workflow: `render_metabolomics()`

One function call per dataset. The Quarto workflow lives only in the package
(`inst/workflows/metabolomics_report.qmd`), so there is nothing to copy into projects.

```r
# separate files
pondeR::render_metabolomics("serum_features.csv", "serum_metadata.csv")

# one workbook, two sheets
pondeR::render_metabolomics("serum.xlsx", features_sheet = "features", metadata_sheet = "metadata",
                            qc_types = "EQC", injection_col = "Injection")

# urine: normalize by specific gravity (numeric metadata column; QC rows may be NA)
pondeR::render_metabolomics("urine.xlsx", features_sheet = "features", metadata_sheet = "metadata",
                            norm_factor_col = "SpecificGravity")

# one table holding metadata + features (extra non-feature columns via metadata_cols)
pondeR::render_metabolomics("neg.xlsx", features_sheet = "metadata_mother",
                            metadata_cols = c("SubjectID", "Replicate"))
```

- **Input:** csv, tsv/txt, xlsx. Samples in rows, features in columns (a features table with
  samples as columns is transposed). Features may carry a `Sample` column or follow the metadata
  row order. Metadata needs sample ID, group (including QC label), batch and injection order
  (names via `sample_id_col`, `group_col`, `batch_col`, `injection_col`).
- **Output:** `results/<name>_report.html` (standalone: tables and figures embedded) and
  `results/<name>_results.rds`, next to the features file (`output_dir` to change).
- **Overriding function arguments** goes through `...`, one named list per function, so the same
  argument name in different functions is never ambiguous:

  ```r
  render_metabolomics("serum.xlsx", features_sheet = "features", metadata_sheet = "metadata",
    run_foldchange   = list(arrange = c("Control", "Case"), up = 2),
    plot_volcano     = list(pval = 0.01),
    run_DIpreprocess = list(missing_threshold = 0.3))
  ```

  Valid names: `run_DIpreprocess`, `plot_beforeafter`, `plot_dist_beforeafter`, `run_pca`,
  `run_pls`, `plot_score`, `run_foldchange`, `run_diff`, `get_volcanodata`, `plot_volcano`,
  `plot_meanmap`, `run_auc`. Typos in function or argument names are rejected with the allowed
  list. Overrides used are listed in the report.
- **Update pondeR** (`devtools::install()`, then restart R) whenever the workflow changes;
  Quarto renders in a separate R process that uses the installed package.
- If rendering fails, rerun with `quiet = FALSE` for Quarto's messages.

## Behaviour to know

- Two groups: OPLS-DA (falls back to PLS-DA if no model), fold change, t-test, volcano, ROC.
- More than two groups: PLS-DA only; ANOVA; ROC for every pair of groups (each with its own
  pairwise PLS-DA, fold change and test). `control_group` applies to two-group designs only.
- Preprocessing runs before rendering, so data problems stop with an error immediately.
- No t-test results, or no features passing thresholds: mean abundance plot instead of volcano.
- ROC uses only features passing the VIP / fold change / p-value filters, so it can be empty.
