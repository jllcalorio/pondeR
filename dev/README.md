# Metabolomics workflow

Files (keep together):

- `metabolomics_workflow.qmd` - main document; edit only its **Configuration** section.
- `_metabolomics_dataset.qmd` - per-dataset analysis, included once per dataset.

## Steps

1. **Install** (once): `devtools::install("<path to pondeR repo>")` so the render uses the
   current pondeR; also needs `dimsprepr`, `gt`, `knitr`, `openxlsx`.
2. **New analysis folder:** copy both `.qmd` files into it, add your data, and open the folder
   in Positron (outputs go to `results/` here).
3. **Prepare each dataset:**
   - `features`: raw intensities, samples in rows, features in columns.
   - `metadata`: same rows/order, with Sample (unique, matches feature rows), Group (incl. QC
     label), Batch, InjectionSequence. Urine also needs a numeric `SpecificGravity` column
     (QC rows may be `NA`).
4. **Configure** in `metabolomics_workflow.qmd`:
   - `use_demo <- FALSE`.
   - Set each dataset's `features` / `metadata` (data frame in session, or `.csv`/`.xlsx` path).
   - Check `common`: column names, `qc_types`, `control_group` (ROC reference level; for
     more than two groups it is ignored and the first group of each pair is the control),
     `fc_up` / `fc_down`.
   - Serum: `normalize_method = "median"`. Urine: `norm_factor_col = "SpecificGravity"`
     (overrides `normalize_method`) and `normalize_qc_method`.
   - Extra `run_DIpreprocess()` arguments go in `extra = list(...)`, e.g.
     `list(min_qc_per_batch = 3)` or `list(rsd_threshold = 0.4)`.
   - More datasets: copy an entry in `datasets`.
5. **Render:** `quarto::quarto_render("metabolomics_workflow.qmd")`. Output: one HTML report
   (section per dataset) and `results/<dataset>/results.rds`.

## Behaviour to know

- Two groups: OPLS-DA (falls back to PLS-DA if no model), fold change, t-test, volcano, ROC.
- More than two groups: PLS-DA only; ANOVA; ROC for every pair of groups (each with its own
  pairwise PLS-DA, fold change and test). OPLS-DA is not run.
- A dataset that fails preprocessing is reported and skipped; others still run. Common causes:
  quality filter removes all features, or too few QCs per batch for drift correction.
- No t-test results, or no features passing thresholds: mean abundance plot instead of volcano.
- ROC uses only features passing the VIP / fold change / p-value filters, so it can be empty.
- The demo serum set has 3 groups and the urine set has 2, to exercise both paths.
