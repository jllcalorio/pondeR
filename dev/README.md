# Metabolomics workflow

The workflow lives **only in the pondeR package** (`inst/workflows/`). A project folder holds
just its data and a small `config.R`; nothing is copied, so updating pondeR updates every project.

## Using it for a project

1. **Update pondeR** whenever the workflow changes: `devtools::install("<path to pondeR repo>")`
   (also needs `dimsprepr`, `gt`, `knitr`, `quarto`, `openxlsx`). Restart R after installing.
2. **Create `config.R`** in the project folder (e.g. `P1_Serum/config.R`), starting from the template:
   `file.copy(system.file("workflows/config_template.R", package = "pondeR"), "config.R")`
3. **Edit `config.R`** (plain R; paths relative to its folder):
   - Load the data and define `datasets`: one named entry per dataset, each with `features`
     (raw intensities, samples in rows) and `metadata` (same rows, same order), plus
     `normalize_method = "median"` (serum) or `norm_factor_col = "SpecificGravity"` +
     `normalize_qc_method` (urine).
   - `common` sets defaults for all datasets: `group_col`, `batch_col`, `injection_col`,
     `qc_types`, `control_group` (ROC reference level), `fc_up`/`fc_down`, `extra`
     (other `run_DIpreprocess()` arguments, e.g. `list(min_qc_per_batch = 3)`).
   - Metadata needs a unique sample ID, the group (including the QC label), batch and injection order.
4. **Render:** `pondeR::render_metabolomics("P1_Serum/config.R")`
   Output goes next to the config: `metabolomics_report.html` and `results/<dataset>/results.rds`.
   Use `quiet = FALSE` to see Quarto's messages if it fails. Several datasets in one config
   give one section each; for another project, make another `config.R`.
   `render_metabolomics()` with no config renders the demo data.

## Behaviour to know

- Two groups: OPLS-DA (falls back to PLS-DA if no model), fold change, t-test, volcano, ROC.
- More than two groups: PLS-DA only; ANOVA; ROC for every pair of groups (each with its own
  pairwise PLS-DA, fold change and test). `control_group` is ignored there (first group of each
  pair is the control).
- A dataset that fails preprocessing is reported and skipped; others still run. Common causes:
  quality filter removes all features, or too few QCs per batch for drift correction.
- No t-test results, or no features passing thresholds: mean abundance plot instead of volcano.
- ROC uses only features passing the VIP / fold change / p-value filters, so it can be empty.
- `run_DIpreprocess` standardizes metadata names (`Group`, `Sample`, `Batch`, `InjectionSequence`);
  your own column names go in the config and are handled for you.
- The demo serum set has 3 groups and the urine set has 2, to exercise both paths.
