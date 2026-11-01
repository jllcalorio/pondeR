# Config for pondeR::render_metabolomics(). Copy into a project folder, edit, then run:
#   pondeR::render_metabolomics("config.R")
# Paths are relative to this file's folder. This is plain R: load and reshape data however you need.

# Optional: override defaults for ALL datasets (any entry can also be set per dataset)
common <- list(
  sample_id_col = "Sample",
  group_col     = "Group",
  batch_col     = "Batch",
  injection_col = "InjectionSequence",
  qc_types      = "QC",             # all QC labels used in the Group column
  control_group = NULL,             # ROC reference level (NULL = first alphabetically)
  fc_up         = 1.5,
  fc_down       = 1 / 1.5,
  extra         = list()            # other run_DIpreprocess() args, e.g. list(min_qc_per_batch = 3)
)

# Required: one named entry per dataset. `features`: raw intensities, samples in rows.
# `metadata`: same rows, in the same order. Each may be a data frame or a .csv/.xlsx path.
dat <- openxlsx::read.xlsx("data-cleaned-v2/my_data.xlsx", "my_sheet")
meta_cols <- c("Sample", "Group", "Batch", "InjectionSequence")  # everything else = features

datasets <- list(
  Serum = list(
    metadata = dat[meta_cols],
    features = dat[setdiff(names(dat), meta_cols)],
    normalize_method = "median"
  )
  # Urine = list(
  #   metadata = ..., features = ...,
  #   norm_factor_col = "SpecificGravity",   # numeric metadata column; overrides normalize_method
  #   normalize_qc_method = "median"         # how to treat QC rows (no SG value)
  # )
)

out_dir <- "results"
