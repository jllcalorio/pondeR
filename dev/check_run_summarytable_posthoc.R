# Self-check for run_summarytable(mark_posthoc_significance = TRUE).
# Run with: Rscript dev/check_run_summarytable_posthoc.R  (from package root)
devtools::load_all(".", quiet = TRUE)

set.seed(4821)
n <- 300
df <- data.frame(Group = sample(c("A", "B", "C"), n, replace = TRUE))
df$Age <- NA_real_
# Large, RNG-robust mean separations: A and C similar, B very different.
df$Age[df$Group == "A"] <- rnorm(sum(df$Group == "A"), 30, 5)
df$Age[df$Group == "B"] <- rnorm(sum(df$Group == "B"), 60, 5)
df$Age[df$Group == "C"] <- rnorm(sum(df$Group == "C"), 31, 5)
df$Status <- ifelse(df$Group == "B",
                     sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.95, 0.05)),
                     sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.05, 0.95)))

tbl <- run_summarytable(
  df,
  split_by                  = "Group",
  add_inferential_pvalues   = TRUE,
  mark_posthoc_significance = TRUE,
  test_type_continuous      = "parametric",
  force_categorical         = "Status"
)
dat <- tbl[["_data"]]

# Age: significant ANOVA (A/B/C differ, B is the outlier) -> B's column
# (stat_2) must carry a post-hoc mark; A and C (similar means) are not
# guaranteed to differ from each other, so only check B's involvement.
age_row <- dat[dat$variable == "Age" & dat$row_type == "label", ]
stopifnot(grepl("#[0-9,]+#$", age_row$stat_2))

# Status: at least one level cell should be flagged, since Group B was
# simulated with a very different Yes/No split than A and C.
status_rows <- dat[dat$variable == "Status" & dat$row_type == "level", ]
stopifnot(any(grepl("#[0-9]+#$", unlist(status_rows[, c("stat_1", "stat_2", "stat_3")]))))

# Footnote block exists and is keyed per-variable.
note <- tbl$`_source_notes`[[1]]
stopifnot(grepl("Post-hoc significance", note), grepl("<b>Age</b>", note), grepl("<b>Status</b>", note))

# Two-group split (t-test): no post-hoc for the continuous variable, but the
# categorical adjusted-residual post-hoc still applies.
df2 <- df[df$Group != "C", ]
tbl2 <- run_summarytable(df2, split_by = "Group", add_inferential_pvalues = TRUE,
                          mark_posthoc_significance = TRUE, force_categorical = "Status")
dat2 <- tbl2[["_data"]]
age_row2 <- dat2[dat2$variable == "Age" & dat2$row_type == "label", ]
stopifnot(!grepl("#", age_row2$stat_1), !grepl("#", age_row2$stat_2))
status_rows2 <- dat2[dat2$variable == "Status" & dat2$row_type == "level", ]
stopifnot(any(grepl("#[0-9]+#$", unlist(status_rows2[, c("stat_1", "stat_2")]))))

# mark_posthoc_significance is a no-op (with a warning) when
# add_inferential_pvalues = FALSE.
w <- tryCatch({
  run_summarytable(df, split_by = "Group", mark_posthoc_significance = TRUE)
  NULL
}, warning = function(w) w)
stopifnot(!is.null(w), grepl("no effect", conditionMessage(w)))

# Regression check: binary Yes/No variables that gtsummary collapses into a
# single "dichotomous" label row (i.e. NOT force_categorical'd) must still
# get marked -- this was the reported bug.
tbl_dichot <- run_summarytable(df, split_by = "Group", add_inferential_pvalues = TRUE,
                                mark_posthoc_significance = TRUE)
dat_dichot <- tbl_dichot[["_data"]]
status_label_row <- dat_dichot[dat_dichot$variable == "Status" & dat_dichot$row_type == "label", ]
stopifnot(nrow(status_label_row) == 1L)
stopifnot(any(grepl("#[0-9]+#$", unlist(status_label_row[, c("stat_1", "stat_2", "stat_3")]))))

# mark_posthoc_significance now defaults to TRUE.
stopifnot(isTRUE(formals(run_summarytable)$mark_posthoc_significance))

cat("All run_summarytable post-hoc checks passed.\n")
