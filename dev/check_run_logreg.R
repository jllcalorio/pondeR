# Minimal self-check for firth_control / firth_plcontrol pass-through in run_logreg().
# Run with: source("dev/check_run_logreg_firth_control.R")

devtools::load_all(".", quiet = TRUE)

set.seed(4821)
dat <- data.frame(
  # complete separation -> a genuine zero cell -> Firth kicks in automatically
  outcome = factor(c(rep("No", 20), rep("Yes", 20))),
  x       = factor(c(rep("a", 20), rep("b", 20)))
)

# Sanity: this combination has a zero cell, so Firth kicks in automatically.
default_fit <- run_logreg(dat, y = "outcome", indep = "x", ref = "No", apply_firth = TRUE)
stopifnot(isTRUE(attr(default_fit, "model_metrics")$`Firth Corrected`))

# A control object with an absurdly low maxit should still be accepted and used
# (result differs from the default-control fit, proving it was actually passed through).
tight_ctrl <- logistf::logistf.control(maxit = 1)
tight_fit <- run_logreg(
  dat, y = "outcome", indep = "x", ref = "No",
  apply_firth = TRUE, firth_control = tight_ctrl
)

stopifnot(isTRUE(attr(tight_fit, "model_metrics")$`Firth Corrected`))
stopifnot(!isTRUE(all.equal(
  default_fit$`Log Odds`[default_fit$Predictor == "xb"],
  tight_fit$`Log Odds`[tight_fit$Predictor == "xb"]
)))

# plcontrol should also be accepted without error.
loose_plctrl <- logistf::logistf.control(maxit = 200)
invisible(run_logreg(
  dat, y = "outcome", indep = "x", ref = "No",
  apply_firth = TRUE, firth_plcontrol = loose_plctrl
))

# has_zero_cells() must also catch numeric 0/1-coded indicator columns, not
# just factor/character predictors -- otherwise dynamic Firth detection
# silently misses separated numeric binary predictors.
dat_num <- data.frame(
  outcome = factor(c(rep("No", 20), rep("Yes", 20))),
  x       = c(rep(0, 20), rep(1, 20))
)
numeric_fit <- run_logreg(dat_num, y = "outcome", indep = "x", ref = "No", apply_firth = TRUE)
stopifnot(isTRUE(attr(numeric_fit, "model_metrics")$`Firth Corrected`))

# stepwise = TRUE should drop the higher-VIF member of a near-collinear pair
# and leave an uncorrelated variable untouched.
set.seed(7734)
n <- 200
z <- rnorm(n)
dat_vif <- data.frame(
  outcome = factor(sample(c("No", "Yes"), n, replace = TRUE)),
  a = z + rnorm(n, sd = 0.05),
  b = z + rnorm(n, sd = 0.05),
  c = rnorm(n)
)
step_fit <- run_logreg(
  dat_vif, y = "outcome", indep = c("a", "b", "c"), ref = "No",
  stepwise = TRUE, step_n = 1, apply_firth = FALSE
)
removed <- attr(step_fit, "stepwise_removed")
stopifnot(nrow(removed) == 1, removed$Removed == "b", removed$VIF > 5)
stopifnot(setequal(step_fit$Predictor, c("(Intercept)", "a", "c")))
stopifnot(all(step_fit$VIF[step_fit$Predictor %in% c("a", "c")] < 5, na.rm = TRUE))

# stepwise = TRUE should be a no-op when no VIF reaches the threshold.
no_collinearity_fit <- run_logreg(
  dat_vif, y = "outcome", indep = c("a", "c"), ref = "No",
  stepwise = TRUE, step_n = 1, apply_firth = FALSE
)
stopifnot(nrow(attr(no_collinearity_fit, "stepwise_removed")) == 0)

# stepwise + iterate is an invalid combination.
stopifnot(inherits(
  tryCatch(
    run_logreg(dat_vif, y = "outcome", indep = c("a", "b"), ref = "No",
               stepwise = TRUE, iterate = TRUE),
    error = function(e) e
  ),
  "error"
))

message("run_logreg stepwise VIF removal: OK")

message("run_logreg firth_control / firth_plcontrol pass-through: OK")
