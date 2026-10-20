# Self-check: plot_volcano with character p-values from run_diff(simplify = TRUE)
# Run with: Rscript dev/check_plot_volcano_characs.R (from package root)
devtools::load_all(".", quiet = TRUE)

success_count <- 0
fail_count <- 0

run_test <- function(name, expr) {
  tryCatch({
    expr
    cat("✓", name, "\n")
    success_count <<- success_count + 1
  }, error = function(e) {
    cat("✗", name, "—", conditionMessage(e), "\n")
    fail_count <<- fail_count + 1
  })
}

cat("=== plot_volcano Character p-value Tests ===\n\n")

cat("Test 1: run_diff output → plot_volcano\n")
set.seed(9001)
n <- 200
df <- data.frame(
  outcome1 = rnorm(n, mean = 0, sd = 1),
  outcome2 = rnorm(n, mean = 0.3, sd = 1),
  group    = rep(c("A", "B"), c(n/2, n/2))
)

diff_res <- run_diff(
  x            = df,
  outcome      = c("outcome1", "outcome2"),
  group        = "group",
  test_type    = "parametric",
  summary_table = TRUE,
  simplify     = TRUE
)

run_test("  p_value is character", {
  stopifnot(is.character(diff_res$summary_table$p_value))
})

run_test("  plot_volcano accepts run_diff output", {
  res <- plot_volcano(x = list(diff = diff_res), y = "fold_change")
  stopifnot(!is.null(res$classified))
  stopifnot(is.numeric(res$classified$p_value_used))
})

cat("\nTest 2: Direct data frame input\n")
volcano_df <- data.frame(
  feature = c("Gene_A", "Gene_B", "Gene_C"),
  fold_change = c(2.5, 0.4, 1.8),
  p_value = c(".012ᵃ", ".456ᵇ", ".089ᶜ")
)

run_test("  plot_volcano converts character p_value to numeric", {
  res <- plot_volcano(x = volcano_df, y = "fold_change", z = "p_value")
  stopifnot(is.numeric(res$classified$p_value))
})

cat("\nTest 3: p_value_adj column\n")
volcano_df2 <- data.frame(
  feature = c("Gene_X", "Gene_Y"),
  fold_change = c(1.5, 0.6),
  p_value_adj = c(".033ᵃ", ".201ᵃ")
)

run_test("  plot_volcano converts character p_value_adj to numeric", {
  res <- plot_volcano(x = volcano_df2, y = "fold_change", z = "p_value_adj")
  stopifnot(is.numeric(res$classified$p_value_adj))
})

cat("\n", success_count, "passed,", fail_count, "failed.\n")
if (fail_count > 0) quit(status = 1)
