make_df <- function() {
  set.seed(4817)
  n <- 400
  g <- rep(c("A", "B"), each = n / 2)
  data.frame(
    grp  = g,
    bin  = c(rbinom(n / 2, 1, 0.3), rbinom(n / 2, 1, 0.7)) |> factor(labels = c("No", "Yes")),
    kcat = sample(c("k1", "k2", "k3", "k4"), n, replace = TRUE)
  )
}

notes <- function(tbl) tbl[["_source_notes"]][[1]]

test_that("2x2: no cell-level post hoc, direction line in footnote", {
  skip_if_not_installed("gt")
  out <- run_summarytable(make_df(), summarize_what = "bin", split_by = "grp",
                          add_inferential_pvalues = TRUE)
  expect_false(grepl("Post-hoc significance", notes(out)))
  expect_match(notes(out), "Yes more common in B", fixed = TRUE)
})

test_that("k x 2: adjusts over one column (k tests, not 2k)", {
  ct <- matrix(c(32, 65, 119, 4, 18, 65, 147, 1), ncol = 2)
  p  <- 2 * pnorm(-abs(suppressWarnings(chisq.test(ct, correct = FALSE))$stdres))
  expect_equal(p[, 1], p[, 2])
  expect_lt(min(p.adjust(p[, 1], "bonferroni")), min(p.adjust(as.vector(p), "bonferroni")))
})

test_that("significant omnibus but no significant cell gets a footnote", {
  skip_if_not_installed("gt")
  ct <- matrix(c(32, 65, 119, 4, 18, 65, 147, 1), ncol = 2)
  d <- data.frame(
    grp = rep(c("A", "B"), colSums(ct)),
    age = c(rep(c("m", "y", "o", "x"), ct[, 1]), rep(c("m", "y", "o", "x"), ct[, 2]))
  )
  out <- run_summarytable(d, summarize_what = "age", split_by = "grp",
                          add_inferential_pvalues = TRUE, p_adjust_method = "bonferroni")
  expect_match(notes(out), "no pairwise comparison/cell remained significant")
})
