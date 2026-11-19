# Inferential and descriptive statistics, checked against base R

iv <- droplevels(iris[iris$Species != "setosa", ])

test_that("run_diff reproduces a two-sample t-test and picks Kruskal/ANOVA for 3 groups", {
  r <- suppressMessages(run_diff(iv, outcome = "Sepal.Length", group = "Species"))
  expect_s3_class(r, "run_diff")
  ref <- t.test(Sepal.Length ~ Species, data = iv, var.equal = grepl("Student", r$test_used))
  expect_equal(unname(r$test_result$p.value), ref$p.value, tolerance = 1e-8)
  expect_equal(r$effect_size$metric, "Cohen's d")

  r3 <- suppressMessages(run_diff(iris, outcome = "Sepal.Length", group = "Species"))
  expect_equal(r3$n_groups, 3)
  expect_lt(r3$test_result$p.value, 1e-10)

  expect_error(run_diff(iv, outcome = "nope", group = "Species"))
})

test_that("run_diff multi-outcome summary_table and its gt rendering", {
  r <- run_diff(iv[, 1:4], iv["Species"], names(iv)[1:4], "Species",
                summary_table = TRUE, verbose = FALSE)
  st <- r$summary_table
  expect_s3_class(st, "run_diff_summary")
  expect_equal(nrow(st), 4)
  expect_true(all(st$p_value_adj >= st$p_value))
  expect_s3_class(st[1:2, ], "run_diff_summary")        # `[` keeps the class
  skip_if_not_installed("gt")
  expect_s3_class(as_gt(st), "gt_tbl")
})

test_that("run_diff handles repeated measures", {
  set.seed(1)
  d <- data.frame(ID = factor(rep(1:12, each = 3)),
                  Time = factor(rep(c("T1", "T2", "T3"), 12)))
  d$HR <- 70 + as.integer(d$Time) * 3 + rnorm(36)
  r <- suppressWarnings(suppressMessages(
    run_diff(d, outcome = "HR", within = "Time", subject_id = "ID")))
  expect_s3_class(r, "run_diff")
  expect_lt(r$test_result$p.value, 0.001)
})

test_that("run_assoc matches chisq.test and run_single_assoc tests goodness of fit", {
  a <- suppressMessages(run_assoc(mtcars, "am", "vs", force_categorical = TRUE))
  expect_s3_class(a, "run_assoc")
  ref <- suppressWarnings(chisq.test(table(mtcars$am, mtcars$vs)))
  expect_equal(a$test_result$p.value, ref$p.value, tolerance = 1e-8)

  s <- run_single_assoc(iris, y = "Species")
  expect_equal(s$summary_table$p_value, 1)              # perfectly balanced
})

test_that("run_correl matches cor() and as_gt renders it", {
  c1 <- run_correl(mtcars[, c("mpg", "wt", "hp")], method = "pearson")
  expect_s3_class(c1, "run_correl")
  r <- c1$correlations
  expect_equal(r$estimate[r$var1 == "mpg" & r$var2 == "wt"],
               cor(mtcars$mpg, mtcars$wt), tolerance = 1e-10)
  expect_equal(nrow(run_correl(mtcars[, 1:4], y = "mpg", method = "pearson")$correlations), 3)
  skip_if_not_installed("gt")
  expect_s3_class(as_gt(c1), "gt_tbl")
  expect_s3_class(as_gt(head(mtcars)), "gt_tbl")
})

test_that("as_gt forwards gtsummary objects to gtsummary::as_gt", {
  skip_if_not_installed("gt")
  tb <- gtsummary::tbl_summary(iris[, c("Sepal.Length", "Species")])
  expect_s3_class(as_gt(tb), "gt_tbl")
})

test_that("run_foldchange computes ratios of group means", {
  fc <- run_foldchange(iv[, 1:4], metadata = iv["Species"], group = "Species")
  st <- fc$summary_table
  row <- st[st$feature == "Petal.Length", ]
  m <- tapply(iv$Petal.Length, iv$Species, mean)
  expect_equal(row$fold_change, unname(m[row$num_group] / m[row$den_group]))
  expect_equal(row$log2_fc, log2(row$fold_change))
})

test_that("run_auc gives AUC 1 for a perfectly separating marker", {
  d <- data.frame(g = rep(0:1, each = 20), a = c(1:20, 31:50), b = rnorm(40))
  r <- run_auc(d, y = "g", plot = "none", compute_ci = FALSE)
  expect_s3_class(r, "run_auc")
  expect_equal(r$auc_table$auc[r$auc_table$predictor == "a"], 1)
  expect_equal(r$n_cases + r$n_controls, 40)
})

test_that("run_logreg matches glm() coefficients", {
  d <- mtcars
  d$am <- factor(d$am, labels = c("auto", "manual"))
  r <- run_logreg(d, y = "am", indep = "wt", ref = "auto", crude = TRUE,
                  apply_firth = FALSE)
  ref <- coef(glm(am ~ wt, family = binomial, data = d))
  expect_equal(r$`Log Odds`[r$Predictor == "wt"], unname(ref["wt"]), tolerance = 1e-6)
  expect_equal(r$OR[r$Predictor == "wt"], exp(unname(ref["wt"])), tolerance = 1e-6)
})

test_that("run_summarytable returns a gt table", {
  skip_if_not_installed("gt")
  expect_s3_class(run_summarytable(iris, split_by = "Species"), "gt_tbl")
})

test_that("run_contentvaliditycvr follows Lawshe's formula", {
  skip_if_not_installed("contentValidity")
  r <- run_contentvaliditycvr(data.frame(i1 = c(1, 1, 1, 2, 1), i2 = c(3, 3, 1, 3, 3)),
                              essential = 1)
  expect_equal(r$CVR, c((4 - 2.5) / 2.5, (1 - 2.5) / 2.5))
})

test_that("run_contentvalidity computes I-CVI", {
  skip_if_not_installed("contentValidity")
  skip_if_not_installed("irrCAC")
  r <- run_contentvalidity(data.frame(i1 = c(4, 3, 4, 4, 3), i2 = c(2, 3, 2, 1, 2)))
  expect_equal(r$items$`I-CVI`, c(1, 0.2))             # share of ratings >= 3
  expect_equal(unname(r$scale["S-CVI/Ave"]), 0.6)
})

test_that("run_anthroindex returns WHO z-scores", {
  skip_if_not_installed("anthro")
  d <- data.frame(sex = c("Male", "Female"), age = c(12, 24),
                  wt = c(9.6, 11.5), ht = c(75.7, 86.4))
  r <- run_anthroindex(d, sex_col = "sex", age_col = "age", age_unit = "month",
                       weight_col = "wt", height_col = "ht",
                       height_measured = "L", return = "wfa")
  expect_true(all(abs(r$z_score) < 0.5))                # WHO medians -> z near 0
  expect_true("Nutritional Status" %in% names(r))
})

test_that("run_mixedmodel recovers a fixed effect", {
  set.seed(42)
  meta <- data.frame(subject_id = factor(rep(1:15, each = 2)),
                     time = rep(c("pre", "post"), 15))
  x <- data.frame(g = rnorm(30, ifelse(meta$time == "post", 6, 5), 0.3),
                  i = rnorm(30, 80, 10))
  m <- suppressMessages(run_mixedmodel(x, meta, verbose = FALSE,
                        formula = ".feature ~ time + (1 | subject_id)"))
  fe <- m$combined_fixed_effects
  est <- fe$estimate[fe$feature == "g" & fe$term == "timepre"]
  expect_equal(est, -1, tolerance = 0.25)
  expect_lt(fe$p_value[fe$feature == "g" & fe$term == "timepre"], 0.001)
})

test_that("run_lefse ranks the planted biomarker first, including a single feature", {
  set.seed(1)
  x <- as.data.frame(matrix(rpois(40 * 5, 20), 40,
                            dimnames = list(NULL, paste0("t", 1:5))))
  meta <- data.frame(Group = rep(c("A", "B"), each = 20))
  x$t1[21:40] <- x$t1[21:40] + 40
  r <- suppressWarnings(run_lefse(x / rowSums(x) * 1e6, meta, verbose = FALSE))
  expect_equal(r$summary_table$feature[1], "t1")
  expect_equal(nrow(r$lda), length(unique(r$lda$feature)))   # one row per feature

  one <- suppressWarnings(run_lefse(data.frame(a = c(1:10, 31:40)),
                                    data.frame(Group = rep(c("A", "B"), each = 10)),
                                    verbose = FALSE))
  expect_equal(one$lda$feature, "a")
})

test_that("get_volcanodata merges fold changes with p-values", {
  fc  <- run_foldchange(iv[, 1:4], metadata = iv["Species"], group = "Species")
  dif <- run_diff(iv[, 1:4], iv["Species"], names(iv)[1:4], "Species",
                  summary_table = TRUE, verbose = FALSE)
  v <- get_volcanodata(fc, dif, filter = FALSE)
  expect_equal(nrow(v), 4)
  expect_error(get_volcanodata(list()), "At least one")
})

test_that("get_groupsizes and get_removed_samples account for missing covariates", {
  meta <- data.frame(G = rep(c("a", "b"), 30), A = c(NA, 2:60))
  gs <- get_groupsizes(meta, pred = "G", not_penalized = "A", cv_folds = 5L)
  expect_equal(sum(gs$n_total), 59)

  rm <- get_removed_samples(data.frame(f = 1:4),
                            data.frame(G = c("a", "b", "a", "b"), A = c(NA, 1, 2, 3)),
                            not_penalized = "A", pred = "G")
  expect_equal(rm$.row, 1)
})

test_that("run_regreg and get_sigfeatures find the informative feature", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  skip_if_not_installed("caret")
  set.seed(3)
  n <- 80
  feat <- as.data.frame(matrix(rnorm(n * 10), n,
                               dimnames = list(NULL, paste0("f", 1:10))))
  meta <- data.frame(Score = 3 * feat$f1 + rnorm(n, sd = 0.3))
  res <- suppressWarnings(suppressMessages(
    run_regreg(feat, metadata = meta, pred = "Score", alpha = 1, seed = 1)))
  sf <- get_sigfeatures(res)
  expect_true("f1" %in% unlist(sf))
})
