# Smoke check for plot_fcmetric(); run after devtools::load_all()
set.seed(4821)
m    <- as.data.frame(matrix(rexp(60 * 12, 0.01) + 1, 60, 12, dimnames = list(NULL, paste0("m", 1:12))))
meta <- data.frame(g = factor(rep(c("A", "B"), 30)))
fc   <- run_foldchange(m, metadata = meta, group = "g")
stopifnot(!"regulation" %in% names(fc$summary_table))

ok <- function(r) inherits(r$plots[[1]], "ggplot")
stopifnot(ok(plot_fcmetric(fc, annotate = c("m1", "m2"), x_cut = 100)))
stopifnot(ok(plot_fcmetric(fc, metric = "sd")), ok(plot_fcmetric(fc, metric = "variance", log2 = FALSE)))

# pooled SD matches a direct computation
sd_direct <- sqrt(((29 * sd(m$m1[meta$g == "A"])^2) + (29 * sd(m$m1[meta$g == "B"])^2)) / 58)
stopifnot(isTRUE(all.equal(plot_fcmetric(fc, metric = "sd")$classified$metric_value[
  fc$summary_table$feature == "m1"], sd_direct)))

diff <- run_diff(m, metadata = meta, outcome = colnames(m), group = "g",
                 summary_table = TRUE, verbose = FALSE)
stopifnot(ok(plot_fcmetric(list(fc = fc, diff = diff), metric = "effect_size")))

# VIP with a stand-in run_pls object
pls <- structure(list(vip_scores = setNames(runif(12, 0, 2), colnames(m))), class = c("run_pls", "list"))
stopifnot(ok(plot_fcmetric(list(fc = fc, pls = pls), metric = "vip", x_cut = 1)))

# missing companion object / bad metric must error
stopifnot(inherits(try(plot_fcmetric(fc, metric = "vip"), silent = TRUE), "try-error"))
stopifnot(inherits(try(plot_fcmetric(fc, metric = "nope"), silent = TRUE), "try-error"))

# data.frame input and group split
d <- data.frame(fc = 2^rnorm(40), ab = rexp(40, 0.01) + 1, g = rep(c("A", "B"), 20))
stopifnot(length(plot_fcmetric(d, "fc", "ab", group = "g")$plots) == 2)
