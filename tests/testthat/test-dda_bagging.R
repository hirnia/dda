## ============================================================================
## testthat file: test-dda_bagging.R
## Tests for dda.bagging(), dda.decisions() and agg.value()
## ============================================================================

# --- Shared setup -----------------------------------------------------------
generate_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- rchisq(n, df = 4) - 4
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e
  data.frame(x = x, y = y)
}

generate_test_data_cov <- function(n = 200, seed = 42) {
  set.seed(seed)
  z <- rnorm(n)
  x <- rchisq(n, df = 4) - 4 + 0.3 * z
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e + 0.2 * z
  data.frame(x = x, y = y, z = z)
}

d     <- generate_test_data()
d_cov <- generate_test_data_cov()

base_indep      <- dda.indep(y ~ x, pred = "x", data = d, B = 10)
base_indep_full <- dda.indep(y ~ x, pred = "x", data = d, B = 10,
                             hetero = TRUE, nlfun = 2, diff = TRUE)
base_indep_cov  <- dda.indep(y ~ x + z, pred = "x", data = d_cov, B = 10)
base_resdist    <- dda.resdist(y ~ x, pred = "x", data = d, B = 10)
base_resdist_pt <- dda.resdist(y ~ x, pred = "x", data = d, B = 10, prob.trans = TRUE)
base_vardist    <- dda.vardist(y ~ x, pred = "x", data = d, B = 10)

run_bag <- function(base, dat, iter = 10, ...) {
  dda.bagging(base, data = dat, iter = iter, progress = FALSE, ...)
}

bag_indep   <- run_bag(base_indep_full, d)
bag_resdist <- run_bag(base_resdist, d)
bag_vardist <- run_bag(base_vardist, d)

## ============================================================================
## 1. Input validation
## ============================================================================

test_that("dda.bagging rejects objects that are not DDA results", {
  expect_error(dda.bagging(list(a = 1), data = d, iter = 5),
               regexp = "must be a dda.indep, dda.resdist, or dda.vardist")
})

test_that("dda.bagging requires a data.frame", {
  expect_error(dda.bagging(base_indep, iter = 5), regexp = "Please provide")
  expect_error(dda.bagging(base_indep, data = NULL, iter = 5), regexp = "Please provide")
})

test_that("dda.bagging rejects an unknown agg_stat", {
  expect_error(dda.bagging(base_indep, data = d, iter = 5, agg_stat = "geometric"),
               regexp = "should be one of")
})

## ============================================================================
## 2. Return structure
## ============================================================================

test_that("classes follow the base DDA object", {
  expect_s3_class(bag_indep,   c("dda_bagging_indep", "dda_bagging"), exact = TRUE)
  expect_s3_class(bag_resdist, c("dda_bagging_resdist", "dda_bagging"), exact = TRUE)
  expect_s3_class(bag_vardist, c("dda_bagging_vardist", "dda_bagging"), exact = TRUE)
})

test_that("output holds all top-level elements", {
  expect_named(bag_indep,
               c("bagged_results", "aggregated_stats", "decisions",
                 "decision_proportions", "ols", "n_valid_iterations",
                 "alpha", "agg_stat", "trim_prob", "win_prob"),
               ignore.order = TRUE)
})

test_that("settings are stored in the object", {
  res <- run_bag(base_indep, d, iter = 5, alpha = 0.10,
                 agg_stat = "winsorized", trim_prob = 0.2, win_prob = 0.05)
  expect_equal(res$alpha, 0.10)
  expect_equal(res$agg_stat, "winsorized")
  expect_equal(res$trim_prob, 0.2)
  expect_equal(res$win_prob, 0.05)
})

test_that("n_valid_iterations matches the stored bootstrap results", {
  expect_true(bag_indep$n_valid_iterations > 0)
  expect_true(bag_indep$n_valid_iterations <= 10)
  expect_length(bag_indep$bagged_results, bag_indep$n_valid_iterations)
  expect_equal(nrow(bag_indep$decisions), bag_indep$n_valid_iterations)
})

test_that("agg_stat is stored for every method", {
  for (method in c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey")) {
    res <- run_bag(base_indep, d, iter = 5, agg_stat = method)
    expect_equal(res$agg_stat, method, info = method)
  }
})

## ============================================================================
## 3. Aggregated statistics
## ============================================================================

test_that("aggregated_stats keeps the class of the base object", {
  expect_s3_class(bag_indep$aggregated_stats, "dda.indep")
  expect_s3_class(bag_resdist$aggregated_stats, "dda.resdist")
  expect_s3_class(bag_vardist$aggregated_stats, "dda.vardist")
  expect_equal(bag_indep$aggregated_stats$var.names, c("y", "x"))
})

test_that("indep aggregates hold HSIC, dCor, BP, nlcor and difference results", {
  agg <- bag_indep$aggregated_stats
  expect_true(is.numeric(agg$hsic.yx$statistic))
  expect_true(is.numeric(agg$hsic.xy$p.value))
  expect_true(is.numeric(agg$distance_cor.dcor_yx$statistic))
  expect_true(is.numeric(agg$distance_cor.dcor_xy$p.value))
  expect_length(agg$breusch_pagan, 4)
  expect_false(is.null(agg$nlcor.yx$t1))
  expect_false(is.null(agg$nlcor.xy$t3))
  expect_equal(dim(agg$out.diff), dim(base_indep_full$out.diff))
})

test_that("resdist and vardist aggregates hold the separate tests", {
  expect_true(is.numeric(bag_resdist$aggregated_stats$agostino$target$statistic))
  expect_true(is.numeric(bag_resdist$aggregated_stats$anscombe$alternative$p.value))
  expect_true(is.numeric(bag_vardist$aggregated_stats$agostino$predictor$statistic))
  expect_true(is.numeric(bag_vardist$aggregated_stats$anscombe$outcome$p.value))
})

test_that("resdist boot.warning follows the aggregated kurtosis signs", {
  agg <- bag_resdist$aggregated_stats
  expect_equal(agg$boot.warning,
               sign(agg$anscombe$alternative$statistic[1]) != sign(agg$anscombe$target$statistic[1]))
})

test_that("trimmed with trim_prob = 0 and winsorized with win_prob = 0 equal the mean", {
  set.seed(101); res_mean <- run_bag(base_indep, d, iter = 5)
  set.seed(101); res_trim <- run_bag(base_indep, d, iter = 5, agg_stat = "trimmed", trim_prob = 0)
  set.seed(101); res_win  <- run_bag(base_indep, d, iter = 5, agg_stat = "winsorized", win_prob = 0)
  expect_equal(res_mean$aggregated_stats$hsic.yx$statistic,
               res_trim$aggregated_stats$hsic.yx$statistic)
  expect_equal(res_mean$aggregated_stats$hsic.yx$statistic,
               res_win$aggregated_stats$hsic.yx$statistic)
})

test_that("agg.value drops non-finite values and handles every method", {
  x <- c(1, 2, 3, 4, 100, NA, Inf)
  expect_equal(agg.value(x, "mean"), mean(c(1, 2, 3, 4, 100)))
  expect_equal(agg.value(x, "median"), 3)
  expect_equal(agg.value(x, "midhinge"), mean(quantile(c(1, 2, 3, 4, 100), c(0.25, 0.75), names = FALSE)))
  expect_true(is.numeric(agg.value(x, "tukey")))
  expect_true(is.na(agg.value(c(NA, Inf), "mean")))
  expect_error(agg.value(x, "geometric"), regexp = "Unknown agg_stat")
})

## ============================================================================
## 4. Decisions
## ============================================================================

test_that("decision_proportions has one row per test and rows sum to 1", {
  props <- bag_indep$decision_proportions
  expect_equal(rownames(props), colnames(bag_indep$decisions))
  expect_equal(colnames(props), c("Target", "Alternative", "Confounding", "Undecided"))
  expect_equal(unname(rowSums(props)), rep(1, nrow(props)), tolerance = 1e-8)
})

test_that("resdist and vardist decisions have no Confounding level", {
  expect_equal(colnames(bag_resdist$decision_proportions), c("Target", "Alternative", "Undecided"))
  expect_equal(colnames(bag_vardist$decision_proportions), c("Target", "Alternative", "Undecided"))
})

test_that("dda.decisions returns the expected tests for each object type", {
  expect_setequal(names(dda.decisions(base_indep_full)),
                  c("hsic", "dcor", "bp", "nlcor", "hsic.diff", "dcor.diff", "mi.diff"))
  expect_true(all(c("agostino", "anscombe", "skewdiff", "kurtdiff") %in%
                    names(dda.decisions(base_vardist))))
  expect_true(all(dda.decisions(base_resdist) %in% c("Target", "Alternative", "Undecided", NA)))
})

test_that("dda.decisions swaps the separate resdist tests under prob.trans = TRUE", {
  obj <- base_resdist_pt
  obj$agostino$target$p.value      <- 0.01
  obj$agostino$alternative$p.value <- 0.50
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Target")
  obj$probtrans <- FALSE
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Alternative")
})

test_that("both significant separate tests are confounding for indep and undecided otherwise", {
  obj <- base_indep
  obj$hsic.yx$p.value <- 0.01
  obj$hsic.xy$p.value <- 0.01
  expect_equal(unname(dda.decisions(obj)["hsic"]), "Confounding")
  obj <- base_vardist
  obj$agostino$outcome$p.value   <- 0.01
  obj$agostino$predictor$p.value <- 0.01
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Undecided")
})

test_that("dda.decisions rejects other objects", {
  expect_error(dda.decisions(list(a = 1)), regexp = "must be a dda.indep")
})

## ============================================================================
## 5. save_file, covariates, prob.trans
## ============================================================================

test_that("save_file writes the dda_bagging object to disk", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  run_bag(base_indep, d, iter = 5, save_file = tmp)
  saved <- readRDS(tmp)
  expect_s3_class(saved, "dda_bagging_indep")
})

test_that("dda.bagging works with covariates", {
  res <- run_bag(base_indep_cov, d_cov, iter = 5)
  expect_s3_class(res, "dda_bagging_indep")
  expect_true("z" %in% colnames(res$ols$target$coef))
  expect_true("z" %in% colnames(res$ols$alternative$coef))
})

test_that("dda.bagging handles prob.trans = TRUE", {
  res <- run_bag(base_resdist_pt, d, iter = 5)
  expect_true(isTRUE(res$aggregated_stats$probtrans))
  expect_true("agostino" %in% rownames(res$decision_proportions))
})

## ============================================================================
## 6. Failed bootstrap samples
## ============================================================================

test_that("failed bootstrap samples are dropped with a warning", {
  n.calls <- 0
  B.fail <- function() {
    n.calls <<- n.calls + 1
    if (n.calls %% 2 == 0) stop("planned failure")
    10
  }
  base <- dda.vardist(y ~ x, pred = "x", data = d, B = B.fail())
  expect_warning(res <- dda.bagging(base, data = d, iter = 4, progress = FALSE),
                 regexp = "2 of 4 bootstrap samples failed and were dropped. Last error: planned failure")
  expect_equal(res$n_valid_iterations, 2)
  expect_equal(nrow(res$decisions), 2)
  expect_equal(nrow(res$ols$target$coef), 2)
})

test_that("dda.bagging reports the error when every bootstrap sample fails", {
  expect_error(dda.bagging(base_vardist, data = data.frame(y = d$y, w = d$x), iter = 2, progress = FALSE),
               regexp = "failed in every bootstrap sample: .*not found")
})
