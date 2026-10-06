## ============================================================================
## testthat file: test-print_dda_bagging.R
## Tests for print.dda_bagging(), reaggregate_bagging() and summary_ols()
## ============================================================================

# --- Shared setup -----------------------------------------------------------
generate_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- rchisq(n, df = 4) - 4
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e
  data.frame(x = x, y = y)
}

d <- generate_test_data()

bag_indep <- dda.bagging(
  dda.indep(y ~ x, pred = "x", data = d, B = 10,
            hetero = TRUE, nlfun = 2, diff = TRUE),
  data = d, iter = 10, progress = FALSE
)
bag_resdist <- dda.bagging(
  dda.resdist(y ~ x, pred = "x", data = d, B = 10),
  data = d, iter = 10, progress = FALSE
)
bag_vardist <- dda.bagging(
  dda.vardist(y ~ x, pred = "x", data = d, B = 10),
  data = d, iter = 10, progress = FALSE, agg_stat = "trimmed"
)

## ============================================================================
## 1. reaggregate_bagging()
## ============================================================================

test_that("reaggregate_bagging returns the object unchanged when agg_stat is NULL", {
  expect_identical(reaggregate_bagging(bag_indep, agg_stat = NULL), bag_indep)
})

test_that("reaggregate_bagging stores the new method", {
  res <- reaggregate_bagging(bag_indep, agg_stat = "median")
  expect_equal(res$agg_stat, "median")
  expect_true(is.numeric(res$aggregated_stats$hsic.yx$statistic))
})

test_that("a new trim_prob alone changes the trimmed aggregate", {
  res <- reaggregate_bagging(bag_vardist, agg_stat = "trimmed", trim_prob = 0.4)
  expect_equal(res$trim_prob, 0.4)
  expected <- agg.value(sapply(bag_vardist$bagged_results, function(b) b$agostino$outcome$statistic[1]),
                        "trimmed", trim_prob = 0.4)
  expect_equal(unname(res$aggregated_stats$agostino$outcome$statistic[1]), expected)
})

test_that("stored trim_prob and win_prob are the defaults when re-aggregating", {
  res <- reaggregate_bagging(bag_vardist, agg_stat = "trimmed")
  expect_equal(res$trim_prob, bag_vardist$trim_prob)
  expect_equal(res$aggregated_stats, bag_vardist$aggregated_stats)
})

test_that("reaggregate_bagging works for resdist objects", {
  res <- reaggregate_bagging(bag_resdist, agg_stat = "tukey")
  expect_equal(res$agg_stat, "tukey")
  expect_s3_class(res$aggregated_stats, "dda.resdist")
})

## ============================================================================
## 2. print.dda_bagging()
## ============================================================================

test_that("print shows the header, sample count and method", {
  expect_output(print(bag_indep), regexp = "BOOTSTRAP AGGREGATED DDA")
  expect_output(print(bag_indep), regexp = "Number of bootstrap samples")
  expect_output(print(bag_indep), regexp = "Aggregation method: mean")
})

test_that("print shows the tests of each object type", {
  expect_output(print(bag_indep), regexp = "HSIC")
  expect_output(print(bag_indep), regexp = "Robust BP-test")
  expect_output(print(bag_indep), regexp = "Non-linear Correlation")
  expect_output(print(bag_resdist), regexp = "Residual Distributions")
  expect_output(print(bag_vardist), regexp = "Variable Distributions")
})

test_that("print respects an agg_stat override", {
  expect_output(print(bag_indep, agg_stat = "median"), regexp = "Aggregation method: median")
})

test_that("print returns the original object invisibly", {
  out <- withVisible(print(bag_indep, agg_stat = "median"))
  expect_false(out$visible)
  expect_identical(out$value, bag_indep)
})

## ============================================================================
## 3. summary_ols()
## ============================================================================

test_that("summary_ols rejects other objects", {
  expect_error(summary_ols(list(a = 1)), regexp = "must be a dda_bagging object")
})

test_that("summary_ols prints both models with R-squared", {
  expect_output(summary_ols(bag_indep), regexp = "OLS Summary: Target Model")
  expect_output(summary_ols(bag_indep), regexp = "OLS Summary: Alternative Model")
  expect_output(summary_ols(bag_indep), regexp = "R-squared")
})

test_that("summary_ols respects an agg_stat override", {
  expect_output(summary_ols(bag_indep, agg_stat = "median"), regexp = "Aggregation method: median")
})

test_that("summary_ols returns coefficient tables invisibly", {
  out <- withVisible(summary_ols(bag_resdist))
  expect_false(out$visible)
  expect_named(out$value, c("target", "alternative"))
  expect_equal(ncol(out$value$target$coefficients), 4)
  expect_equal(rownames(out$value$target$coefficients), c("(Intercept)", "x"))
})

test_that("summary_ols uses the stored trim_prob", {
  a <- capture.output(r1 <- summary_ols(bag_vardist))
  a <- capture.output(r2 <- summary_ols(bag_vardist, trim_prob = 0.4))
  expect_equal(r1$target$coefficients[, 1],
               apply(bag_vardist$ols$target$coef, 2, agg.value, "trimmed", bag_vardist$trim_prob))
  expect_equal(r2$target$coefficients[, 1],
               apply(bag_vardist$ols$target$coef, 2, agg.value, "trimmed", 0.4))
})
