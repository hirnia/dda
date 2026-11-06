## ============================================================================
## testthat file: test-print_dda_bagging.R
## Tests for print.dda.bagging(), reaggregate.bagging() and ols.summary()
## ============================================================================

# --- Shared setup -----------------------------------------------------------
generate.test.data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- rchisq(n, df = 4) - 4
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e
  data.frame(x = x, y = y)
}

d <- generate.test.data()

bag.indep <- dda.bagging(
  dda.indep(y ~ x, pred = "x", data = d, B = 10,
            hetero = TRUE, nlfun = 2, diff = TRUE),
  data = d, iter = 10, progress = FALSE
)
bag.resdist <- dda.bagging(
  dda.resdist(y ~ x, pred = "x", data = d, B = 10),
  data = d, iter = 10, progress = FALSE
)
bag.vardist <- dda.bagging(
  dda.vardist(y ~ x, pred = "x", data = d, B = 10),
  data = d, iter = 10, progress = FALSE, agg.stat = "trimmed"
)

## ============================================================================
## 1. reaggregate.bagging()
## ============================================================================

test_that("reaggregate.bagging returns the object unchanged when agg.stat is NULL", {
  expect_identical(reaggregate.bagging(bag.indep, agg.stat = NULL), bag.indep)
})

test_that("reaggregate.bagging stores the new method", {
  res <- reaggregate.bagging(bag.indep, agg.stat = "median")
  expect_equal(res$agg.stat, "median")
  expect_true(is.numeric(res$aggregated.stats$hsic.yx$statistic))
})

test_that("a new trim.prob alone changes the trimmed aggregate", {
  res <- reaggregate.bagging(bag.vardist, agg.stat = "trimmed", trim.prob = 0.4)
  expect_equal(res$trim.prob, 0.4)
  expected <- agg.value(sapply(bag.vardist$bagged.results, function(b) b$agostino$outcome$statistic[1]),
                        "trimmed", trim.prob = 0.4)
  expect_equal(unname(res$aggregated.stats$agostino$outcome$statistic[1]), expected)
})

test_that("stored trim.prob and win.prob are the defaults when re-aggregating", {
  res <- reaggregate.bagging(bag.vardist, agg.stat = "trimmed")
  expect_equal(res$trim.prob, bag.vardist$trim.prob)
  expect_equal(res$aggregated.stats, bag.vardist$aggregated.stats)
})

test_that("reaggregate.bagging works for resdist objects", {
  res <- reaggregate.bagging(bag.resdist, agg.stat = "tukey")
  expect_equal(res$agg.stat, "tukey")
  expect_s3_class(res$aggregated.stats, "dda.resdist")
})

## ============================================================================
## 2. print.dda.bagging()
## ============================================================================

test_that("print shows the header, sample count and method", {
  expect_output(print(bag.indep), regexp = "BOOTSTRAP AGGREGATED DDA")
  expect_output(print(bag.indep), regexp = "Number of bootstrap samples")
  expect_output(print(bag.indep), regexp = "Aggregation method: mean")
})

test_that("print shows the tests of each object type", {
  expect_output(print(bag.indep), regexp = "HSIC")
  expect_output(print(bag.indep), regexp = "Robust BP-test")
  expect_output(print(bag.indep), regexp = "Non-linear Correlation")
  expect_output(print(bag.resdist), regexp = "Residual Distributions")
  expect_output(print(bag.vardist), regexp = "Variable Distributions")
})

test_that("print respects an agg.stat override", {
  expect_output(print(bag.indep, agg.stat = "median"), regexp = "Aggregation method: median")
})

test_that("print returns the original object invisibly", {
  out <- withVisible(print(bag.indep, agg.stat = "median"))
  expect_false(out$visible)
  expect_identical(out$value, bag.indep)
})

## ============================================================================
## 3. ols.summary()
## ============================================================================

test_that("ols.summary rejects other objects", {
  expect_error(ols.summary(list(a = 1)), regexp = "must be a dda.bagging object")
})

test_that("ols.summary prints both models with R-squared", {
  expect_output(ols.summary(bag.indep), regexp = "OLS Summary: Target Model")
  expect_output(ols.summary(bag.indep), regexp = "OLS Summary: Alternative Model")
  expect_output(ols.summary(bag.indep), regexp = "R-squared")
})

test_that("ols.summary respects an agg.stat override", {
  expect_output(ols.summary(bag.indep, agg.stat = "median"), regexp = "Aggregation method: median")
})

test_that("ols.summary returns coefficient tables invisibly", {
  out <- withVisible(ols.summary(bag.resdist))
  expect_false(out$visible)
  expect_named(out$value, c("target", "alternative"))
  expect_equal(ncol(out$value$target$coefficients), 4)
  expect_equal(rownames(out$value$target$coefficients), c("(Intercept)", "x"))
})

test_that("ols.summary uses the stored trim.prob", {
  a <- capture.output(r1 <- ols.summary(bag.vardist))
  a <- capture.output(r2 <- ols.summary(bag.vardist, trim.prob = 0.4))
  expect_equal(r1$target$coefficients[, 1],
               apply(bag.vardist$ols$target$coef, 2, agg.value, "trimmed", bag.vardist$trim.prob))
  expect_equal(r2$target$coefficients[, 1],
               apply(bag.vardist$ols$target$coef, 2, agg.value, "trimmed", 0.4))
})
