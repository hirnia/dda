## ============================================================================
## testthat file: test-summary_dda_bagging.R
## Tests for summary.dda.bagging() and round_preserve_sum()
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
  data = d, iter = 10, progress = FALSE
)

## ============================================================================
## 1. round_preserve_sum()
## ============================================================================

test_that("round_preserve_sum keeps the sum at 1", {
  expect_equal(sum(round_preserve_sum(c(1/3, 1/3, 1/3), 2)), 1)
  expect_equal(sum(round_preserve_sum(c(0.999, 0.001, 0), 2)), 1)
  expect_equal(sum(round_preserve_sum(c(0.125, 0.375, 0.5), 1)), 1)
})

test_that("round_preserve_sum handles zeros and missing values", {
  expect_equal(round_preserve_sum(c(0, 0, 0), 2), c(0, 0, 0))
  expect_equal(round_preserve_sum(c(NaN, NaN), 2), c(NaN, NaN))
})

## ============================================================================
## 2. summary.dda.bagging()
## ============================================================================

test_that("summary prints the header and footnote for each object type", {
  expect_output(summary(bag.indep), regexp = "Independence Properties")
  expect_output(summary(bag.resdist), regexp = "Residual Distributions")
  expect_output(summary(bag.vardist), regexp = "Variable Distributions")
  expect_output(summary(bag.indep), regexp = "Target is x -> y")
  expect_output(summary(bag.indep), regexp = "Alternative is y -> x")
})

test_that("summary returns the decision proportions invisibly", {
  out <- withVisible(summary(bag.indep))
  expect_false(out$visible)
  expect_identical(out$value, bag.indep$decision.proportions)
})

test_that("Confounding appears only for independence properties", {
  expect_output(summary(bag.indep), regexp = "Confounding")
  expect_false(any(grepl("Confounding", capture.output(summary(bag.vardist)))))
})

test_that("show filters the reported statistics", {
  out <- withVisible(summary(bag.indep, show = "hsic"))$value
  expect_setequal(rownames(out), c("hsic", "hsic.diff"))
  expect_output(summary(bag.indep, show = "bp"), regexp = "Robust BP-test")
  expect_output(summary(bag.indep, show = "nlcor"), regexp = "Non-linear Correlation")

  out <- withVisible(summary(bag.vardist, show = c("skew", "coskew")))$value
  expect_true(all(rownames(out) %in% c("agostino", "skewdiff", "cor12diff", "RHS")))
  expect_false(any(grepl("Anscombe-Glynn", capture.output(summary(bag.vardist, show = "skew")))))
})

test_that("every show option runs on a matching object", {
  for (s in c("hsic", "dcor", "mi", "bp", "nlcor")) {
    expect_no_error(capture.output(summary(bag.indep, show = s)))
  }
  for (s in c("skew", "kurt", "coskew", "cokurt")) {
    expect_no_error(capture.output(summary(bag.vardist, show = s)))
    expect_no_error(capture.output(summary(bag.resdist, show = s)))
  }
})

test_that("unknown or unavailable show options give an error", {
  expect_error(summary(bag.indep, show = "not.a.stat"), regexp = "Unknown statistic")
  expect_error(summary(bag.vardist, show = "hsic"), regexp = "None of the requested")
})

test_that("rows follow the order of the base print methods", {
  expect_equal(rownames(bag.vardist$decision.proportions),
               c("agostino", "anscombe", "skewdiff", "kurtdiff", "cor12diff", "RHS", "cor13diff", "RCC", "Rtanh"))
  expect_equal(rownames(bag.resdist$decision.proportions),
               c("agostino", "anscombe", "skewdiff.test", "kurtdiff.test", "skewdiff", "kurtdiff",
                 "cor12diff", "RHS3", "cor13diff", "RHS4", "RCC"))
  expect_equal(rownames(bag.indep$decision.proportions),
               c("hsic", "dcor", "bp", "nlcor", "hsic.diff", "dcor.diff", "mi.diff"))
})

test_that("digits sets the decimal places", {
  expect_output(summary(bag.indep, digits = 3), regexp = "[0-9]\\.[0-9]{3}")
})
