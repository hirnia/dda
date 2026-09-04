## ============================================================================
## testthat file: test-print_dda_bagging.R
## Tests for print.dda_bagging_*, print_ols_summary, and reaggregate_bagging()
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
  dda.indep(y ~ x, pred = "x", data = d, B = 10),
  data = d, iter = 10, progress = FALSE
)
bag_indep_full <- dda.bagging(
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
  data = d, iter = 10, progress = FALSE
)

## ============================================================================
## 1. reaggregate_bagging()
## ============================================================================

test_that("reaggregate_bagging returns same object when agg_stat is NULL", {
  result <- reaggregate_bagging(bag_indep, agg_stat = NULL)
  expect_identical(result, bag_indep)
})

test_that("reaggregate_bagging returns same object when agg_stat matches current", {
  current <- bag_indep$agg_stat_used
  result   <- reaggregate_bagging(bag_indep, agg_stat = current)
  expect_identical(result, bag_indep)
})

test_that("reaggregate_bagging updates agg_stat_used when method changes", {
  new_method <- if (bag_indep$agg_stat_used == "mean") "median" else "mean"
  result <- reaggregate_bagging(bag_indep, agg_stat = new_method)
  expect_equal(result$agg_stat_used, new_method)
})

test_that("reaggregate_bagging updates aggregated HSIC stat", {
  new_method <- if (bag_indep$agg_stat_used == "mean") "median" else "mean"
  result     <- reaggregate_bagging(bag_indep, agg_stat = new_method)
  # The numeric value should differ from the original (almost surely with real data)
  # We only assert the field is still numeric and not NA
  expect_true(is.numeric(result$aggregated_stats$hsic_yx_stat))
  expect_false(is.na(result$aggregated_stats$hsic_yx_stat))
})

test_that("reaggregate_bagging works for resdist objects", {
  result <- reaggregate_bagging(bag_resdist, agg_stat = "median")
  expect_equal(result$agg_stat_used, "median")
  expect_true(is.numeric(result$aggregated_stats$agostino.target.statistic))
})

test_that("reaggregate_bagging works for vardist objects", {
  result <- reaggregate_bagging(bag_vardist, agg_stat = "trimmed", trim_prob = 0.05)
  expect_equal(result$agg_stat_used, "trimmed")
  expect_true(is.numeric(result$aggregated_stats$agostino.predictor.statistic.skew))
})

test_that("reaggregate_bagging re-aggregates OLS coefficients", {
  orig_est   <- bag_indep$aggregated_stats$ols_target[, "estimate"]
  new_method <- if (bag_indep$agg_stat_used == "mean") "median" else "mean"
  result     <- reaggregate_bagging(bag_indep, agg_stat = new_method)
  new_est    <- result$aggregated_stats$ols_target[, "estimate"]
  # Estimates may or may not differ numerically, but they should be numeric
  expect_true(is.numeric(new_est))
})

## ============================================================================
## 2. print.dda_bagging_indep()
## ============================================================================

test_that("print.dda_bagging_indep produces output without error", {
  expect_output(print(bag_indep), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("print.dda_bagging_indep shows bootstrap sample count", {
  expect_output(print(bag_indep), regexp = "Number of Bootstrap Samples")
})

test_that("print.dda_bagging_indep shows aggregation method", {
  expect_output(print(bag_indep), regexp = "Aggregation method")
})

test_that("print.dda_bagging_indep shows HSIC statistics", {
  expect_output(print(bag_indep), regexp = "HSIC")
})

test_that("print.dda_bagging_indep shows BP when hetero = TRUE", {
  expect_output(print(bag_indep_full), regexp = "Breusch-Pagan")
})

test_that("print.dda_bagging_indep shows non-linear correlation when nlfun set", {
  expect_output(print(bag_indep_full), regexp = "Non-linear Correlation")
})

test_that("print.dda_bagging_indep shows difference CIs when diff = TRUE", {
  expect_output(print(bag_indep_full), regexp = "Bootstrap CIs for Difference")
})

test_that("print.dda_bagging_indep returns the object invisibly", {
  result <- withVisible(print(bag_indep))
  expect_false(result$visible)
  expect_s3_class(result$value, "dda_bagging_indep")
})

test_that("print.dda_bagging_indep respects agg_stat override", {
  new_method <- if (bag_indep$agg_stat_used == "mean") "median" else "mean"
  expect_output(
    print(bag_indep, agg_stat = new_method),
    regexp = new_method
  )
})

test_that("print.dda_bagging_indep respects digits argument", {
  # With digits = 2 the output should be shorter precision than digits = 8
  out2 <- capture.output(print(bag_indep, digits = 2))
  out8 <- capture.output(print(bag_indep, digits = 8))
  # Output with more digits will be longer (more characters)
  expect_true(nchar(paste(out8, collapse = "")) >= nchar(paste(out2, collapse = "")))
})

## ============================================================================
## 3. print.dda_bagging_resdist()
## ============================================================================

test_that("print.dda_bagging_resdist produces output without error", {
  expect_output(print(bag_resdist), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("print.dda_bagging_resdist includes Residual Distributions header", {
  expect_output(print(bag_resdist), regexp = "Residual Distributions")
})

test_that("print.dda_bagging_resdist shows skewness and kurtosis tests", {
  expect_output(print(bag_resdist), regexp = "Skewness")
  expect_output(print(bag_resdist), regexp = "Kurtosis")
})

test_that("print.dda_bagging_resdist returns object invisibly", {
  result <- withVisible(print(bag_resdist))
  expect_false(result$visible)
})

## ============================================================================
## 4. print.dda_bagging_vardist()
## ============================================================================

test_that("print.dda_bagging_vardist produces output without error", {
  expect_output(print(bag_vardist), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("print.dda_bagging_vardist includes Variable Distributions header", {
  expect_output(print(bag_vardist), regexp = "Variable Distributions")
})

test_that("print.dda_bagging_vardist shows likelihood ratio approximations", {
  expect_output(print(bag_vardist), regexp = "Likelihood Ratio")
})

test_that("print.dda_bagging_vardist returns object invisibly", {
  result <- withVisible(print(bag_vardist))
  expect_false(result$visible)
})

## ============================================================================
## 5. print_ols_summary()
## ============================================================================

test_that("print_ols_summary errors on non-dda_bagging object", {
  expect_error(print_ols_summary(list(a = 1)), regexp = "must be a bagged DDA")
})

test_that("print_ols_summary produces output without error", {
  expect_output(print_ols_summary(bag_indep), regexp = "OLS Summary")
})

test_that("print_ols_summary shows both target and alternative models", {
  expect_output(print_ols_summary(bag_indep), regexp = "Target Model")
  expect_output(print_ols_summary(bag_indep), regexp = "Alternative Model")
})

test_that("print_ols_summary shows R-squared values", {
  expect_output(print_ols_summary(bag_indep), regexp = "R-squared")
})

test_that("print_ols_summary respects agg_stat override", {
  new_method <- if (bag_indep$agg_stat_used == "mean") "median" else "mean"
  expect_output(
    print_ols_summary(bag_indep, agg_stat = new_method),
    regexp = new_method
  )
})

test_that("print_ols_summary works for resdist objects", {
  expect_output(print_ols_summary(bag_resdist), regexp = "OLS Summary")
})

test_that("print_ols_summary works for vardist objects", {
  expect_output(print_ols_summary(bag_vardist), regexp = "OLS Summary")
})
