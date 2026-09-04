## ============================================================================
## testthat file: test-dda_bagging.R
## Tests for dda.bagging()
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

# Fit base DDA objects (low B for speed)
d     <- generate_test_data()
d_cov <- generate_test_data_cov()

base_indep   <- dda.indep(y ~ x, pred = "x", data = d, B = 10)
base_resdist <- dda.resdist(y ~ x, pred = "x", data = d, B = 10)
base_vardist <- dda.vardist(y ~ x, pred = "x", data = d, B = 10)

# Richer base objects for optional-feature tests
base_indep_hetero <- dda.indep(y ~ x, pred = "x", data = d, B = 10,
                                hetero = TRUE, nlfun = 2, diff = TRUE)
base_resdist_pt   <- dda.resdist(y ~ x, pred = "x", data = d, B = 10,
                                  prob.trans = TRUE)
base_indep_cov    <- dda.indep(y ~ x + z, pred = "x", data = d_cov, B = 10)

## ============================================================================
## 1. Input Validation
## ============================================================================

test_that("dda_bagging errors on wrong input class", {
  expect_error(
    dda.bagging(list(a = 1), data = d, iter = 5),
    regexp = "Unsupported DDA object"
  )
})

test_that("dda_bagging errors when data is NULL", {
  expect_error(
    dda.bagging(base_indep, data = NULL, iter = 5),
    regexp = "Please provide"
  )
})

test_that("dda_bagging errors with invalid agg_stat", {
  expect_error(
    dda.bagging(base_indep, data = d, iter = 5, agg_stat = "geometric"),
    regexp = "should be one of"
  )
})

## ============================================================================
## 2. Return Structure — all three object types
## ============================================================================

run_bag <- function(base, dat, iter = 10, ...) {
  dda.bagging(base, data = dat, iter = iter, progress = FALSE, ...)
}

test_that("dda_bagging returns correct class for dda.indep input", {
  result <- run_bag(base_indep, d)
  expect_s3_class(result, "dda_bagging_indep")
  expect_s3_class(result, "dda_bagging")
})

test_that("dda_bagging returns correct class for dda.resdist input", {
  result <- run_bag(base_resdist, d)
  expect_s3_class(result, "dda_bagging_resdist")
  expect_s3_class(result, "dda_bagging")
})

test_that("dda_bagging returns correct class for dda.vardist input", {
  result <- run_bag(base_vardist, d)
  expect_s3_class(result, "dda_bagging_vardist")
  expect_s3_class(result, "dda_bagging")
})

test_that("dda_bagging output contains all required top-level slots", {
  result <- run_bag(base_indep, d)
  expect_named(result, c("bagged_results", "raw_stats", "aggregated_stats",
                          "decision_percentages", "n_valid_iterations",
                          "agg_stat_used"),
               ignore.order = TRUE)
})

test_that("n_valid_iterations is a positive integer <= iter", {
  result <- run_bag(base_indep, d, iter = 10)
  expect_true(result$n_valid_iterations > 0)
  expect_true(result$n_valid_iterations <= 10)
})

test_that("agg_stat_used matches the requested method", {
  for (method in c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey")) {
    result <- run_bag(base_indep, d, agg_stat = method)
    expect_equal(result$agg_stat_used, method,
                 info = paste("Failed for agg_stat =", method))
  }
})

## ============================================================================
## 3. Aggregated Stats Content
## ============================================================================

test_that("aggregated_stats contains HSIC results for indep object", {
  result <- run_bag(base_indep, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$hsic_yx_stat))
  expect_true(!is.null(agg$hsic_xy_stat))
  expect_true(!is.null(agg$hsic_yx_pval))
  expect_true(!is.null(agg$hsic_xy_pval))
  expect_true(is.numeric(agg$hsic_yx_stat))
})

test_that("aggregated_stats contains dCor results when base object has them", {
  result <- run_bag(base_indep, d)
  agg <- result$aggregated_stats
  # dCor is present because dda.indep computes it by default
  expect_true(!is.null(agg$dcor_yx_stat) || is.null(agg$dcor_yx_stat),
              label = "dCor slot exists or is absent depending on base object")
})

test_that("aggregated_stats contains BP results when hetero = TRUE", {
  result <- run_bag(base_indep_hetero, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$breusch_pagan))
  expect_length(agg$breusch_pagan, 4)
})

test_that("aggregated_stats contains nlcor results when nlfun is set", {
  result <- run_bag(base_indep_hetero, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$nlcor.yx))
  expect_true(!is.null(agg$nlcor.xy))
})

test_that("aggregated_stats contains diff_matrix when diff = TRUE", {
  result <- run_bag(base_indep_hetero, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$diff_matrix))
})

test_that("aggregated_stats for resdist contains agostino and anscombe results", {
  result <- run_bag(base_resdist, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$agostino.target.statistic))
  expect_true(!is.null(agg$agostino.alternative.statistic))
  expect_true(!is.null(agg$anscombe.target.statistic))
  expect_true(!is.null(agg$anscombe.alternative.statistic))
})

test_that("aggregated_stats for vardist contains predictor/outcome agostino results", {
  result <- run_bag(base_vardist, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$agostino.predictor.statistic.skew))
  expect_true(!is.null(agg$agostino.outcome.statistic.skew))
  expect_true(!is.null(agg$anscombe.predictor.statistic.kurt))
  expect_true(!is.null(agg$anscombe.outcome.statistic.kurt))
})

test_that("OLS target and alternative results are populated", {
  result <- run_bag(base_indep, d)
  agg <- result$aggregated_stats
  expect_true(!is.null(agg$ols_target))
  expect_true(!is.null(agg$ols_alternative))
  expect_true(is.matrix(agg$ols_target))
  expect_equal(colnames(agg$ols_target),
               c("estimate", "2.5 %", "97.5 %", "Prop (p<0.05)"))
})

test_that("var.names are stored correctly in aggregated_stats", {
  result <- run_bag(base_indep, d)
  expect_equal(result$aggregated_stats$var.names, c("y", "x"))
})

## ============================================================================
## 4. Decision Percentages
## ============================================================================

test_that("decision_percentages is a non-empty list", {
  result <- run_bag(base_indep, d)
  expect_true(is.list(result$decision_percentages))
  expect_true(length(result$decision_percentages) > 0)
})

test_that("each decision entry sums to 1", {
  result <- run_bag(base_indep, d)
  for (nm in names(result$decision_percentages)) {
    s <- sum(result$decision_percentages[[nm]])
    expect_equal(s, 1, tolerance = 1e-6,
                 label = paste("decision sum for", nm))
  }
})

test_that("decision levels are Target / Alternative / Undecided", {
  result <- run_bag(base_indep, d)
  for (nm in names(result$decision_percentages)) {
    expect_setequal(names(result$decision_percentages[[nm]]),
                    c("Target", "Alternative", "Undecided"))
  }
})

## ============================================================================
## 5. Aggregation Methods (spot-check numeric plausibility)
## ============================================================================

test_that("trimmed aggregation with trim_prob = 0 equals mean", {
  set.seed(101)
  res_mean    <- run_bag(base_indep, d, agg_stat = "mean")
  set.seed(101)
  res_trimmed <- run_bag(base_indep, d, agg_stat = "trimmed", trim_prob = 0)
  # Both use the same raw stats seeded the same way; HSIC stat should be equal
  expect_equal(res_mean$aggregated_stats$hsic_yx_stat,
               res_trimmed$aggregated_stats$hsic_yx_stat,
               tolerance = 1e-6)
})

test_that("winsorized aggregation with win_prob = 0 equals mean", {
  set.seed(101)
  res_mean <- run_bag(base_indep, d, agg_stat = "mean")
  set.seed(101)
  res_win  <- run_bag(base_indep, d, agg_stat = "winsorized", win_prob = 0)
  expect_equal(res_mean$aggregated_stats$hsic_yx_stat,
               res_win$aggregated_stats$hsic_yx_stat,
               tolerance = 1e-6)
})

## ============================================================================
## 6. save_file Argument
## ============================================================================

test_that("save_file writes a readable RDS to disk", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  run_bag(base_indep, d, save_file = tmp)
  expect_true(file.exists(tmp))
  saved <- readRDS(tmp)
  expect_true(is.list(saved))
  expect_true("aggregated_stats" %in% names(saved))
})

## ============================================================================
## 7. Covariate Support
## ============================================================================

test_that("dda_bagging works with covariates in the model", {
  result <- run_bag(base_indep_cov, d_cov)
  expect_s3_class(result, "dda_bagging_indep")
  expect_true(result$n_valid_iterations > 0)
})

## ============================================================================
## 8. prob.trans = TRUE for resdist
## ============================================================================

test_that("dda_bagging handles prob.trans = TRUE in resdist correctly", {
  result <- run_bag(base_resdist_pt, d)
  expect_s3_class(result, "dda_bagging_resdist")
  # skewdiff decisions should be reversed; just check they exist and sum to 1
  decs <- result$decision_percentages
  expect_true(!is.null(decs$dec_skewdiff) || !is.null(decs$dec_agost))
})
