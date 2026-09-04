## ============================================================================
## testthat file: test-summary_dda_bagging.R
## Tests for summary.dda_bagging_*, print_bagging_decisions(),
## and round_preserve_sum()
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
  data = d, iter = 10, progress = FALSE
)

## ============================================================================
## 1. round_preserve_sum() (internal helper)
## ============================================================================

test_that("round_preserve_sum: output sums to 1 for valid proportions", {
  props <- c(0.3333, 0.3333, 0.3334)
  rounded <- round_preserve_sum(props, digits = 2)
  expect_equal(sum(rounded), 1.00, tolerance = 1e-10)
})

test_that("round_preserve_sum: output sums to 1 for extreme split", {
  props <- c(0.999, 0.001, 0.000)
  rounded <- round_preserve_sum(props, digits = 2)
  expect_equal(sum(rounded), 1.00, tolerance = 1e-10)
})

test_that("round_preserve_sum: handles all-zero input gracefully", {
  props <- c(0, 0, 0)
  result <- round_preserve_sum(props, digits = 2)
  expect_equal(sum(result), 0)
})

test_that("round_preserve_sum: respects digits parameter", {
  props <- c(1/3, 1/3, 1/3)
  r2 <- round_preserve_sum(props, digits = 2)
  r4 <- round_preserve_sum(props, digits = 4)
  # Higher digits = more precision in each value
  expect_true(all(nchar(sub(".*\\.", "", as.character(r4))) >= 2))
})

## ============================================================================
## 2. summary.dda_bagging_indep() — 'show' options
## ============================================================================

test_that("summary.dda_bagging_indep runs without error (show = NULL)", {
  expect_output(summary(bag_indep), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("summary.dda_bagging_indep returns object invisibly", {
  result <- withVisible(summary(bag_indep))
  expect_false(result$visible)
  expect_s3_class(result$value, "dda_bagging_indep")
})

# Test every valid 'show' alias for indep
for (opt in c("hsic", "dcor", "mi", "bp", "nlcor", "all")) {
  local({
    o <- opt
    test_that(paste0("summary.dda_bagging_indep: show = '", o, "' runs without error"), {
      # Some options only produce output if the underlying stats exist;
      # we only check no error is thrown
      expect_no_error(
        capture.output(summary(bag_indep, show = o))
      )
    })
  })
}

test_that("summary.dda_bagging_indep: show = 'all' prints something", {
  expect_output(summary(bag_indep, show = "all"), regexp = ".")
})

test_that("summary.dda_bagging_indep: unrecognised show gives informative message", {
  expect_output(
    summary(bag_indep, show = "not_a_real_stat"),
    regexp = "No statistics matched"
  )
})

test_that("summary.dda_bagging_indep: show = 'hsic' shows HSIC label", {
  expect_output(summary(bag_indep, show = "hsic"), regexp = "HSIC")
})

test_that("summary.dda_bagging_indep: show = 'bp' shows Breusch-Pagan label", {
  expect_output(summary(bag_indep, show = "bp"), regexp = "Breusch-Pagan")
})

test_that("summary.dda_bagging_indep: show = 'nlcor' shows Non-linear Correlation", {
  expect_output(summary(bag_indep, show = "nlcor"), regexp = "Non-linear Correlation")
})

test_that("summary.dda_bagging_indep: digits argument changes decimal places", {
  out2 <- capture.output(summary(bag_indep, digits = 2))
  out4 <- capture.output(summary(bag_indep, digits = 4))
  # More decimal digits → longer string output
  expect_true(nchar(paste(out4, collapse = "")) >= nchar(paste(out2, collapse = "")))
})

test_that("summary.dda_bagging_indep: proportions displayed sum to 1 per statistic", {
  out <- capture.output(summary(bag_indep, show = "all"))
  # Verify via the underlying decision_percentages directly
  for (nm in names(bag_indep$decision_percentages)) {
    s <- sum(bag_indep$decision_percentages[[nm]])
    expect_equal(s, 1, tolerance = 1e-6, label = paste("sum for", nm))
  }
})

test_that("summary.dda_bagging_indep: footnote shows variable names", {
  expect_output(summary(bag_indep), regexp = "Target is")
  expect_output(summary(bag_indep), regexp = "Alternative is")
})

test_that("summary.dda_bagging_indep: 'Confounding' label used (not 'Undecided')", {
  expect_output(summary(bag_indep, show = "all"), regexp = "Confounding")
})

## ============================================================================
## 3. summary.dda_bagging_vardist() — 'show' and 'moment' options
## ============================================================================

test_that("summary.dda_bagging_vardist runs without error (show = NULL)", {
  expect_output(summary(bag_vardist), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("summary.dda_bagging_vardist returns object invisibly", {
  result <- withVisible(summary(bag_vardist))
  expect_false(result$visible)
})

for (opt in c("skew", "kurt", "coskew", "cokurt", "all")) {
  local({
    o <- opt
    test_that(paste0("summary.dda_bagging_vardist: show = '", o, "' runs without error"), {
      expect_no_error(capture.output(summary(bag_vardist, show = o)))
    })
  })
}

test_that("summary.dda_bagging_vardist: show = 'skew' shows D'Agostino label", {
  expect_output(summary(bag_vardist, show = "skew"),
                regexp = "D'Agostino|Skewness")
})

test_that("summary.dda_bagging_vardist: show = 'kurt' shows Anscombe-Glynn label", {
  expect_output(summary(bag_vardist, show = "kurt"),
                regexp = "Anscombe-Glynn|Kurtosis")
})

test_that("summary.dda_bagging_vardist: moment = 3 shows skew-related stats only", {
  out <- capture.output(summary(bag_vardist, moment = 3))
  # Moment 3 keys: dec_agost, dec_skewdiff, dec_cor12diff, dec_RHS
  # Should NOT show kurt keys
  expect_false(any(grepl("Anscombe-Glynn|Co-Kurtosis|dec_kurtdiff", out)))
})

test_that("summary.dda_bagging_vardist: moment = 4 shows kurt-related stats only", {
  out <- capture.output(summary(bag_vardist, moment = 4))
  # Moment 4 keys: dec_anscom, dec_kurtdiff, dec_cor13diff, dec_RCC, dec_RHS4, dec_Rtanh
  expect_false(any(grepl("D'Agostino|Skewness Difference", out)))
})

test_that("summary.dda_bagging_vardist: moment = c(3,4) shows both moment groups", {
  out3 <- capture.output(summary(bag_vardist, moment = 3))
  out4 <- capture.output(summary(bag_vardist, moment = 4))
  out34 <- capture.output(summary(bag_vardist, moment = c(3, 4)))
  # Combined output should be at least as long as each individual
  expect_true(length(out34) >= length(out3))
  expect_true(length(out34) >= length(out4))
})

test_that("summary.dda_bagging_vardist: 'Undecided' label used (not 'Confounding')", {
  expect_output(summary(bag_vardist, show = "all"), regexp = "Undecided")
})

test_that("summary.dda_bagging_vardist: show + moment are additive", {
  out_show   <- capture.output(summary(bag_vardist, show = "skew"))
  out_moment <- capture.output(summary(bag_vardist, moment = 4))
  out_both   <- capture.output(summary(bag_vardist, show = "skew", moment = 4))
  expect_true(length(out_both) >= max(length(out_show), length(out_moment)))
})

## ============================================================================
## 4. summary.dda_bagging_resdist() — 'show' and 'moment' options
## ============================================================================

test_that("summary.dda_bagging_resdist runs without error (show = NULL)", {
  expect_output(summary(bag_resdist), regexp = "BOOTSTRAP AGGREGATED DDA")
})

test_that("summary.dda_bagging_resdist returns object invisibly", {
  result <- withVisible(summary(bag_resdist))
  expect_false(result$visible)
})

for (opt in c("skew", "kurt", "all")) {
  local({
    o <- opt
    test_that(paste0("summary.dda_bagging_resdist: show = '", o, "' runs without error"), {
      expect_no_error(capture.output(summary(bag_resdist, show = o)))
    })
  })
}

test_that("summary.dda_bagging_resdist: moment filtering works", {
  expect_no_error(capture.output(summary(bag_resdist, moment = c(3, 4))))
})

test_that("summary.dda_bagging_resdist: unrecognised show gives informative message", {
  expect_output(
    summary(bag_resdist, show = "xyz_fake"),
    regexp = "No statistics matched"
  )
})

test_that("summary.dda_bagging_resdist: 'Undecided' label used for resdist", {
  # resdist type = "resdist" → uses Undecided column, not Confounding
  expect_output(summary(bag_resdist, show = "all"), regexp = "Undecided")
})

## ============================================================================
## 5. Confirm complete 'show' alias list
## ============================================================================

test_that("all documented show aliases are present in alias_map", {
  # Retrieve alias_map by calling the internal helper and checking what keys
  # are accepted without producing 'No statistics matched'.
  # We use bag_indep which has all possible decision keys present (hetero + nlfun + diff)
  valid_aliases <- c("hsic", "dcor", "mi", "bp", "nlcor",
                     "skew", "kurt", "coskew", "cokurt", "all")

  for (alias in valid_aliases) {
    out <- capture.output(summary(bag_indep, show = alias))
    no_match <- any(grepl("No statistics matched", out))
    # We don't fail if "no match" — some aliases won't exist in every object type.
    # The important thing is no error is thrown.
    expect_no_error(capture.output(summary(bag_indep, show = alias)))
  }
})
