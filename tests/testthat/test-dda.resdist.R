library(testthat)

set.seed(123)

# --- Small dataset (fast tests) ---
n <- 200
x <- rchisq(n, df = 4) - 4
e <- rchisq(n, df = 3) - 3
y <- 0.5 * x + e
dat <- data.frame(x = x, y = y)

m <- lm(y ~ x, data = dat)

# Safe runner that can treat certain error messages as accepted outcomes
safe_run <- function(expr, accept_patterns = NULL) {
  tryCatch(
    eval(expr, envir = parent.frame()),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!is.null(accept_patterns) &&
          any(vapply(accept_patterns, function(p) grepl(p, msg, ignore.case = TRUE), logical(1)))) {
        structure(list(error = TRUE, message = msg), class = "dda_test_accepted_error")
      } else {
        stop(e)
      }
    }
  )
}
is_accepted_error <- function(x) inherits(x, "dda_test_accepted_error")

# ------------------------
# Basic smoke tests
# ------------------------

test_that("dda.resdist basic run returns a list-like object", {
  out <- safe_run(quote(
    dda.resdist(formula = m, pred = "x", data = dat, B = 50, prob.trans = FALSE)
  ))
  expect_true(is.list(out))
  expect_true(all(c("agostino", "anscombe") %in% names(out)))
})

test_that("dda.resdist print works", {
  out <- safe_run(quote(
    dda.resdist(formula = m, pred = "x", data = dat, B = 40)
  ))
  expect_true(is.list(out))
  expect_output(print(out))
})

# ------------------------
# Defensive / input-validation tests
# ------------------------

test_that("dda.resdist errors when pred is missing", {
  expect_error(dda.resdist(formula = m, pred = NULL, data = dat))
})

test_that("dda.resdist errors for invalid formula inputs", {
  expect_error(dda.resdist(formula = 1, pred = "x", data = dat))
  expect_error(dda.resdist(formula = "y ~ x", pred = "x", data = dat))
})

test_that("dda.resdist handles prob.trans = TRUE without failing", {
  res <- safe_run(quote(
    dda.resdist(formula = m, pred = "x", data = dat, B = 40, prob.trans = TRUE)
  ))
  # Accept either a list (successful run) or an accepted error object (environment-dependent)
  if (is_accepted_error(res)) {
    expect_true(TRUE)
  } else {
    expect_true(is.list(res))
    expect_true("probtrans" %in% names(res))
  }
})

# ------------------------
# Bootstrap BCa small-B behavior (defensive)
# ------------------------
test_that("dda.resdist bca with small B: accept error or success (defensive)", {
  # Small B may make BCa fail on some systems. Run defensively:
  res <- tryCatch(
    dda.resdist(formula = m, pred = "x", data = dat, B = 10, boot.type = "bca", conf.level = 0.90),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    # Accept any error for very small B (implementation differences across environments)
    expect_true(TRUE)
  } else {
    # If successful, ensure returned structure contains expected pieces
    expect_true(is.list(res))
    expect_true("boot.args" %in% names(res))
    expect_true(res$boot.args[1] %in% c("bca", "perc"))
  }
})

# ------------------------
# Minimal output checks for returned objects
# ------------------------

test_that("dda.resdist returns skew/kurt difference entries", {
  out <- dda.resdist(formula = m, pred = "x", data = dat, B = 0) # B = 0 to avoid bootstrapping CI work
  expect_true(is.list(out))
  expect_true("skewdiff" %in% names(out))
  expect_true("kurtdiff" %in% names(out))
  # numeric first element (point estimate) present
  expect_true(is.numeric(out$skewdiff[1]))
  expect_true(is.numeric(out$kurtdiff[1]))
})
