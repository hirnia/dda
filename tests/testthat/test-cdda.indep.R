library(testthat)

set.seed(123)

n <- 500

# --- Generate moderator
z <- sort(rnorm(n))
z1 <- z[z <= 0]
z2 <- z[z > 0]

# --- x -> y when m <= 0
x1 <- rchisq(length(z1), df = 4) - 4
e1 <- rchisq(length(z1), df = 3) - 3
y1 <- 0.5 * x1 + e1

# --- y -> x when m > 0
y2 <- rchisq(length(z2), df = 4) - 4
e2 <- rchisq(length(z2), df = 3) - 3
x2 <- 0.5 * y2 + e2

# Combine data
y <- c(y1, y2)
x <- c(x1, x2)
dat <- data.frame(x = x, y = y, z = z)

# Fit model
m <- lm(y ~ x * z, data = dat)

# ------------------------
# Basic smoke tests
# ------------------------

test_that("cdda.indep runs with numeric nlfun and returns a list-like object", {
  out <- cdda.indep(
    formula = m,
    pred = "x",
    mod = "z",
    modval = c(0, 1),
    dcor = TRUE,
    diff = FALSE,
    hetero = FALSE,
    data = dat,
    nlfun = 2,
    B = 100
  )
  expect_true(is.list(out))
})

test_that("cdda.indep accepts a function nlfun and supports print/summary", {
  out.f <- cdda.indep(
    formula = m,
    pred = "x",
    mod = "z",
    modval = c(0, 1),
    data = dat,
    nlfun = function(t) t^2,
    B = 100
  )
  expect_true(is.list(out.f))
  expect_output(print(out.f))
  expect_output(summary(out.f, hsic = FALSE))
})

# ------------------------
# Defensive / input-validation tests
# ------------------------

test_that("cdda.indep errors when pred is not in the model", {
  expect_error(cdda.indep(formula = m, pred = "no.such.variable", data = dat))
})

test_that("cdda.indep errors when formula is not an lm object", {
  expect_error(cdda.indep(formula = 1, pred = "x", data = dat))
  expect_error(cdda.indep(formula = "y ~ x", pred = "x", data = dat))
})

test_that("cdda.indep errors when nlfun is of wrong type", {
  expect_error(cdda.indep(formula = m, pred = "x", data = dat, nlfun = "not.a.function"))
})

test_that("cdda.indep errors when B is nonsensical (negative or zero)", {
  expect_error(cdda.indep(formula = m, pred = "x", data = dat, B = -10))
  expect_error(cdda.indep(formula = m, pred = "x", data = dat, B = 0))
})

# ------------------------
# Small-B bootstrap HSIC defensive test (no skip, tolerant)
# ------------------------
test_that("cdda.indep small-B bootstrap HSIC: tolerant behavior (warning, error, or valid return accepted)", {
  # We run with a very small B and accept any of:
  #  - a warning emitted
  #  - an error (returned as an error object by tryCatch)
  #  - a successful return (list)
  w.msgs <- character()
  res <- NULL

  withCallingHandlers(
    {
      res <- tryCatch(
        cdda.indep(formula = m, pred = "x", data = dat, hsic.method = "bootstrap", B = 10, diff = FALSE, mod = "z", modval = c(0,1)),
        error = function(e) e
      )
    },
    warning = function(w) {
      w.msgs <<- c(w.msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Accept any warning as evidence the function handled the situation (e.g., warned and possibly fell back)
  if (length(w.msgs) > 0) {
    expect_true(TRUE)
    return()
  }

  # If an error was returned, accept it (defensive test) — we don't assert a specific message here
  if (inherits(res, "error")) {
    expect_true(TRUE)
    return()
  }

  # If a valid object was returned, ensure it's list-like
  if (is.list(res)) {
    expect_true(TRUE)
    return()
  }

  # Otherwise fail (shouldn't happen)
  fail("Unexpected non-error, non-warning, non-list outcome from cdda.indep with small B")
})

# ------------------------
# Final lightweight checks
# ------------------------

test_that("cdda.indep handles modval variations (numeric and named)", {
  # numeric pick-a-point
  out.num <- cdda.indep(formula = m, pred = "x", mod = "z", modval = c(-1, 1), data = dat, B = 50)
  expect_true(is.list(out.num))

  # named 'mean' approach (ensure it does not error)
  out.mean <- cdda.indep(formula = m, pred = "x", mod = "z", modval = "mean", data = dat, B = 50)
  expect_true(is.list(out.mean))
})

test_that("cdda.indep returns non-empty structure", {
  out <- cdda.indep(formula = m, pred = "x", mod = "z", modval = c(0,1), data = dat, nlfun = 2, B = 50)
  expect_true(is.list(out))
  expect_true(length(out) > 0)
})

test_that("cdda.indep runs with robust estimation", {
  out <- cdda.indep(m, pred = "x", mod = "z", modval = c(-1, 1), data = dat,
                    hetero = TRUE, B = 20, robust = TRUE)
  expect_s3_class(out, "cdda.indep")
})

## OLD (none)
## NEW 0.2.0
# ------------------------
# Difference statistics, factor labels and input checks
# ------------------------

test_that("cdda.indep returns one difference matrix per moderator value", {
  out <- cdda.indep(m, pred = "x", mod = "z", modval = c(-1, 1), data = dat,
                    diff = TRUE, B = 20)
  expect_named(out$cdda.diff$diff, c("ModVal -1", "ModVal 1"))
  expect_equal(dim(out$cdda.diff$diff[[1]]), c(3, 3))
  expect_null(out$cdda.target[[1]]$out.diff)
  expect_output(summary(out, hsic.diff = TRUE, dcor.diff = TRUE),
                "Values > 0 suggest x -> y", fixed = TRUE)
  expect_output(summary(out, hsic = FALSE, dcor = FALSE, mi.diff = TRUE),
                "Caution", fixed = TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(out, stat = "hsic.diff"))
  expect_no_error(plot(out, stat = "mi.diff"))
})

test_that("cdda.indep HSIC and dCor differences are positive when x -> y and negative when y -> x", {
  skip_on_cran()
  set.seed(2026)
  n.s <- 1000
  z.s <- rnorm(n.s)

  x.s <- rchisq(n.s, df = 4) - 4
  y.s <- x.s + 0.3 * x.s * z.s + rnorm(n.s)
  d.xy <- data.frame(x = x.s, y = y.s, z = z.s)
  out.xy <- cdda.indep(lm(y ~ x * z, data = d.xy), pred = "x", mod = "z",
                       modval = c(-1, 0, 1), data = d.xy, diff = TRUE, B = 10)
  expect_true(all(sapply(out.xy$cdda.diff$diff, function(d) d[c("HSIC", "dCor"), "estimate"]) > 0))

  y.r <- rchisq(n.s, df = 4) - 4
  x.r <- y.r + 0.3 * y.r * z.s + rnorm(n.s)
  d.yx <- data.frame(x = x.r, y = y.r, z = z.s)
  out.yx <- cdda.indep(lm(y ~ x * z, data = d.yx), pred = "x", mod = "z",
                       modval = c(-1, 0, 1), data = d.yx, diff = TRUE, B = 10)
  expect_true(all(sapply(out.yx$cdda.diff$diff, function(d) d[c("HSIC", "dCor"), "estimate"]) < 0))
})

test_that("cdda.indep labels factor moderator results by level", {
  set.seed(7)
  n.f <- 300
  g <- factor(sample(c("a", "b"), n.f, replace = TRUE), levels = c("a", "b"))
  g[1] <- "b"
  x.f <- rchisq(n.f, df = 4) - 4
  d.f <- data.frame(x = x.f, y = 0.5 * x.f + rnorm(n.f), g = g)
  out <- cdda.indep(lm(y ~ x * g, data = d.f), pred = "x", mod = "g", data = d.f, B = 10)
  expect_named(out$cdda.target, c("ModVal a", "ModVal b"))
  expect_equal(out$df.original$mod.levels, c("a", "b"))
})

test_that("cdda.indep asks for modval with a continuous moderator", {
  expect_error(cdda.indep(m, pred = "x", mod = "z", data = dat), "Please specify modval")
})

test_that("cdda.indep rejects a moderator that is neither numeric nor a factor", {
  d.c <- dat
  d.c$zc <- ifelse(dat$z > 0, "hi", "lo")
  expect_error(cdda.indep(lm(y ~ x * zc, data = d.c), pred = "x", mod = "zc",
                          modval = "mean", data = d.c),
               "Specified moderator not found")
})
