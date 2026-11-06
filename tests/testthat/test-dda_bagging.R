## ============================================================================
## testthat file: test-dda_bagging.R
## Tests for dda.bagging(), dda.decisions() and agg.value()
## ============================================================================

# --- Shared setup -----------------------------------------------------------
generate.test.data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- rchisq(n, df = 4) - 4
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e
  data.frame(x = x, y = y)
}

generate.test.data.cov <- function(n = 200, seed = 42) {
  set.seed(seed)
  z <- rnorm(n)
  x <- rchisq(n, df = 4) - 4 + 0.3 * z
  e <- rchisq(n, df = 3) - 3
  y <- 0.5 * x + e + 0.2 * z
  data.frame(x = x, y = y, z = z)
}

d     <- generate.test.data()
d.cov <- generate.test.data.cov()

base.indep      <- dda.indep(y ~ x, pred = "x", data = d, B = 10)
base.indep.full <- dda.indep(y ~ x, pred = "x", data = d, B = 10,
                             hetero = TRUE, nlfun = 2, diff = TRUE)
base.indep.cov  <- dda.indep(y ~ x + z, pred = "x", data = d.cov, B = 10)
base.resdist    <- dda.resdist(y ~ x, pred = "x", data = d, B = 10)
base.resdist.pt <- dda.resdist(y ~ x, pred = "x", data = d, B = 10, prob.trans = TRUE)
base.vardist    <- dda.vardist(y ~ x, pred = "x", data = d, B = 10)

run.bag <- function(base, dat, iter = 10, ...) {
  dda.bagging(base, data = dat, iter = iter, progress = FALSE, ...)
}

bag.indep   <- run.bag(base.indep.full, d)
bag.resdist <- run.bag(base.resdist, d)
bag.vardist <- run.bag(base.vardist, d)

## ============================================================================
## 1. Input validation
## ============================================================================

test_that("dda.bagging rejects objects that are not DDA results", {
  expect_error(dda.bagging(list(a = 1), data = d, iter = 5),
               regexp = "must be a dda.indep, dda.resdist, or dda.vardist")
})

test_that("dda.bagging requires a data.frame", {
  expect_error(dda.bagging(base.indep, iter = 5), regexp = "Please provide")
  expect_error(dda.bagging(base.indep, data = NULL, iter = 5), regexp = "Please provide")
})

test_that("dda.bagging rejects an unknown agg.stat", {
  expect_error(dda.bagging(base.indep, data = d, iter = 5, agg.stat = "geometric"),
               regexp = "should be one of")
})

## ============================================================================
## 2. Return structure
## ============================================================================

test_that("classes follow the base DDA object", {
  expect_s3_class(bag.indep,   c("dda.bagging.indep", "dda.bagging"), exact = TRUE)
  expect_s3_class(bag.resdist, c("dda.bagging.resdist", "dda.bagging"), exact = TRUE)
  expect_s3_class(bag.vardist, c("dda.bagging.vardist", "dda.bagging"), exact = TRUE)
})

test_that("output holds all top-level elements", {
  expect_named(bag.indep,
               c("bagged.results", "aggregated.stats", "decisions",
                 "decision.proportions", "ols", "n.valid.iterations",
                 "alpha", "agg.stat", "trim.prob", "win.prob"),
               ignore.order = TRUE)
})

test_that("settings are stored in the object", {
  res <- run.bag(base.indep, d, iter = 5, alpha = 0.10,
                 agg.stat = "winsorized", trim.prob = 0.2, win.prob = 0.05)
  expect_equal(res$alpha, 0.10)
  expect_equal(res$agg.stat, "winsorized")
  expect_equal(res$trim.prob, 0.2)
  expect_equal(res$win.prob, 0.05)
})

test_that("n.valid.iterations matches the stored bootstrap results", {
  expect_true(bag.indep$n.valid.iterations > 0)
  expect_true(bag.indep$n.valid.iterations <= 10)
  expect_length(bag.indep$bagged.results, bag.indep$n.valid.iterations)
  expect_equal(nrow(bag.indep$decisions), bag.indep$n.valid.iterations)
})

test_that("agg.stat is stored for every method", {
  for (method in c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey")) {
    res <- run.bag(base.indep, d, iter = 5, agg.stat = method)
    expect_equal(res$agg.stat, method, info = method)
  }
})

## ============================================================================
## 3. Aggregated statistics
## ============================================================================

test_that("aggregated.stats keeps the class of the base object", {
  expect_s3_class(bag.indep$aggregated.stats, "dda.indep")
  expect_s3_class(bag.resdist$aggregated.stats, "dda.resdist")
  expect_s3_class(bag.vardist$aggregated.stats, "dda.vardist")
  expect_equal(bag.indep$aggregated.stats$var.names, c("y", "x"))
})

test_that("indep aggregates hold HSIC, dCor, BP, nlcor and difference results", {
  agg <- bag.indep$aggregated.stats
  expect_true(is.numeric(agg$hsic.yx$statistic))
  expect_true(is.numeric(agg$hsic.xy$p.value))
  expect_true(is.numeric(agg$distance.cor.dcor.yx$statistic))
  expect_true(is.numeric(agg$distance.cor.dcor.xy$p.value))
  expect_length(agg$breusch.pagan, 4)
  expect_false(is.null(agg$nlcor.yx$t1))
  expect_false(is.null(agg$nlcor.xy$t3))
  expect_equal(dim(agg$out.diff), dim(base.indep.full$out.diff))
})

test_that("resdist and vardist aggregates hold the separate tests", {
  expect_true(is.numeric(bag.resdist$aggregated.stats$agostino$target$statistic))
  expect_true(is.numeric(bag.resdist$aggregated.stats$anscombe$alternative$p.value))
  expect_true(is.numeric(bag.vardist$aggregated.stats$agostino$predictor$statistic))
  expect_true(is.numeric(bag.vardist$aggregated.stats$anscombe$outcome$p.value))
})

test_that("resdist boot.warning follows the aggregated kurtosis signs", {
  agg <- bag.resdist$aggregated.stats
  expect_equal(agg$boot.warning,
               sign(agg$anscombe$alternative$statistic[1]) != sign(agg$anscombe$target$statistic[1]))
})

test_that("trimmed with trim.prob = 0 and winsorized with win.prob = 0 equal the mean", {
  set.seed(101); res.mean <- run.bag(base.indep, d, iter = 5)
  set.seed(101); res.trim <- run.bag(base.indep, d, iter = 5, agg.stat = "trimmed", trim.prob = 0)
  set.seed(101); res.win  <- run.bag(base.indep, d, iter = 5, agg.stat = "winsorized", win.prob = 0)
  expect_equal(res.mean$aggregated.stats$hsic.yx$statistic,
               res.trim$aggregated.stats$hsic.yx$statistic)
  expect_equal(res.mean$aggregated.stats$hsic.yx$statistic,
               res.win$aggregated.stats$hsic.yx$statistic)
})

test_that("agg.value drops non-finite values and handles every method", {
  x <- c(1, 2, 3, 4, 100, NA, Inf)
  expect_equal(agg.value(x, "mean"), mean(c(1, 2, 3, 4, 100)))
  expect_equal(agg.value(x, "median"), 3)
  expect_equal(agg.value(x, "midhinge"), mean(quantile(c(1, 2, 3, 4, 100), c(0.25, 0.75), names = FALSE)))
  expect_true(is.numeric(agg.value(x, "tukey")))
  expect_true(is.na(agg.value(c(NA, Inf), "mean")))
  expect_error(agg.value(x, "geometric"), regexp = "Unknown agg.stat")
})

test_that("p-values are pooled with the harmonic mean p-value", {
  p <- c()
  for (fit in bag.vardist$bagged.results) p <- c(p, fit$agostino$outcome$p.value)
  expect_equal(bag.vardist$aggregated.stats$agostino$outcome$p.value, agg.p(p))

  p <- c()
  for (fit in bag.resdist$bagged.results) p <- c(p, fit$skewdiff[[3]])
  expect_equal(bag.resdist$aggregated.stats$skewdiff[[3]], agg.p(p))

  p <- c()
  for (fit in bag.indep$bagged.results) p <- c(p, fit$nlcor.yx$t1[[4]])
  expect_equal(bag.indep$aggregated.stats$nlcor.yx$t1[[4]], agg.p(p))
})

test_that("agg.stat leaves the pooled p-values unchanged", {
  set.seed(101); res.mean <- run.bag(base.indep, d, iter = 5)
  set.seed(101); res.med  <- run.bag(base.indep, d, iter = 5, agg.stat = "median")
  expect_equal(res.mean$aggregated.stats$hsic.yx$p.value,
               res.med$aggregated.stats$hsic.yx$p.value)
})

test_that("agg.p drops non-finite values", {
  expect_true(is.na(agg.p(c(NA, Inf))))
  expect_equal(agg.p(c(0.01, 0.2, NA)), agg.p(c(0.01, 0.2)))
})

test_that("agg.p returns zero when a p-value is zero", {
  expect_equal(agg.p(c(0, 0.2, 0.5)), 0)
})

## ============================================================================
## 4. Decisions
## ============================================================================

test_that("decision.proportions has one row per test and rows sum to 1", {
  props <- bag.indep$decision.proportions
  expect_equal(rownames(props), colnames(bag.indep$decisions))
  expect_equal(colnames(props), c("Target", "Alternative", "Confounding", "Undecided"))
  expect_equal(unname(rowSums(props)), rep(1, nrow(props)), tolerance = 1e-8)
})

test_that("resdist and vardist decisions have no Confounding level", {
  expect_equal(colnames(bag.resdist$decision.proportions), c("Target", "Alternative", "Undecided"))
  expect_equal(colnames(bag.vardist$decision.proportions), c("Target", "Alternative", "Undecided"))
})

test_that("dda.decisions returns the expected tests for each object type", {
  expect_setequal(names(dda.decisions(base.indep.full)),
                  c("hsic", "dcor", "bp", "nlcor", "hsic.diff", "dcor.diff", "mi.diff"))
  expect_true(all(c("agostino", "anscombe", "skewdiff", "kurtdiff") %in%
                    names(dda.decisions(base.vardist))))
  expect_true(all(dda.decisions(base.resdist) %in% c("Target", "Alternative", "Undecided", NA)))
})

test_that("dda.decisions swaps the separate resdist tests under prob.trans = TRUE", {
  obj <- base.resdist.pt
  obj$agostino$target$p.value      <- 0.01
  obj$agostino$alternative$p.value <- 0.50
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Target")
  obj$probtrans <- FALSE
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Alternative")
})

test_that("both significant separate tests are confounding for indep and undecided otherwise", {
  obj <- base.indep
  obj$hsic.yx$p.value <- 0.01
  obj$hsic.xy$p.value <- 0.01
  expect_equal(unname(dda.decisions(obj)["hsic"]), "Confounding")
  obj <- base.vardist
  obj$agostino$outcome$p.value   <- 0.01
  obj$agostino$predictor$p.value <- 0.01
  expect_equal(unname(dda.decisions(obj)["agostino"]), "Undecided")
})

test_that("resdist difference tests decide from p and the sign of the difference", {
  obj <- base.resdist
  obj$skewdiff[1] <- 0.4
  obj$skewdiff["p.value"] <- 0.01
  expect_equal(unname(dda.decisions(obj)["skewdiff.test"]), "Target")
  obj$skewdiff[1] <- -0.4
  expect_equal(unname(dda.decisions(obj)["skewdiff.test"]), "Alternative")
  obj$skewdiff["p.value"] <- 0.20
  expect_equal(unname(dda.decisions(obj)["skewdiff.test"]), "Undecided")
})

test_that("non-linear correlation decisions use alpha / 3", {
  obj <- base.indep.full
  for (k in c("t1", "t2", "t3")) {
    obj$nlcor.yx[[k]][4] <- 0.50
    obj$nlcor.xy[[k]][4] <- 0.03
  }
  expect_equal(unname(dda.decisions(obj, alpha = 0.05)["nlcor"]), "Undecided")
  for (k in c("t1", "t2", "t3")) obj$nlcor.xy[[k]][4] <- 0.01
  expect_equal(unname(dda.decisions(obj, alpha = 0.05)["nlcor"]), "Target")
})

test_that("dda.decisions rejects other objects", {
  expect_error(dda.decisions(list(a = 1)), regexp = "must be a dda.indep")
})

## ============================================================================
## 5. save.file, covariates, prob.trans
## ============================================================================

test_that("save.file writes the dda.bagging object to disk", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  run.bag(base.indep, d, iter = 5, save.file = tmp)
  saved <- readRDS(tmp)
  expect_s3_class(saved, "dda.bagging.indep")
})

test_that("dda.bagging works with covariates", {
  res <- run.bag(base.indep.cov, d.cov, iter = 5)
  expect_s3_class(res, "dda.bagging.indep")
  expect_true("z" %in% colnames(res$ols$target$coef))
  expect_true("z" %in% colnames(res$ols$alternative$coef))
})

test_that("dda.bagging handles prob.trans = TRUE", {
  res <- run.bag(base.resdist.pt, d, iter = 5)
  expect_true(isTRUE(res$aggregated.stats$probtrans))
  expect_true("agostino" %in% rownames(res$decision.proportions))
})

## ============================================================================
## 6. Failed bootstrap samples
## ============================================================================

test_that("failed bootstrap samples are dropped with a warning", {
  n.calls <- 0
  vardist <- dda.vardist
  local_mocked_bindings(dda.vardist = function(...) {
    n.calls <<- n.calls + 1
    if (n.calls %% 2 == 0) stop("planned failure")
    vardist(...)
  })
  expect_warning(res <- dda.bagging(base.vardist, data = d, iter = 4, progress = FALSE),
                 regexp = "2 of 4 bootstrap samples failed and were dropped. Last error: planned failure")
  expect_equal(res$n.valid.iterations, 2)
  expect_equal(nrow(res$decisions), 2)
  expect_equal(nrow(res$ols$target$coef), 2)
})

test_that("dda.bagging reports the error when every bootstrap sample fails", {
  expect_error(dda.bagging(base.vardist, data = data.frame(y = d$y, w = d$x), iter = 2, progress = FALSE),
               regexp = "failed in every bootstrap sample: .*not found")
})

## ============================================================================
## 7. inner.B and parallelize
## ============================================================================

test_that("inner.B sets B of the DDA call within each bootstrap sample", {
  res <- run.bag(base.vardist, d, iter = 3, inner.B = 15)
  expect_equal(as.numeric(res$bagged.results[[1]]$boot.args[3]), 15)
})

test_that("arguments of the original call are evaluated when bagging starts", {
  B.local <- 12
  base <- dda.vardist(y ~ x, pred = "x", data = d, B = B.local)
  res <- dda.bagging(base, data = d, iter = 2, progress = FALSE)
  expect_equal(as.numeric(res$bagged.results[[1]]$boot.args[3]), 12)
})

test_that("parallel bagging returns the same structure as serial bagging", {
  skip_on_cran()
  set.seed(1)
  res <- run.bag(base.vardist, d, iter = 4, parallelize = TRUE, cores = 2)
  expect_s3_class(res, "dda.bagging.vardist")
  expect_equal(res$n.valid.iterations, 4)
  expect_equal(nrow(res$ols$target$coef), 4)
})
