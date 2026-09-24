#' @title Direction Dependence Analysis: Residual Distributions
#'
#' @description \code{dda.resdist} evaluates patterns of asymmetry of error
#' distributions of causally competing models (\code{y ~ x} vs. \code{x ~ y}).
#' \code{print} returns DDA test statistics associated with \code{dda.resdist}
#' objects.
#'
#' @param formula Symbolic formula of the target model to be tested or an
#'   \code{lm} object.
#' @param pred Variable name of the predictor which serves as the outcome in
#'   the alternative model.
#' @param data An optional data frame containing the variables in the model
#'   (by default variables are taken from the environment which
#'   \code{dda.resdist} is called from).
#' @param B Number of bootstrap samples.
#' @param boot.type A vector of character strings representing the type of
#'   bootstrap confidence intervals. Must be one of the two values
#'   \code{c("perc", "bca")}; \code{boot.type = "perc"} is the default.
#' @param prob.trans A logical value indicating whether a probability integral
#'   transformation should be performed prior computation of skewness and
#'   kurtosis tests.
#' @param conf.level Confidence level for bootstrap confidence intervals.
#' @param robust A logical value indicating whether Siegel's (1982) repeated
#'   median estimator should be used for model estimation. If
#'   \code{robust = TRUE} repeated median estimation is applied for the
#'   causally competing models, otherwise ordinary least squares (OLS)
#'   estimation is used.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return An object of class \code{dda.resdist} containing the results of
#'   direction dependence tests of error distributions.
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction
#'   Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK:
#'   Cambridge University Press.
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' result <- dda.resdist(y ~ x, pred = "x", data = d,
#'             B = 50, conf.level = 0.90, prob.trans = FALSE, robust = FALSE)
#'
#' print(result)
#'
#' ## --- prob.trans = TRUE maps the predictor onto the marginal of the outcome,
#' ## --- so it is informative when the outcome is itself non-normal.
#'
#' d.skew <- data.frame(x, y = 0.5 * x + (rchisq(n, df = 3) - 3))
#'
#' dda.resdist(y ~ x, pred = "x", data = d.skew,
#'             B = 50, conf.level = 0.90, prob.trans = TRUE)
#'
#' \dontrun{
#'
#' result <- dda.resdist(y ~ x, pred = "x", data = d, B = 500,
#'   boot.type = "bca", conf.level = 0.95, prob.trans = FALSE)
#'
#' print(result)
#'
#' ## --- Siegel (1982) repeated median estimation
#'
#' dda.resdist(y ~ x, pred = "x", data = d, B = 500, robust = TRUE)
#' }
#'
#' @export
#' @rdname dda.resdist
dda.resdist <- function(formula,
                        pred = NULL,
                        data = list(),
                        B = 200,
                        boot.type = "perc",
                        prob.trans = FALSE,
                        conf.level = 0.95,
                        robust = FALSE
) {
  # --- Helper functions ---
  mysd <- function(x) sqrt(sum((x - mean(x))^2) / length(x))

  cor.ij <- function(x, y, i = 1, j = 1) {
    n <- length(x)
    mx <- mean(x)
    my <- mean(y)
    Cov <- sum((x - mx)^i * (y - my)^j) / n
    denom <- (mysd(x)^i * mysd(y)^j)
    if (denom == 0) return(NA_real_)
    Cov / denom
  }

  # Safe probability integral transform used in original code
  prob.int <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    fval <- ecdf(y)(y)
    # create monotone mapping of fval -> y by averaging y for same fval
    df <- data.frame(f = fval, y = y)
    df2 <- aggregate(y ~ f, data = df, FUN = mean)
    quant_y <- approxfun(df2$f, df2$y, rule = 2)
    cum_x <- ecdf(x)
    x_trans <- quant_y(cum_x(x))
    # if outside range -> use min/max
    x_trans[is.na(x_trans) & cum_x(x) <= min(df2$f, na.rm=TRUE)] <- min(df2$y, na.rm = TRUE)
    x_trans[is.na(x_trans) & cum_x(x) >= max(df2$f, na.rm=TRUE)] <- max(df2$y, na.rm = TRUE)
    list(x = x_trans, y = y)
  }

  # --- Bootstrap statistic function (resamples the residual pairs) ---
  boot.diff <- function(dat, indices, prob.trans) {
    d <- dat[indices, , drop = FALSE]
    # choose whether to use transformed residuals (if present) or the raw residuals
    if (isTRUE(prob.trans) &&
        all(c("alternative.trans", "target.trans") %in% names(d)) &&
        !all(is.na(d$alternative.trans)) && !all(is.na(d$target.trans))) {
      x <- as.vector(scale(d$alternative.trans))
      y <- as.vector(scale(d$target.trans))
    } else {
      x <- as.vector(scale(d$alternative))
      y <- as.vector(scale(d$target))
    }
    # Degenerate checks
    if (anyNA(x) || anyNA(y) || length(unique(x)) < 3 || length(unique(y)) < 3) return(rep(NA_real_, 7))

    # differences in marginal moments
    skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
    kurt.diff <- (moments::kurtosis(x) - 3)^2 - (moments::kurtosis(y) - 3)^2

    # joint statistics (mirroring original calculations)
    cor12.diff <- (cor.ij(y, x, 2, 1)^2) - (cor.ij(y, x, 1, 2)^2)
    cor13.diff <- ((cor.ij(y, x, 3, 1)^2) - (cor.ij(y, x, 1, 3)^2)) * sign(moments::kurtosis(y) - 3)

    xx <- sign(moments::skewness(x)) * x
    yy <- sign(moments::skewness(y)) * y

    safe_cor_xx_yy <- tryCatch(cor(xx, yy), error = function(e) NA_real_)
    safe_cor_x_y <- tryCatch(cor(x, y), error = function(e) NA_real_)

    RHS3 <- NA_real_
    if (!is.na(safe_cor_xx_yy)) {
      RHS3 <- safe_cor_xx_yy * mean((yy^2 * xx) - (yy * xx^2))
    }
    RHS4 <- NA_real_
    if (!is.na(safe_cor_x_y)) {
      RHS4 <- safe_cor_x_y * mean((y^3 * x) - (y * x^3)) * sign(moments::kurtosis(y) - 3)
    }
    C1 <- mean(y^3 * x) - 3 * safe_cor_x_y * var(y)
    C2 <- mean(y * x^3) - 3 * safe_cor_x_y * var(x)
    RCC <- (C1 + C2) * (C1 - C2)

    c(skew.diff, kurt.diff, cor12.diff, cor13.diff, RHS3, RCC, RHS4)
  }

  # --- Normality test helpers ---
  skew.diff.test <- function(x, y){
    agostino.zvalue <- function(x){
      n  <- length(x)
      s3 <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
      y <- s3 * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
      b2 <- 3 * (n^2 + 27*n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
      w <- -1 + sqrt(2 * (b2 - 1))
      d <- 1/sqrt(log(sqrt(w)))
      a <- sqrt(2/(w - 1))
      z <- d * log(y/a + sqrt((y/a)^2 + 1))
      return(z)
    }
    zval <- (agostino.zvalue(x) - agostino.zvalue(y))/sqrt(2 - 2*cor(x,y)^3)
    pval <- (1 - pnorm(abs(zval))) * 2
    return(list(z.value = zval, p.value = pval))
  }

  kurt.diff.test <- function(x, y){
    anscombe.zvalue <- function(x){
      n   <- length(x)
      b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
      eb2 <- 3 * (n - 1)/(n + 1)
      vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
      m3  <- (6 * (n^2 - 5*n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
      a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
      xx  <- (b - eb2)/sqrt(vb2)
      cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
      z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))
      return(z)
    }
    zval <- (anscombe.zvalue(x) - anscombe.zvalue(y))/sqrt(2 - 2*cor(x,y)^4)
    pval <- (1 - pnorm(abs(zval))) * 2
    return(list(z.value = zval, p.value = pval))
  }

  myanscombe.test <- function(x, alternative = c("two.sided", "less", "greater")){
    n   <- length(x)
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    s <- match.arg(alternative)
    alter <- switch(s, two.sided = 0, less = 1, greater = 2)
    b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    eb2 <- 3 * (n - 1)/(n + 1)
    vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
    m3  <- (6 * (n^2 - 5*n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
    a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
    xx  <- (b - eb2)/sqrt(vb2)
    cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
    z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))
    pval <- pnorm(z, lower.tail = FALSE)
    if (alter == 0) { pval <- 2 * pval; if (pval > 1) pval <- 2 - pval; alt <- "kurtosis is not equal to 3" }
    else if (alter == 1) { alt <- "kurtosis is greater than 3" }
    else { pval <- 1 - pval; alt <- "kurtosis is lower than 3" }
    RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval,
                 alternative = alt, method = "Anscombe-Glynn kurtosis test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }

  # --- Input checks ---
  if (is.null(pred)) stop("Tentative predictor is missing.")
  if (!is.numeric(B) || B < 0) stop("Number of resamples 'B' must be non-negative.")
  if (!is.numeric(conf.level) || conf.level < 0 || conf.level > 1) stop("'conf.level' must be between 0 and 1")
  if (!boot.type %in% c("bca", "perc")) stop("Unknown argument in boot.type; choose 'bca' or 'perc'.")

  # --- Prepare outcome, predictor, and model matrix for covariates ---
  if (!inherits(formula, "formula")) {
    # accept a fitted lm object (mirrors dda.indep)
    X <- if (is.matrix(formula$x)) formula$x else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y)) formula$y else model.response(model.frame(formula))
    delete.pred <- which(colnames(X) == pred)
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")
    x <- X[, delete.pred]
    X <- X[, -delete.pred, drop = FALSE]
    if (!is.matrix(X)) X <- as.matrix(X)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    delete.pred <- which(colnames(X) == pred)
    if (length(delete.pred) == 0) stop("Specified predictor not found in the target model.")
    x <- X[, delete.pred]
    X <- X[, -delete.pred, drop = FALSE]
    if (!is.matrix(X)) X <- as.matrix(X)
  }

  # Compute residuals of y ~ covariates and x ~ covariates
  ry <- tryCatch(as.vector(scale(lm.fit(X, y)$residuals)), error = function(e) stop("Failed to compute residuals for response."))
  rx <- tryCatch(as.vector(scale(lm.fit(X, x)$residuals)), error = function(e) stop("Failed to compute residuals for predictor."))

  resid_df <- data.frame(rx = rx, ry = ry)

  # fit target and alternative regressions (on residualized variables)
  # Uses 'RobustLinearReg' namespace if robust=TRUE, assume package is installed/loaded or import handled
  tar <- tryCatch(if (robust) RobustLinearReg::siegel_regression(ry ~ rx, data = resid_df) else lm(ry ~ rx, data = resid_df), error = function(e) NULL)
  alt <- tryCatch(if (robust) RobustLinearReg::siegel_regression(rx ~ ry, data = resid_df) else lm(rx ~ ry, data = resid_df), error = function(e) NULL)
  if (is.null(tar) || is.null(alt)) stop("Model fitting failed on residuals. Consider using robust = FALSE or checking data.")

  dat <- data.frame(alternative = as.vector(scale(resid(alt))),
                    target = as.vector(scale(resid(tar))),
                    pred.adj = rx,
                    out.adj = ry)

  if (isTRUE(prob.trans)) {
    # perform probability integral transform on the adjusted predictor & outcome
    trans <- prob.int(rx, ry)
    rx.trans <- trans$x
    ry.trans <- trans$y
    tar.trans <- tryCatch(if (robust) RobustLinearReg::siegel_regression(ry.trans ~ rx.trans, data = data.frame(ry.trans = ry.trans, rx.trans = rx.trans)) else lm(ry.trans ~ rx.trans), error = function(e) NULL)
    alt.trans <- tryCatch(if (robust) RobustLinearReg::siegel_regression(rx.trans ~ ry.trans, data = data.frame(rx.trans = rx.trans, ry.trans = ry.trans)) else lm(rx.trans ~ ry.trans), error = function(e) NULL)
    dat$alternative.trans <- if (!is.null(alt.trans)) as.vector(scale(resid(alt.trans))) else rep(NA_real_, nrow(dat))
    dat$target.trans      <- if (!is.null(tar.trans)) as.vector(scale(resid(tar.trans))) else rep(NA_real_, nrow(dat))
  }

  # --- Run separate normality tests ---
  if (isTRUE(prob.trans)) {
    agostino.out <- apply(dat[c("alternative.trans", "target.trans")], 2, moments::agostino.test)
    names(agostino.out) <- c("alternative", "target")
    agostino.out <- lapply(agostino.out, unclass)
    anscombe.out <- apply(dat[c("alternative.trans", "target.trans")], 2, myanscombe.test)
    names(anscombe.out) <- c("alternative", "target")
    anscombe.out <- lapply(anscombe.out, unclass)
  } else {
    agostino.out <- apply(dat[c("alternative", "target")], 2, moments::agostino.test)
    agostino.out <- lapply(agostino.out, unclass)
    anscombe.out <- apply(dat[c("alternative", "target")], 2, myanscombe.test)
    anscombe.out <- lapply(anscombe.out, unclass)
  }
  # prune some fields to match original output
  agostino.out$alternative[3:5] <- NULL
  agostino.out$target[3:5] <- NULL
  anscombe.out$alternative[3:5] <- NULL
  anscombe.out$target[3:5] <- NULL

  output <- list(agostino.out, anscombe.out)
  names(output) <- c("agostino", "anscombe")
  # adjust statistic to be kurtosis excess if needed, though anscombe test already returns kurtosis.
  # original code does -3 adjustment
  output$anscombe$target$statistic[1] <- output$anscombe$target$statistic[1] - 3
  output$anscombe$alternative$statistic[1] <- output$anscombe$alternative$statistic[1] - 3

  # --- Run asymptotic difference tests ---
  if (isTRUE(prob.trans)) {
    output <- c(output,
                list(skewdiff = unlist(skew.diff.test(dat$alternative.trans, dat$target.trans))),
                list(kurtdiff = unlist(kurt.diff.test(dat$alternative.trans, dat$target.trans)))
    )
    output$skewdiff <- c(moments::skewness(dat$alternative.trans)^2 - moments::skewness(dat$target.trans)^2, output$skewdiff)
    output$kurtdiff <- c((moments::kurtosis(dat$alternative.trans)-3)^2 - (moments::kurtosis(dat$target.trans)-3)^2, output$kurtdiff)
  } else {
    output <- c(output,
                list(skewdiff = unlist(skew.diff.test(dat$alternative, dat$target))),
                list(kurtdiff = unlist(kurt.diff.test(dat$alternative, dat$target)))
    )
    output$skewdiff <- c(moments::skewness(dat$alternative)^2 - moments::skewness(dat$target)^2, output$skewdiff)
    output$kurtdiff <- c((moments::kurtosis(dat$alternative)-3)^2 - (moments::kurtosis(dat$target)-3)^2, output$kurtdiff)
  }

  # --- Run bootstrap confidence intervals ---
  boot.warning <- FALSE
  if (B > 0) {
    suppressWarnings({
      boot.res <- tryCatch(boot::boot(dat, boot.diff, R = B, prob.trans = prob.trans), error = function(e) NULL)
    })

    if (is.null(boot.res)) {
      boot.warning <- "Bootstrap failed entirely (boot::boot returned NULL)."
      output$boot.args <- c(boot.type, conf.level, B)
      output$boot.warning <- boot.warning
    } else {
      empinf_ok <- TRUE
      if (boot.type == "bca") {
        empinf <- tryCatch(boot::empinf(boot.res), error = function(e) NA)
        if (any(is.na(empinf))) {
          empinf_ok <- FALSE
          boot.warning <- "Bootstrap failed: empinf() contains NA. BCa CIs not available; computed percentiles when possible."
        }
      }

      stat_names <- c("skew.diff", "kurt.diff", "cor12.diff", "cor13.diff", "RHS3", "RCC", "RHS4")
      boot.out <- vector("list", length = length(stat_names))
      names(boot.out) <- stat_names

      for (i in seq_along(stat_names)) {
        t_vec <- boot.res$t[, i]
        t0 <- boot.res$t0[i]
        t_non_na <- t_vec[!is.na(t_vec)]

        if (length(t_non_na) < 2 || sd(t_non_na) == 0 || is.na(t0)) {
          boot.out[[i]] <- list(type = boot.type, t0 = t0, lower = NA_real_, upper = NA_real_, t = t_vec)
          next
        }

        if (boot.type == "perc" || !empinf_ok) {
          alpha <- (1 - conf.level)/2
          lower <- as.numeric(stats::quantile(t_non_na, probs = alpha, na.rm = TRUE, names = FALSE))
          upper <- as.numeric(stats::quantile(t_non_na, probs = 1 - alpha, na.rm = TRUE, names = FALSE))
          boot.out[[i]] <- list(type = "perc", t0 = t0, lower = lower, upper = upper, t = t_vec)
        } else { # BCa requested and empinf OK
          ci_try <- tryCatch(boot::boot.ci(boot.res, conf = conf.level, type = "bca", t0 = t0, t = t_vec), error = function(e) NULL, warning = function(w) NULL)
          if (is.null(ci_try)) {
            alpha <- (1 - conf.level)/2
            lower <- as.numeric(stats::quantile(t_non_na, probs = alpha, na.rm = TRUE, names = FALSE))
            upper <- as.numeric(stats::quantile(t_non_na, probs = 1 - alpha, na.rm = TRUE, names = FALSE))
            boot.out[[i]] <- list(type = "perc_fallback", t0 = t0, lower = lower, upper = upper, t = t_vec)
          } else {
            ci_vals <- tryCatch({
              if (!is.null(ci_try$bca)) {
                as.numeric(ci_try$bca[4:5])
              } else if (!is.null(ci_try$percent)) {
                as.numeric(ci_try$percent[4:5])
              } else {
                rep(NA_real_, 2)
              }
            }, error = function(e) rep(NA_real_, 2))
            boot.out[[i]] <- list(type = "bca", t0 = t0, lower = ci_vals[1], upper = ci_vals[2], t = t_vec)
          }
        }
      }

      # Build the CI results
      ci.skewdiff  <- c(lower = boot.out$skew.diff$lower, upper = boot.out$skew.diff$upper)
      ci.kurtdiff  <- c(lower = boot.out$kurt.diff$lower, upper = boot.out$kurt.diff$upper)
      ci.cor12diff <- c(lower = boot.out$cor12.diff$lower, upper = boot.out$cor12.diff$upper)
      ci.cor13diff <- c(lower = boot.out$cor13.diff$lower, upper = boot.out$cor13.diff$upper)
      ci.RHS3      <- c(lower = boot.out$RHS3$lower, upper = boot.out$RHS3$upper)
      ci.RCC       <- c(lower = boot.out$RCC$lower, upper = boot.out$RCC$upper)
      ci.RHS4      <- c(lower = boot.out$RHS4$lower, upper = boot.out$RHS4$upper)

      output$skewdiff <- c(output$skewdiff, ci.skewdiff)
      output$kurtdiff <- c(output$kurtdiff, ci.kurtdiff)
      output <- c(output,
                  list(cor12diff = c(boot.res$t0[3], ci.cor12diff)),
                  list(cor13diff = c(boot.res$t0[4], ci.cor13diff)),
                  list(RHS3 = c(boot.res$t0[5], ci.RHS3)),
                  list(RCC = c(boot.res$t0[6], ci.RCC)),
                  list(RHS4 = c(boot.res$t0[7], ci.RHS4)),
                  list(boot.args = c(boot.type, conf.level, B)),
                  list(boot.warning = if (is.character(boot.warning)) boot.warning else FALSE)
      )
    }
  }

  response.name <- all.vars(formula(formula))[1]
  output <- c(output, list(var.names = c(response.name, pred), probtrans = prob.trans))

  # --- CRITICAL FIX: Add call_info so dda_bagging can find data ---
  call_info <- list(
    "function_call" = match.call(),
    "function_name" = "dda.resdist",
    "all_args" = as.list(match.call())[-1],
    "formula" = formula,
    "data_name" = deparse(substitute(data)),
    "original_data" = if(missing(data) || is.null(data)) NULL else data
  )
  output <- c(output, list(call_info = call_info))

  class(output) <- "dda.resdist"
  return(output)
}

#' @title Print Method for \code{dda.resdist} Objects
#'
#' @param x An object of class \code{dda.resdist} when using \code{print}.
#' @param ... Additional arguments to be passed to the method.
#'
#' @examples
#' \dontrun{
#' print(result)
#' }
#'
#' @export
#' @rdname dda.resdist
#' @method print dda.resdist
print.dda.resdist <- function(x, ...){
  object <- x

  varnames <- object$var.names
  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n")
  cat("Skewness and kurtosis tests:", "\n")

  sigtests <- rbind( c(object[[1]]$target$statistic, object[[1]]$target$p.value, object[[1]]$alternative$statistic, object[[1]]$alternative$p.value),
                     c(object[[2]]$target$statistic, object[[2]]$target$p.value, object[[2]]$alternative$statistic, object[[2]]$alternative$p.value)
  )
  sigtests <- round(sigtests, 4)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
  print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  if(is.null(object$boot.args)){
    cat("\n")
    cat("Skewness and kurtosis difference tests:", "\n")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  if(!is.null(object$boot.args)){
    ci.level <- as.numeric(object$boot.args[2]) * 100
    cat("\n")

    # Print skewness and kurtosis results

    if(object$boot.args[1] == "bca")  cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "BCa bootstrap CIs:", "\n", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "Percentile bootstrap CIs:", "\n", "\n", sep = "")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")

    # Print co-skewness and co-kurtosis results

    if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for joint higher moment differences:", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for joint higher moment differences:", "\n", sep = "")

    jointtests <- rbind(object$cor12diff, object$RHS3, object$cor13diff, object$RHS4, object$RCC)
    jointtests <- round(jointtests, 4)
    rownames(jointtests) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
    colnames(jointtests) <- c("estimate", "lower", "upper")
    print.default(format( jointtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    # hoctests <- rbind(object$cor12diff, object$cor13diff)
    # hoctests <- round(hoctests, 4)
    # rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]" )
    # colnames(hoctests) <- c("estimate", "lower", "upper")
    # print.default(format( hoctests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    # cat("\n")

    # # Print LR approximative measures

    # if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")
    # if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")

    # LRtests <- rbind(object$RHS3, object$RHS4, object$RCC )
    # LRtests <- round(LRtests, 4)
    # rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (co-kurtosis)", "Chen-Chan (co-kurtosis)")
    # colnames(LRtests) <- c("estimate", "lower", "upper")
    # print.default(format( LRtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    cat("\n")
    cat(paste("Number of resamples:", object$boot.args[3]))
    cat("\n")
  }

  cat("---")
  cat("\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
  cat("\n")
  if(isFALSE(object$probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], sep = " "))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and", "\n", "     co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], sep = " "))
  }
  cat("\n")
  if(object$boot.warning) { cat("Warning: Excess-kurtosis values of residuals have unequal signs", "\n", "        Also compute Co-Kurtosis and Hyvarinen-Smith Co-Kurtosis for", varnames[1], "->", varnames[2], "\n") }
	cat("\n")
  }


