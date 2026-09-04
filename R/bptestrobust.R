#' Robust Breusch-Pagan Test
#'
#' Modified Breusch-Pagan test for heteroskedasticity, compatible with formula or fitted models from lm/mblm.
#' @param formula A regression formula (e.g., y ~ x).
#' @param data Optional: data.frame containing the variables in the formula.
#' @param studentize Logical, if TRUE uses squared studentized residuals, else uses squared raw residuals.
#' @param ... Further arguments, ignored.
#' @return A list with test statistic, df, p-value, and method.
#'
#' @noRd
bptestrobust <- function(formula, data = NULL, studentize = FALSE, ...) {
  # Handle formula or fitted model
  if (inherits(formula, "formula")) {
    mf <- model.frame(formula, data)
    y <- model.response(mf)
    X <- model.matrix(formula, mf)
    # Fit robust or OLS model
    if (!is.null(data) && ("mblm" %in% class(data))) {
      fit <- data
    } else {
      fit <- lm(formula, data = data)
    }
  } else if (inherits(formula, "lm") || inherits(formula, "mblm")) {
    fit <- formula
    y <- fit$model[[1]]
    X <- model.matrix(fit)
  } else {
    stop("First argument must be a formula or fitted model")
  }
  # Residuals
  res <- residuals(fit)
  if (studentize) {
    res <- res / sqrt(sum(res^2) / length(res))
  }
  # Auxiliary regression: regress squared residuals on regressors
  aux_fit <- lm(res^2 ~ X[, -1, drop = FALSE]) # Remove intercept
  bp_stat <- summary(aux_fit)$r.squared * length(res)
  df <- ncol(X) - 1
  pval <- pchisq(bp_stat, df, lower.tail = FALSE)
  structure(list(
    statistic = bp_stat,
    parameter = df,
    p.value = pval,
    method = sprintf("Robust Breusch-Pagan test (studentize=%s)", studentize)
  ), class = "htest")
}
