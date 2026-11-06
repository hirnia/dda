#' Robust Breusch-Pagan Test
#'
#' Breusch-Pagan test for heteroskedasticity based on the residuals of a fitted model (e.g., Siegel repeated median regression).
#' @param fit A fitted model with \code{residuals} and \code{model.matrix} methods.
#' @param studentize Logical, if TRUE computes Koenker's studentized statistic, else the original Breusch-Pagan statistic.
#' @param ... Further arguments, ignored.
#' @return A list with test statistic, df, p-value, and method.
#'
#' @noRd
bptestrobust <- function(fit, studentize = FALSE, ...) {
  res <- residuals(fit)
  X <- model.matrix(fit)
  # Auxiliary regression: regress squared residuals on regressors
  aux.fit <- lm(res^2 ~ X[, -1, drop = FALSE])
  if (studentize) {
    bp.stat <- length(res) * summary(aux.fit)$r.squared
  } else {
    bp.stat <- sum((fitted(aux.fit) - mean(res^2))^2) / (2 * mean(res^2)^2)
  }
  df <- ncol(X) - 1
  pval <- pchisq(bp.stat, df, lower.tail = FALSE)
  structure(list(
    statistic = bp.stat,
    parameter = df,
    p.value = pval,
    method = sprintf("Robust Breusch-Pagan test (studentize=%s)", studentize)
  ), class = "htest")
}
