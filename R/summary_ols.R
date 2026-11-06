#' @title OLS Summary of Bootstrap Aggregated DDA Objects
#'
#' @description \code{ols.summary} returns aggregated ordinary least squares
#' (OLS) regression summaries for the causally competing target and
#' alternative models of a bootstrap aggregated DDA object. Coefficients are
#' aggregated across bootstrap samples with the method used in
#' \code{dda.bagging} or the method given in \code{agg.stat}.
#'
#' @param object An object of class \code{dda.bagging} obtained from
#'   \code{dda.bagging}.
#' @param agg.stat Character. Specifies the method used for aggregating
#'   coefficients and R-squared values across bootstrap samples. Must be one
#'   of the following specifications \code{c("mean", "median", "trimmed",
#'   "winsorized", "midhinge", "tukey")}. If \code{NULL} (default), the
#'   function uses the method applied within \code{dda.bagging}.
#' @param trim.prob Numeric. Proportion of observations to be trimmed from
#'   each side of the sampling distribution when \code{agg.stat = "trimmed"}
#'   (default: the value used in \code{dda.bagging}).
#' @param win.prob Numeric. Proportion of observations to be winsorized from
#'   each side of the sampling distribution when
#'   \code{agg.stat = "winsorized"} (default: the value used in
#'   \code{dda.bagging}).
#' @param digits Integer. Number of digits used for rounding (default: 4).
#' @param ... Additional arguments to be passed to the function.
#'
#' @details For each coefficient the table reports the aggregated estimate,
#'   the percentile interval of the estimates across bootstrap samples at the
#'   \code{alpha} level used in \code{dda.bagging}, and the proportion of
#'   bootstrap samples in which the coefficient is significant at that
#'   level.
#'
#' @return Invisibly returns a list with the coefficient tables and
#'   aggregated R-squared values of the target and alternative models.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' base.model <- dda.vardist(y ~ x, pred = "x", data = d, B = 10)
#' bagged <- dda.bagging(base.model, data = d, iter = 10, progress = FALSE)
#'
#' ols.summary(bagged)
#' ols.summary(bagged, agg.stat = "median")
#'
#' @export
ols.summary <- function(object,
                        agg.stat = NULL,
                        trim.prob = object$trim.prob,
                        win.prob = object$win.prob,
                        digits = 4,
                        ...){

  if (!inherits(object, "dda.bagging")) stop("object must be a dda.bagging object.")
  if (is.null(agg.stat)) agg.stat <- object$agg.stat

  alpha <- object$alpha
  output <- list()

  cat("\n")
  cat(paste("Aggregation method:", agg.stat), "\n")
  cat(paste("Number of bootstrap samples:", object$n.valid.iterations), "\n")

  for (model in c("target", "alternative")) {

    coefs <- object$ols[[model]]$coef
    pvals <- object$ols[[model]]$p.value
    r2    <- object$ols[[model]]$r.squared

    tab <- matrix(NA, nrow = ncol(coefs), ncol = 4)
    for (j in 1:ncol(coefs)) {
      tab[j, 1]   <- agg.value(coefs[, j], agg.stat, trim.prob, win.prob)
      tab[j, 2:3] <- quantile(coefs[, j], probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
      tab[j, 4]   <- mean(pvals[, j] < alpha, na.rm = TRUE)
    }
    rownames(tab) <- colnames(coefs)
    colnames(tab) <- c("estimate", paste(100 * alpha / 2, "%"), paste(100 * (1 - alpha / 2), "%"),
                       paste("Prop(p < ", alpha, ")", sep = ""))

    rsq     <- agg.value(r2[, "r.squared"], agg.stat, trim.prob, win.prob)
    adj.rsq <- agg.value(r2[, "adj.r.squared"], agg.stat, trim.prob, win.prob)

    if (model == "target")      cat("\n", "OLS Summary: Target Model", "\n", sep = "")
    if (model == "alternative") cat("\n", "OLS Summary: Alternative Model", "\n", sep = "")
    cat(deparse(object$ols[[model]]$formula), "\n", "\n")

    print.default(format(round(tab, digits), nsmall = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
    cat(paste("R-squared: ", round(rsq, digits), ", Adjusted R-squared: ", round(adj.rsq, digits), sep = ""), "\n")

    output[[model]] <- list(coefficients = tab, r.squared = c(r.squared = rsq, adj.r.squared = adj.rsq))
  }

  invisible(output)
}
