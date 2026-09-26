#' @title OLS Summary of Bootstrap Aggregated DDA Objects
#'
#' @description \code{summary_ols} returns aggregated ordinary least squares
#' (OLS) regression summaries for the causally competing target and
#' alternative models of a bootstrap aggregated DDA object. Coefficients are
#' aggregated across bootstrap samples with the method used in
#' \code{dda.bagging} or the method given in \code{agg_stat}.
#'
#' @param object An object of class \code{dda_bagging} obtained from
#'   \code{dda.bagging}.
#' @param agg_stat Character. Method used to aggregate coefficients and
#'   R-squared values across bootstrap samples. Must be one of
#'   \code{c("mean", "median", "trimmed", "winsorized", "midhinge",
#'   "tukey")}. If \code{NULL}, the method used in \code{dda.bagging} is
#'   kept.
#' @param trim_prob Numeric. Proportion of observations trimmed from each
#'   side of the sampling distribution when \code{agg_stat = "trimmed"}.
#' @param win_prob Numeric. Proportion of observations winsorized on each
#'   side of the sampling distribution when \code{agg_stat = "winsorized"}.
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
#' base_model <- dda.vardist(y ~ x, pred = "x", data = d, B = 10)
#' bagged <- dda.bagging(base_model, data = d, iter = 10, progress = FALSE)
#'
#' summary_ols(bagged)
#' summary_ols(bagged, agg_stat = "median")
#'
#' @export
summary_ols <- function(object,
                        agg_stat = NULL,
                        trim_prob = object$trim_prob,
                        win_prob = object$win_prob,
                        digits = 4,
                        ...){

  if (!inherits(object, "dda_bagging")) stop("object must be a dda_bagging object.")
  if (is.null(agg_stat)) agg_stat <- object$agg_stat

  alpha <- object$alpha
  output <- list()

  cat("\n")
  cat(paste("Aggregation method:", agg_stat), "\n")
  cat(paste("Number of bootstrap samples:", object$n_valid_iterations), "\n")

  for (model in c("target", "alternative")) {

    coefs <- object$ols[[model]]$coef
    pvals <- object$ols[[model]]$p.value
    r2    <- object$ols[[model]]$r.squared

    tab <- matrix(NA, nrow = ncol(coefs), ncol = 4)
    for (j in 1:ncol(coefs)) {
      tab[j, 1]   <- agg.value(coefs[, j], agg_stat, trim_prob, win_prob)
      tab[j, 2:3] <- quantile(coefs[, j], probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
      tab[j, 4]   <- mean(pvals[, j] < alpha, na.rm = TRUE)
    }
    rownames(tab) <- colnames(coefs)
    colnames(tab) <- c("estimate", paste(100 * alpha / 2, "%"), paste(100 * (1 - alpha / 2), "%"),
                       paste("Prop(p < ", alpha, ")", sep = ""))

    rsq     <- agg.value(r2[, "r.squared"], agg_stat, trim_prob, win_prob)
    adj.rsq <- agg.value(r2[, "adj.r.squared"], agg_stat, trim_prob, win_prob)

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
