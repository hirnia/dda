## OLD
## #' @details For each regression coefficient the table reports the aggregated
## #'   point estimate, the bootstrap percentile confidence interval of the ...
## NEW 0.2.0
## Edited in the roxygen block below: @details in WW's wording.

#' @title OLS Summary of Bootstrap Aggregated DDA Objects
#'
#' @description \code{ols.summary} returns aggregated ordinary least squares
#' (OLS) regression summaries for the causally competing target and
#' alternative models of a bootstrap aggregated DDA object. Coefficients are
#' aggregated across bootstrap samples with the method used in
#' \code{dda.bagging} or the method given in \code{agg.stat}.
#' \code{print} returns the aggregated OLS summaries of \code{ols.summary}
#' objects.
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
#' @param x An object of class \code{ols.summary} when using \code{print}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @details For each coefficient the table reports the aggregated point
#'   estimates, the bootstrap percentile confidence intervals for the estimates
#'   (confidence level \code{1 - alpha}, with \code{alpha} as used in
#'   \code{dda.bagging}), and the proportion of bootstrap samples in which the
#'   regression coefficients are significant at the \code{alpha} level.
#'
#' @return An object of class \code{ols.summary} containing the coefficient
#'   tables and aggregated R-squared values of the target and alternative
#'   models:
#'   \item{target}{Formula, coefficient table, and aggregated R-squared values
#'     of the target model.}
#'   \item{alternative}{Formula, coefficient table, and aggregated R-squared
#'     values of the alternative model.}
#'   \item{agg.stat}{Aggregation method.}
#'   \item{n.valid.iterations}{Number of bootstrap samples with a DDA result.}
#'   \item{alpha}{Significance level used in \code{dda.bagging}.}
#'   \item{digits}{Number of digits used by \code{print}.}
#'
#' @seealso \code{\link{dda.bagging}}, \code{\link{summary.dda.bagging}}
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
#'
#' ## --- coefficient tables of the target model
#'
#' ols <- ols.summary(bagged, agg.stat = "median")
#' ols$target$coefficients
#'
#' @export
#' @rdname ols.summary
ols.summary <- function(object, ...){
  UseMethod("ols.summary")
}


#' @export
#' @rdname ols.summary
#' @method ols.summary default
ols.summary.default <- function(object, ...){
  stop("object must be a dda.bagging object.")
}


#' @export
#' @rdname ols.summary
#' @method ols.summary dda.bagging
ols.summary.dda.bagging <- function(object,
                                    agg.stat = NULL,
                                    trim.prob = object$trim.prob,
                                    win.prob = object$win.prob,
                                    digits = 4,
                                    ...){

  if (is.null(agg.stat)) agg.stat <- object$agg.stat

  alpha <- object$alpha
  output <- list()

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

    output[[model]] <- list(formula = object$ols[[model]]$formula,
                            coefficients = tab,
                            r.squared = c(r.squared = rsq, adj.r.squared = adj.rsq))
  }

  output$agg.stat <- agg.stat
  output$n.valid.iterations <- object$n.valid.iterations
  output$alpha <- alpha
  output$digits <- digits

  class(output) <- "ols.summary"
  return(output)
}


#' @export
#' @rdname ols.summary
#' @method print ols.summary
print.ols.summary <- function(x, digits = x$digits, ...){

  cat("\n")
  cat(paste("Aggregation method:", x$agg.stat), "\n")
  cat(paste("Number of bootstrap samples:", x$n.valid.iterations), "\n")

  for (model in c("target", "alternative")) {

    if (model == "target")      cat("\n", "OLS Summary: Target Model", "\n", sep = "")
    if (model == "alternative") cat("\n", "OLS Summary: Alternative Model", "\n", sep = "")
    cat(deparse(x[[model]]$formula), "\n", "\n")

    tab <- x[[model]]$coefficients
    print.default(format(round(tab, digits), nsmall = digits), print.gap = 2L, quote = FALSE)
    cat("\n")

    rsq <- x[[model]]$r.squared
    cat(paste("R-squared: ", round(rsq["r.squared"], digits), ", Adjusted R-squared: ", round(rsq["adj.r.squared"], digits), sep = ""), "\n")
  }

  invisible(x)
}
