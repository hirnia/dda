#' @title Print Summary for OLS Models of Bootstrap Aggregated DDA Objects
#'
#' @description \code{summary_ols} returns aggregated ordinary least
#' squares (OLS) regression summaries for the causally competing target and
#' alternative models of a bootstrap aggregated DDA object.
#' Regression coefficients and standard errors are aggregated across bootstrap
#' samples using the method specified in \code{dda.bagging()} or overridden
#' via \code{agg.stat}.
#'
#' @param object Output from \code{dda.bagging()}.
#' @param agg_stat Character. Specifies the method used for aggregating test
#'   statistics and coefficients across bootstrap samples. Must be one of the
#'   following specifications \code{c("mean", "median", "trimmed",
#'   "winsorized", "midhinge", "tukey")}. If \code{NULL}, the function
#'   uses the method applied with \code{dda.bagging()}.
#' @param trim_prob Numeric. Proportion of observations to be trimmed on each
#'   side of the sampling distribution when \code{agg.stat = "trimmed"}
#'   (default: 0.10).
#' @param win_prob Numeric. Proportion of observations to be winsorized on
#'   each side of the sampling distribution when
#'   \code{agg.stat = "winsorized"} (default: 0.10).
#' @param digits Integer. Number of digits used for rounding.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' base_model <- dda.indep(y ~ x, pred = "x", data = d, B = 20)
#' bagged <- dda.bagging(base_model, data = d, iter = 5,
#'                        progress = FALSE)
#'
#' # Print aggregated OLS coefficients for target and alternative models
#' summary_ols(bagged)
#'
#' \donttest{
#' # Override aggregation method
#' summary_ols(bagged, agg_stat = "median")
#' }
#' @export

summary_ols <- function(object,
                              agg_stat = NULL,
                              trim_prob = 0.10,
                              win_prob = 0.10,
                              digits = 4,
                              ...) {

  if (!inherits(object, "dda_bagging")) {
    stop("Object must be a bagged DDA result.")
  }

  # Ensure dynamic re-aggregation of coefficients occurs with split trim/win probs
  object <- reaggregate_bagging(object, agg_stat, trim_prob, win_prob)
  stats <- object$aggregated_stats
  raw <- object$raw_stats

  if (is.null(stats$ols_target)) {
    cat("No OLS summary available in this object.\n")
    return(invisible(NULL))
  }

  # Local helper for R-squared aggregation
  current_agg <- if (!is.null(agg_stat)) agg_stat else object$agg_stat_used

  agg_helper <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    switch(current_agg,
           "mean" = mean(x),
           "median" = median(x),
           "trimmed" = mean(x, trim = trim_prob),
           "winsorized" = {
             q_low <- quantile(x, probs = win_prob, na.rm = TRUE, names = FALSE)
             q_high <- quantile(x, probs = 1 - win_prob, na.rm = TRUE, names = FALSE)
             x[x < q_low] <- q_low
             x[x > q_high] <- q_high
             mean(x)
           },
           "midhinge" = {
             q <- quantile(x, probs = c(0.25, 0.75), names = FALSE)
             mean(q)
           },
           "tukey" = {
             q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
             (q[1] + 2*q[2] + q[3]) / 4
           }
    )
  }

  # Print Aggregation Method
  cat("\nAggregation method:", current_agg, "\n\n")

  # --- Target Model Print ---
  cat("OLS Summary: Target Model\n")
  # Fixed: Use round() for explicit decimal places instead of format()
  print.default(round(stats$ols_target, digits = digits), print.gap = 2L, quote = FALSE)

  if (!is.null(raw$ols_tar_rsq)) {
    r2 <- agg_helper(raw$ols_tar_rsq[, 1])
    adj_r2 <- agg_helper(raw$ols_tar_rsq[, 2])
    cat(sprintf("\nR-squared: %.*f, Adjusted R-squared: %.*f\n", digits, r2, digits, adj_r2))
  }

  cat("\n---\n\n")

  # --- Alternative Model Print ---
  cat("OLS Summary: Alternative Model\n")
  # Fixed: Use round() for explicit decimal places instead of format()
  print.default(round(stats$ols_alternative, digits = digits), print.gap = 2L, quote = FALSE)

  if (!is.null(raw$ols_alt_rsq)) {
    r2_alt <- agg_helper(raw$ols_alt_rsq[, 1])
    adj_r2_alt <- agg_helper(raw$ols_alt_rsq[, 2])
    cat(sprintf("\nR-squared: %.*f, Adjusted R-squared: %.*f\n", digits, r2_alt, digits, adj_r2_alt))
  }

  invisible(object)
}
