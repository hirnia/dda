#' @title Summary Method for Bootstrap Aggregated DDA Objects
#'
#' @description \code{summary} returns the proportion of model selection
#' decisions (target model, alternative model, confounding, or undecided)
#' across the bootstrap samples of a bootstrap aggregated Direction
#' Dependence Analysis (DDA) object.
#'
#' @param object An object of class \code{dda_bagging} obtained from
#'   \code{dda.bagging}.
#' @param show Character vector specifying the DDA statistics to report.
#'   Accepts \code{c("hsic", "dcor", "mi", "bp", "nlcor")} for independence
#'   properties and \code{c("skew", "kurt", "coskew", "cokurt")} for variable
#'   and residual distributions. If \code{NULL} (default), all available
#'   statistics are reported.
#' @param digits Integer. Number of decimal places for proportions
#'   (default: 2).
#' @param ... Additional arguments to be passed to the function.
#'
#' @details Confounding is a possible decision only for the separate
#'   independence tests of \code{dda.indep}. The decision rules are listed in
#'   the details of \code{\link{dda.bagging}}.
#'
#' @return Invisibly returns the table of decision proportions.
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
#' ## proportion of model selection decisions across bootstrap samples
#' summary(bagged)
#'
#' ## third moment statistics only
#' summary(bagged, show = c("skew", "coskew"))
#'
#' @export
#' @rdname summary.dda_bagging
#' @method summary dda_bagging
summary.dda_bagging <- function(object, show = NULL, digits = 2, ...){

  props <- object$decision_proportions
  varnames <- object$aggregated_stats$var.names

  groups <- list(hsic   = c("hsic", "hsic.diff"),
                 dcor   = c("dcor", "dcor.diff"),
                 mi     = "mi.diff",
                 bp     = "bp",
                 nlcor  = "nlcor",
                 skew   = c("agostino", "skewdiff"),
                 kurt   = c("anscombe", "kurtdiff"),
                 coskew = c("cor12diff", "RHS", "RHS3"),
                 cokurt = c("cor13diff", "RCC", "RHS4", "Rtanh"))

  labels <- c(hsic      = "HSIC",
              dcor      = "dCor",
              bp        = "Robust Breusch-Pagan",
              nlcor     = "Non-linear Correlation",
              hsic.diff = "HSIC Difference",
              dcor.diff = "dCor Difference",
              mi.diff   = "MI Difference",
              agostino  = "Separate D'Agostino Tests",
              anscombe  = "Separate Anscombe-Glynn Tests",
              skewdiff  = "Skewness Difference",
              kurtdiff  = "Kurtosis Difference",
              cor12diff = "Co-Skewness Difference",
              cor13diff = "Co-Kurtosis Difference",
              RHS       = "Hyvarinen-Smith Co-Skewness Difference",
              RHS3      = "Hyvarinen-Smith Co-Skewness Difference",
              RHS4      = "Hyvarinen-Smith Co-Kurtosis Difference",
              RCC       = "Chen-Chan Co-Kurtosis Difference",
              Rtanh     = "Hyvarinen-Smith tanh Difference")

  if (!is.null(show)) {
    keep <- c()
    for (s in show) {
      if (!s %in% names(groups)) stop(paste("Unknown statistic in show:", s))
      keep <- c(keep, groups[[s]])
    }
    props <- props[rownames(props) %in% keep, , drop = FALSE]
  }

  if (nrow(props) == 0) stop("None of the requested statistics are available in this object.")

  out <- props
  for (r in 1:nrow(out)) out[r, ] <- round_preserve_sum(props[r, ], digits)
  rownames(out) <- labels[rownames(props)]

  if (inherits(object, "dda_bagging_indep"))   type <- "Independence Properties"
  if (inherits(object, "dda_bagging_resdist")) type <- "Residual Distributions"
  if (inherits(object, "dda_bagging_vardist")) type <- "Variable Distributions"

  cat("\n")
  cat(paste("BOOTSTRAP AGGREGATED DDA:", type), "\n")
  cat(paste("Number of bootstrap samples:", object$n_valid_iterations), "\n")
  cat(paste("Proportion of model selection decisions (alpha = ", object$alpha, "):", sep = ""), "\n", "\n")

  print.default(format(out, nsmall = digits), print.gap = 2L, quote = FALSE)

  cat("---")
  cat("\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
  cat("\n")

  invisible(props)
}


#' @title Round Proportions so They Sum to One
#'
#' @description Rounds a vector of proportions down to \code{digits}
#'   decimals and adds the missing units to the entries with the largest
#'   remainders (largest remainder method).
#'
#' @keywords internal
#' @noRd
round_preserve_sum <- function(x, digits = 2){

  if (any(is.na(x))) return(round(x, digits))

  scaled  <- x * 10^digits
  rounded <- floor(scaled)
  missing <- round(sum(scaled) - sum(rounded))

  add <- order(scaled - rounded, decreasing = TRUE)[seq_len(missing)]
  rounded[add] <- rounded[add] + 1

  rounded / 10^digits
}
