#' @title Bootstrap Aggregated Direction Dependence Analysis (DDA)
#'
#' @description
#' \code{dda.bagging} performs bootstrap aggregation (bagging) on an existing
#' Direction Dependence Analysis (DDA) object to evaluate the stability of
#' direction dependence decisions. \code{print} returns the aggregated DDA
#' test statistics.
#'
#' @param dda_result An output object from \code{dda.indep},
#'   \code{dda.resdist}, or \code{dda.vardist}.
#' @param data A \code{data.frame} containing all variables used in the
#'   original DDA model (outcome, predictor, and any covariates).
#' @param iter Number of bootstrap samples (default: 100).
#' @param alpha Numeric. Significance level used for causal model selection
#'   (default: 0.05).
#' @param agg_stat Character. Method used to aggregate test statistics,
#'   p-values, and coefficients across bootstrap samples. Must be one of
#'   \code{c("mean", "median", "trimmed", "winsorized", "midhinge",
#'   "tukey")}. In \code{print}, \code{NULL} (default) keeps the method
#'   used in \code{dda.bagging}.
#' @param trim_prob Numeric. Proportion of observations trimmed from each
#'   side of the sampling distribution when \code{agg_stat = "trimmed"}
#'   (default: 0.10).
#' @param win_prob Numeric. Proportion of observations winsorized on each
#'   side of the sampling distribution when \code{agg_stat = "winsorized"}
#'   (default: 0.10).
# inner_B is disabled for now. To restore it, uncomment this block, the
# argument in the signature and its use before the bootstrap loop.
# #' @param inner_B Optional positive integer. Number of inner bootstrap
# #'   resamples (\code{B}) passed to each DDA call. \code{NULL} (default)
# #'   keeps the \code{B} used in the original DDA call.
#' @param progress Logical. Whether to display a progress bar (default:
#'   \code{TRUE}).
#' @param save_file Character. Optional file path used to save the result as
#'   an R data file (e.g., \code{"results.rds"}).
#'
#' @details
#' In each of the \code{iter} bootstrap samples, \code{n} observations are
#' drawn with replacement from \code{data}, the original DDA function is
#' called again with its original arguments, and the causally competing OLS
#' models are fitted. Test statistics and p-values are aggregated across
#' bootstrap samples with \code{agg_stat}. An aggregated p-value summarizes
#' the distribution of p-values across bootstrap samples; it is not a pooled
#' test of a single null hypothesis.
#'
#' Within each bootstrap sample every DDA test leads to a model selection
#' decision. The target model is \code{x -> y} and the alternative model is
#' \code{y -> x}; \code{p_yx} is the p-value obtained under the target model
#' and \code{p_xy} the p-value obtained under the alternative model.
#'
#' Separate normality tests (\code{dda.vardist}: outcome and predictor;
#' \code{dda.resdist}: target and alternative residuals):
#' \itemize{
#'   \item \code{p_yx > alpha} and \code{p_xy <= alpha}: target model
#'   \item \code{p_yx <= alpha} and \code{p_xy > alpha}: alternative model
#'   \item otherwise: undecided
#' }
#' For \code{dda.resdist} with \code{prob.trans = TRUE} the target and
#' alternative residual p-values swap roles in these rules
#' (\code{p_yx <= alpha} and \code{p_xy > alpha}: target model).
#'
#' Separate independence tests (\code{dda.indep}: HSIC, dCor, robust
#' Breusch-Pagan, and non-linear correlation using the smallest of the three
#' p-values):
#' \itemize{
#'   \item \code{p_yx > alpha} and \code{p_xy <= alpha}: target model
#'   \item \code{p_yx <= alpha} and \code{p_xy > alpha}: alternative model
#'   \item \code{p_yx <= alpha} and \code{p_xy <= alpha}: confounding
#'   \item \code{p_yx > alpha} and \code{p_xy > alpha}: undecided
#' }
#'
#' Difference statistics use their bootstrap confidence intervals. An
#' interval above zero speaks for the target model, an interval below zero
#' speaks for the alternative model, and an interval containing zero is
#' undecided. As of version 0.2.0 this holds for \code{dda.resdist} under
#' both \code{prob.trans = FALSE} and \code{prob.trans = TRUE}. Difference
#' statistics without a bootstrap interval (\code{B = 0}) are not decided.
#'
#' Run time grows with \code{iter} times the resampling budget (\code{B}) of
#' the original DDA call.
#'
#' @return An object of class \code{dda_bagging} (with subclass
#'   \code{dda_bagging_indep}, \code{dda_bagging_resdist}, or
#'   \code{dda_bagging_vardist}) containing
#'   \item{bagged_results}{The DDA results of the bootstrap samples.}
#'   \item{aggregated_stats}{A DDA object of the original class holding the
#'     aggregated test statistics and p-values.}
#'   \item{decisions}{A matrix of model selection decisions with one row per
#'     bootstrap sample and one column per test.}
#'   \item{decision_proportions}{A matrix of decision proportions with one
#'     row per test.}
#'   \item{ols}{Coefficients, p-values, and R-squared values of the target
#'     and alternative OLS models in each bootstrap sample.}
#'   \item{n_valid_iterations}{Number of bootstrap samples with a DDA
#'     result.}
#'
#' @references
#' Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis:
#' Foundations and Statistical Methods}. Cambridge, UK: Cambridge University
#' Press.
#'
#' @seealso \code{\link{dda.indep}}, \code{\link{dda.vardist}},
#'   \code{\link{dda.resdist}}, \code{\link{summary.dda_bagging}},
#'   \code{\link{summary_ols}}
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' ## --- Fit a base DDA independence model
#'
#' base_model <- dda.indep(y ~ x, pred = "x", data = d, B = 20, hetero = TRUE)
#'
#' ## --- Bootstrap aggregation of the base model
#'
#' bagged <- dda.bagging(base_model, data = d, iter = 5, progress = FALSE)
#' # Note: n, B and iter are kept small here to lower computation time.
#'
#' print(bagged)
#' summary(bagged, show = c("hsic", "dcor"))
#'
#' \dontrun{
#' ## --- Realistic settings; run time is substantial
#'
#' base_model <- dda.indep(y ~ x, pred = "x", data = d, B = 500,
#'   hetero = TRUE, nlfun = 2, diff = TRUE)
#'
#' bagged <- dda.bagging(base_model, data = d, iter = 200,
#'   agg_stat = "trimmed", trim_prob = 0.05)
#'
#' print(bagged)
#' print(bagged, agg_stat = "median")
#' summary(bagged, show = c("hsic", "dcor", "bp"))
#' summary_ols(bagged)
#' }
#'
#' @export
#' @rdname dda.bagging
dda.bagging <- function(dda_result,
                        data,
                        iter = 100,
                        alpha = 0.05,
                        agg_stat = c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey"),
                        trim_prob = 0.10,
                        win_prob = 0.10,
                        # inner_B = NULL,
                        progress = TRUE,
                        save_file = NULL
                        ){

  ### --- check input

  if (!inherits(dda_result, c("dda.indep", "dda.resdist", "dda.vardist")))
    stop("dda_result must be a dda.indep, dda.resdist, or dda.vardist object.")
  if (missing(data) || !is.data.frame(data)) stop("Please provide the data.frame used to fit dda_result.")
  agg_stat <- match.arg(agg_stat)

  ### --- target and alternative model formulas

  y.name <- dda_result$var.names[1]  # tentative outcome
  x.name <- dda_result$var.names[2]  # tentative predictor

  formula.tar <- formula(dda_result$call_info$formula)
  formula.alt <- stats::update(formula.tar, paste(x.name, "~ . -", x.name, "+", y.name))

  ### --- original DDA call, evaluated again in every bootstrap sample

  boot.call <- dda_result$call_info$function_call
  boot.call$formula <- formula.tar
  boot.call$data <- quote(boot.data)
  # if (!is.null(inner_B)) boot.call$B <- inner_B

  # boot.data is stored here, and all other arguments of the original call are
  # found where dda.bagging is init called
  boot.env <- new.env(parent = parent.frame())

  ### --- bootstrap loop

  nobs <- nrow(data)
  results <- vector("list", iter)
  ols.tar.coef <- ols.tar.p <- ols.tar.r2 <- NULL
  ols.alt.coef <- ols.alt.p <- ols.alt.r2 <- NULL

  if (progress) pb <- txtProgressBar(min = 0, max = iter, style = 3)

  for (i in 1:iter) {

    boot.env$boot.data <- data[sample(1:nobs, nobs, replace = TRUE), ]

    fit <- try(eval(boot.call, boot.env), silent = TRUE)
    if (progress) setTxtProgressBar(pb, i)
    if (inherits(fit, "try-error")) next

    fit$call_info <- NULL
    results[[i]] <- fit

    s.tar <- summary(lm(formula.tar, data = boot.env$boot.data))
    s.alt <- summary(lm(formula.alt, data = boot.env$boot.data))

    ols.tar.coef <- rbind(ols.tar.coef, s.tar$coefficients[, 1])
    ols.tar.p    <- rbind(ols.tar.p,    s.tar$coefficients[, 4])
    ols.tar.r2   <- rbind(ols.tar.r2,   c(r.squared = s.tar$r.squared, adj.r.squared = s.tar$adj.r.squared))
    ols.alt.coef <- rbind(ols.alt.coef, s.alt$coefficients[, 1])
    ols.alt.p    <- rbind(ols.alt.p,    s.alt$coefficients[, 4])
    ols.alt.r2   <- rbind(ols.alt.r2,   c(r.squared = s.alt$r.squared, adj.r.squared = s.alt$adj.r.squared))
  }

  if (progress) close(pb)

  results <- results[!sapply(results, is.null)]
  n.valid <- length(results)
  if (n.valid == 0) stop("The DDA call failed in every bootstrap sample.")

  ### --- aggregated test statistics

  agg <- bag.aggregate(results, agg_stat, trim_prob, win_prob)

  ### --- model selection decisions and their proportions

  tests <- names(dda.decisions(results[[1]], alpha))
  decisions <- matrix(NA, nrow = n.valid, ncol = length(tests), dimnames = list(NULL, tests))
  for (i in 1:n.valid) decisions[i, ] <- dda.decisions(results[[i]], alpha)[tests]

  if (inherits(dda_result, "dda.indep")) {
    levs <- c("Target", "Alternative", "Confounding", "Undecided")
  } else {
    levs <- c("Target", "Alternative", "Undecided")
  }

  props <- matrix(NA, nrow = ncol(decisions), ncol = length(levs),
                  dimnames = list(colnames(decisions), levs))
  for (test in colnames(decisions)) {
    dec <- decisions[, test]
    dec <- dec[!is.na(dec)]
    for (lev in levs) props[test, lev] <- mean(dec == lev)
  }

  ### --- output

  output <- list(bagged_results = results,
                 aggregated_stats = agg,
                 decisions = decisions,
                 decision_proportions = props,
                 ols = list(target = list(formula = formula.tar, coef = ols.tar.coef,
                                          p.value = ols.tar.p, r.squared = ols.tar.r2),
                            alternative = list(formula = formula.alt, coef = ols.alt.coef,
                                               p.value = ols.alt.p, r.squared = ols.alt.r2)),
                 n_valid_iterations = n.valid,
                 alpha = alpha,
                 agg_stat = agg_stat,
                 trim_prob = trim_prob,
                 win_prob = win_prob)

  type <- sub("dda.", "", class(dda_result)[1], fixed = TRUE)
  class(output) <- c(paste0("dda_bagging_", type), "dda_bagging")

  if (!is.null(save_file)) saveRDS(output, file = save_file)

  return(output)
}


#' @title Print Method for \code{dda_bagging} Objects
#'
#' @param x An object of class \code{dda_bagging} when using \code{print}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @examples
#' print(bagged)
#'
#' @export
#' @rdname dda.bagging
#' @method print dda_bagging
print.dda_bagging <- function(x, agg_stat = NULL, trim_prob = x$trim_prob, win_prob = x$win_prob, ...){

  object <- reaggregate_bagging(x, agg_stat, trim_prob, win_prob)

  cat("\n")
  cat("BOOTSTRAP AGGREGATED DDA", "\n")
  cat(paste("Number of bootstrap samples:", object$n_valid_iterations), "\n")
  cat(paste("Aggregation method:", object$agg_stat), "\n")
  cat("Test statistics and p-values are aggregated across bootstrap samples.", "\n")

  print(object$aggregated_stats)

  invisible(x)
}


#' @title Re-aggregate a Bootstrap Aggregated DDA Object
#'
#' @description Recomputes the aggregated test statistics of a
#'   \code{dda_bagging} object with a different aggregation method. The
#'   bootstrap results are reused; no resampling is done. With
#'   \code{agg_stat = NULL} the object is returned unchanged.
#'
#' @keywords internal
#' @noRd
reaggregate_bagging <- function(object, agg_stat = NULL, trim_prob = object$trim_prob, win_prob = object$win_prob){

  if (is.null(agg_stat)) return(object)
  agg_stat <- match.arg(agg_stat, c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey"))

  object$aggregated_stats <- bag.aggregate(object$bagged_results, agg_stat, trim_prob, win_prob)
  object$agg_stat  <- agg_stat
  object$trim_prob <- trim_prob
  object$win_prob  <- win_prob
  object
}


#' @title Aggregate DDA Results Across Bootstrap Samples
#'
#' @description Returns a DDA object of the same class as the bootstrap
#'   results in which every test statistic and p-value is the aggregate of
#'   that value across bootstrap samples. The object prints with the print
#'   method of the original DDA function.
#'
#' @keywords internal
#' @noRd
bag.aggregate <- function(results, agg_stat, trim_prob = 0.10, win_prob = 0.10){

  first <- results[[1]]

  # aggregate of one stored value; index is its position in a DDA object,
  # e.g. c("hsic.yx", "p.value") or list("breusch_pagan", 2, "p.value")
  bag <- function(index){
    out <- first
    for (k in index) out <- out[[k]]
    values <- matrix(NA, nrow = length(results), ncol = length(out))
    for (i in 1:length(results)) {
      value <- results[[i]]
      for (k in index) value <- value[[k]]
      if (!is.null(value)) values[i, ] <- as.numeric(value)
    }
    for (j in 1:length(out)) out[j] <- agg.value(values[, j], agg_stat, trim_prob, win_prob)
    out
  }

  agg <- list()

  if (inherits(first, "dda.vardist")) {

    for (test in c("agostino", "anscombe")) {
      for (m in c("predictor", "outcome")) {
        agg[[test]][[m]] <- list(statistic = bag(c(test, m, "statistic")),
                                 p.value = bag(c(test, m, "p.value")))
      }
    }

    for (s in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) agg[[s]] <- bag(s)

    agg$boot.args <- first$boot.args
    agg$boot.warning <- sign(agg$anscombe$predictor$statistic[1]) != sign(agg$anscombe$outcome$statistic[1])
  }

  if (inherits(first, "dda.resdist")) {

    for (test in c("agostino", "anscombe")) {
      for (m in c("target", "alternative")) {
        agg[[test]][[m]] <- list(statistic = bag(c(test, m, "statistic")),
                                 p.value = bag(c(test, m, "p.value")))
      }
    }

    agg$skewdiff <- bag("skewdiff")
    agg$kurtdiff <- bag("kurtdiff")

    if (!is.null(first$cor12diff)) {
      for (s in c("cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) agg[[s]] <- bag(s)
      agg$boot.args <- first$boot.args
      agg$boot.warning <- FALSE
    }

    agg$probtrans <- first$probtrans
  }

  if (inherits(first, "dda.indep")) {

    for (m in c("hsic.yx", "hsic.xy")) {
      agg[[m]] <- list(statistic = bag(c(m, "statistic")),
                       p.value = bag(c(m, "p.value")))
    }
    agg$hsic.method <- first$hsic.method

    for (m in c("distance_cor.dcor_yx", "distance_cor.dcor_xy")) {
      agg[[m]] <- list(statistic = bag(c(m, "statistic")),
                       p.value = bag(c(m, "p.value")))
    }

    if (!is.null(first$breusch_pagan)) {
      agg$breusch_pagan <- list()
      for (k in 1:4) {
        agg$breusch_pagan[[k]] <- list(statistic = bag(list("breusch_pagan", k, "statistic")),
                                       parameter = bag(list("breusch_pagan", k, "parameter")),
                                       p.value = bag(list("breusch_pagan", k, "p.value")))
      }
    }

    if (!is.null(first$nlfun)) {
      for (m in c("nlcor.yx", "nlcor.xy")) {
        agg[[m]] <- list(t1 = bag(c(m, "t1")),
                         t2 = bag(c(m, "t2")),
                         t3 = bag(c(m, "t3")),
                         func = first[[m]]$func)
      }
      agg$nlfun <- first$nlfun
    }

    if (!is.null(first$out.diff)) {
      agg$out.diff <- bag("out.diff")
      agg$boot.args <- first$boot.args
      agg$boot.warning <- FALSE
    }
  }

  agg$var.names <- first$var.names
  class(agg) <- class(first)
  agg
}


#' @title Aggregate a Numeric Vector
#'
#' @description Summarizes the values of one statistic across bootstrap
#'   samples. Non-finite values are dropped.
#'
#' @keywords internal
#' @noRd
agg.value <- function(x, agg_stat = "mean", trim_prob = 0.10, win_prob = 0.10){

  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA)

  if (agg_stat == "mean")     return(mean(x))
  if (agg_stat == "median")   return(median(x))
  if (agg_stat == "trimmed")  return(mean(x, trim = trim_prob))
  if (agg_stat == "midhinge") return(mean(quantile(x, probs = c(0.25, 0.75), names = FALSE)))

  if (agg_stat == "winsorized") {
    q <- quantile(x, probs = c(win_prob, 1 - win_prob), names = FALSE)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(mean(x))
  }

  if (agg_stat == "tukey") {
    q <- quantile(x, probs = c(0.25, 0.50, 0.75), names = FALSE)
    return((q[1] + 2 * q[2] + q[3]) / 4)
  }

  stop("Unknown agg_stat. Choose 'mean', 'median', 'trimmed', 'winsorized', 'midhinge', or 'tukey'.")
}


#' @title Model Selection Decisions for a DDA Object
#'
#' @description Translates the tests stored in a \code{dda.indep},
#'   \code{dda.resdist}, or \code{dda.vardist} object into model selection
#'   decisions (\code{"Target"}, \code{"Alternative"}, \code{"Confounding"},
#'   or \code{"Undecided"}). The rules are listed in the details of
#'   \code{dda.bagging}.
#'
#' @param dda_result An output object from \code{dda.indep},
#'   \code{dda.resdist}, or \code{dda.vardist}.
#' @param alpha Numeric. Significance level (default: 0.05).
#'
#' @return A named character vector with one decision per test.
#'
#' @keywords internal
#' @noRd
dda.decisions <- function(dda_result, alpha = 0.05){

  # separate tests: p.yx from the target model, p.xy from the alternative model
  decide.p <- function(p.yx, p.xy, both.significant = "Undecided"){
    if (is.na(p.yx) || is.na(p.xy)) return(NA)
    if (p.yx >  alpha & p.xy <= alpha) return("Target")
    if (p.yx <= alpha & p.xy >  alpha) return("Alternative")
    if (p.yx <= alpha & p.xy <= alpha) return(both.significant)
    return("Undecided")
  }

  # difference statistics: bootstrap CI with elements named lower and upper
  decide.ci <- function(ci){
    lower <- ci["lower"]
    upper <- ci["upper"]
    if (is.na(lower) || is.na(upper)) return(NA)
    if (lower > 0 & upper > 0) return("Target")
    if (lower < 0 & upper < 0) return("Alternative")
    return("Undecided")
  }

  obj <- dda_result
  dec <- c()

  if (inherits(obj, "dda.vardist")) {

    dec["agostino"] <- decide.p(obj$agostino$outcome$p.value, obj$agostino$predictor$p.value)
    dec["anscombe"] <- decide.p(obj$anscombe$outcome$p.value, obj$anscombe$predictor$p.value)

    for (s in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) dec[s] <- decide.ci(obj[[s]])

  } else if (inherits(obj, "dda.resdist")) {

    # under prob.trans = TRUE the target and alternative p-values swap roles
    if (isFALSE(obj$probtrans)) {
      dec["agostino"] <- decide.p(obj$agostino$target$p.value, obj$agostino$alternative$p.value)
      dec["anscombe"] <- decide.p(obj$anscombe$target$p.value, obj$anscombe$alternative$p.value)
    } else if (isTRUE(obj$probtrans)) {
      dec["agostino"] <- decide.p(obj$agostino$alternative$p.value, obj$agostino$target$p.value)
      dec["anscombe"] <- decide.p(obj$anscombe$alternative$p.value, obj$anscombe$target$p.value)
    } else stop("The prob.trans setting of the dda.resdist object is missing.")

    for (s in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if ("lower" %in% names(obj[[s]])) dec[s] <- decide.ci(obj[[s]])
    }

  } else if (inherits(obj, "dda.indep")) {

    dec["hsic"] <- decide.p(obj$hsic.yx$p.value, obj$hsic.xy$p.value, "Confounding")
    dec["dcor"] <- decide.p(obj$distance_cor.dcor_yx$p.value, obj$distance_cor.dcor_xy$p.value, "Confounding")

    # robust Breusch-Pagan tests are elements 2 (target) and 4 (alternative)
    if (!is.null(obj$breusch_pagan)) {
      dec["bp"] <- decide.p(obj$breusch_pagan[[2]]$p.value, obj$breusch_pagan[[4]]$p.value, "Confounding")
    }

    # smallest p-value of the three non-linear correlation tests
    if (!is.null(obj$nlcor.yx)) {
      p.yx <- min(obj$nlcor.yx$t1[4], obj$nlcor.yx$t2[4], obj$nlcor.yx$t3[4])
      p.xy <- min(obj$nlcor.xy$t1[4], obj$nlcor.xy$t2[4], obj$nlcor.xy$t3[4])
      dec["nlcor"] <- decide.p(p.yx, p.xy, "Confounding")
    }

    if (!is.null(obj$out.diff)) {
      dec["hsic.diff"] <- decide.ci(obj$out.diff["HSIC", ])
      dec["dcor.diff"] <- decide.ci(obj$out.diff["dCor", ])
      dec["mi.diff"]   <- decide.ci(obj$out.diff["MI", ])
    }

  } else stop("dda_result must be a dda.indep, dda.resdist, or dda.vardist object.")

  dec
}
