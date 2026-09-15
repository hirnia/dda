#' @title Model Selection Decisions for a Fitted DDA Object
#'
#' @description
#' \code{dda.decisions} translates the tests stored in a fitted DDA object
#' (\code{dda.indep}, \code{dda.resdist}, or \code{dda.vardist}) into causal
#' model selection decisions. Decisions based on significance tests compare
#' p-values against \code{alpha}; decisions based on bootstrap confidence
#' intervals check whether the interval excludes zero. The same rules are
#' used inside \code{dda.bagging}, so simulation code can call
#' \code{dda.decisions} directly instead of re-implementing decision logic.
#'
#' @param dda_result An output object from \code{dda.indep},
#'   \code{dda.resdist}, or \code{dda.vardist}.
#' @param alpha Numeric. Significance level used for causal model selection
#'   (default: 0.05).
#'
#' @details
#' Throughout, the target model is \code{x -> y} and the alternative model is
#' \code{y -> x}; \code{p_yx} denotes a p-value obtained under the target
#' model and \code{p_xy} a p-value obtained under the alternative model.
#'
#' Separate independence tests (\code{dda.indep}: HSIC, dCor, robust
#' Breusch-Pagan, non-linear correlation):
#' \itemize{
#'   \item \code{p_yx > alpha} and \code{p_xy <= alpha}: \code{"Target"}
#'   \item \code{p_yx <= alpha} and \code{p_xy > alpha}: \code{"Alternative"}
#'   \item both \code{<= alpha}: \code{"Confounding"}
#'   \item both \code{> alpha}: \code{"Undecided"}
#' }
#' Difference statistics (HSIC, dCor, MI) use their bootstrap confidence
#' intervals: an interval above zero speaks for the target model, below zero
#' for the alternative model, otherwise the decision is \code{"Undecided"}.
#'
#' Separate normality tests on residuals (\code{dda.resdist}, D'Agostino and
#' Anscombe-Glynn) depend on the \code{prob.trans} setting of the fitted
#' object. Under \code{prob.trans = FALSE}:
#' \itemize{
#'   \item \code{p_yx > alpha} and \code{p_xy <= alpha}: \code{"Target"}
#'   \item \code{p_yx <= alpha} and \code{p_xy > alpha}: \code{"Alternative"}
#'   \item otherwise: \code{"Undecided"}
#' }
#' Under \code{prob.trans = TRUE} the roles are reversed
#' (\code{p_yx <= alpha} and \code{p_xy > alpha} speaks for the target
#' model). Skewness and kurtosis difference intervals are likewise reversed
#' under \code{prob.trans = TRUE}; co-skewness and co-kurtosis intervals are
#' never reversed.
#'
#' Separate normality tests on observed variables (\code{dda.vardist}):
#' \itemize{
#'   \item outcome p \code{> alpha} and predictor p \code{<= alpha}:
#'     \code{"Target"}
#'   \item outcome p \code{<= alpha} and predictor p \code{> alpha}:
#'     \code{"Alternative"}
#'   \item otherwise: \code{"Undecided"}
#' }
#' All higher-moment difference intervals in \code{dda.vardist} use the
#' interval-excludes-zero rule with values above zero speaking for the
#' target model.
#'
#' @return A named character vector of decisions, one element per test
#'   available in the object. Values are \code{"Target"},
#'   \code{"Alternative"}, \code{"Undecided"}, \code{"Confounding"} (separate
#'   independence tests only), or \code{NA} when a test result is missing.
#'
#' @seealso \code{\link{dda.bagging}}
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rchisq(n, df = 3) - 3
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' fit <- dda.vardist(y ~ x, pred = "x", data = d, B = 50)
#' dda.decisions(fit, alpha = 0.05)
#'
#' @export
dda.decisions <- function(dda_result, alpha = 0.05) {

  if (!inherits(dda_result, c("dda.indep", "dda.resdist", "dda.vardist"))) {
    stop("Unsupported DDA object. Must be dda.indep, dda.resdist, or dda.vardist.")
  }
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)

  # first element as numeric, NA if missing
  p1 <- function(x) {
    x <- suppressWarnings(as.numeric(unlist(x)[1]))
    if (length(x) == 0) NA_real_ else x
  }

  # last two elements of a result vector are its confidence bounds;
  # min_len guards against objects fitted with B = 0 (no bootstrap CI)
  ci_tail <- function(v, min_len = 3) {
    v <- suppressWarnings(as.numeric(v))
    if (length(v) < min_len) return(c(NA_real_, NA_real_))
    c(v[length(v) - 1], v[length(v)])
  }

  # separate-test decision from target and alternative model p-values;
  # confounding = TRUE adds the both-significant category used for
  # independence tests
  dec_sep <- function(p_tar, p_alt, confounding = FALSE) {
    if (is.na(p_tar) || is.na(p_alt)) return(NA_character_)
    if (p_tar >  alpha && p_alt <= alpha) return("Target")
    if (p_tar <= alpha && p_alt >  alpha) return("Alternative")
    if (confounding && p_tar <= alpha && p_alt <= alpha) return("Confounding")
    "Undecided"
  }

  # confidence interval decision; flip = TRUE reverses the direction
  # (used for skewness/kurtosis differences under prob.trans = TRUE)
  dec_ci <- function(lower, upper, flip = FALSE) {
    if (is.na(lower) || is.na(upper)) return(NA_character_)
    if (!flip) {
      if (lower > 0 && upper > 0) return("Target")
      if (lower < 0 && upper < 0) return("Alternative")
    } else {
      if (lower < 0 && upper < 0) return("Target")
      if (lower > 0 && upper > 0) return("Alternative")
    }
    "Undecided"
  }

  obj <- dda_result
  out <- character(0)

  if (inherits(obj, "dda.indep")) {

    if (!is.null(obj$hsic.yx)) {
      out["hsic"] <- dec_sep(p1(obj$hsic.yx$p.value), p1(obj$hsic.xy$p.value),
                             confounding = TRUE)
    }

    dcor_yx <- if (!is.null(obj$distance_cor.dcor_yx)) obj$distance_cor.dcor_yx else obj$dcor.yx
    dcor_xy <- if (!is.null(obj$distance_cor.dcor_xy)) obj$distance_cor.dcor_xy else obj$dcor.xy
    if (!is.null(dcor_yx) && !is.null(dcor_xy)) {
      out["dcor"] <- dec_sep(p1(dcor_yx$p.value), p1(dcor_xy$p.value),
                             confounding = TRUE)
    }

    if (!is.null(obj$breusch_pagan) && length(obj$breusch_pagan) >= 4) {
      # elements 2 and 4 hold the robust Breusch-Pagan tests
      out["dec_bp"] <- dec_sep(p1(obj$breusch_pagan[[2]]$p.value),
                               p1(obj$breusch_pagan[[4]]$p.value),
                               confounding = TRUE)
    }

    if (!is.null(obj$nlcor.yx)) {
      # minimum p-value across the three non-linear transformations
      p_yx <- suppressWarnings(min(c(p1(obj$nlcor.yx$t1[4]),
                                     p1(obj$nlcor.yx$t2[4]),
                                     p1(obj$nlcor.yx$t3[4])), na.rm = TRUE))
      p_xy <- suppressWarnings(min(c(p1(obj$nlcor.xy$t1[4]),
                                     p1(obj$nlcor.xy$t2[4]),
                                     p1(obj$nlcor.xy$t3[4])), na.rm = TRUE))
      if (!is.finite(p_yx)) p_yx <- NA_real_
      if (!is.finite(p_xy)) p_xy <- NA_real_
      out["dec_nl.min"] <- dec_sep(p_yx, p_xy, confounding = TRUE)
    }

    if (!is.null(obj$out.diff)) {
      dm <- as.matrix(obj$out.diff)
      nc <- ncol(dm)
      if (nc >= 2) {
        if (nrow(dm) >= 1) out["diff_hsic"] <- dec_ci(dm[1, nc - 1], dm[1, nc])
        if (nrow(dm) >= 2) out["diff_dcor"] <- dec_ci(dm[2, nc - 1], dm[2, nc])
        if (nrow(dm) >= 3) out["diff_mi"]   <- dec_ci(dm[3, nc - 1], dm[3, nc])
      }
    }
  }

  if (inherits(obj, "dda.resdist")) {

    prob_trans <- obj$probtrans

    p_skew_tar <- p1(obj$agostino$target$p.value)
    p_skew_alt <- p1(obj$agostino$alternative$p.value)
    p_kurt_tar <- p1(obj$anscombe$target$p.value)
    p_kurt_alt <- p1(obj$anscombe$alternative$p.value)

    if (isTRUE(prob_trans)) {
      # non-normal true error: non-normal target-model residuals together
      # with normal-looking alternative-model residuals speak for the target
      out["dec_agost"]  <- dec_sep(p_skew_alt, p_skew_tar)
      out["dec_anscom"] <- dec_sep(p_kurt_alt, p_kurt_tar)
    } else if (isFALSE(prob_trans)) {
      # normal true error, non-normal predictor: normal-looking target-model
      # residuals together with non-normal alternative-model residuals speak
      # for the target
      out["dec_agost"]  <- dec_sep(p_skew_tar, p_skew_alt)
      out["dec_anscom"] <- dec_sep(p_kurt_tar, p_kurt_alt)
    } else {
      stop("'prob.trans' status of the dda.resdist object could not be determined.")
    }

    flip <- isTRUE(prob_trans)
    ci <- ci_tail(obj$skewdiff, min_len = 5)
    out["dec_skewdiff"] <- dec_ci(ci[1], ci[2], flip = flip)
    ci <- ci_tail(obj$kurtdiff, min_len = 5)
    out["dec_kurtdiff"] <- dec_ci(ci[1], ci[2], flip = flip)

    for (k in c("cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if (!is.null(obj[[k]])) {
        ci <- ci_tail(obj[[k]], min_len = 3)
        out[paste0("dec_", k)] <- dec_ci(ci[1], ci[2])
      }
    }
  }

  if (inherits(obj, "dda.vardist")) {

    p_skew_out  <- p1(obj$agostino$outcome$p.value)
    p_skew_pred <- p1(obj$agostino$predictor$p.value)
    p_kurt_out  <- p1(obj$anscombe$outcome$p.value)
    p_kurt_pred <- p1(obj$anscombe$predictor$p.value)

    # a normal-looking outcome together with a non-normal predictor speaks
    # for the target model
    out["dec_agost"]  <- dec_sep(p_skew_out, p_skew_pred)
    out["dec_anscom"] <- dec_sep(p_kurt_out, p_kurt_pred)

    for (k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) {
      if (!is.null(obj[[k]])) {
        ci <- ci_tail(obj[[k]], min_len = 3)
        out[paste0("dec_", k)] <- dec_ci(ci[1], ci[2])
      }
    }
  }

  out
}

#' @title Bootstrap Aggregated Direction Dependence Analysis (DDA)
#'
#' @description
#' \code{dda.bagging} performs bootstrap aggregation (bagging) on an existing
#' Direction Dependence Analysis (DDA) object to test the stability and
#' robustness of direction dependence decisions.
#'
#' @param dda_result An output object from any base DDA function (e.g.,
#'   \code{dda.indep}, \code{dda.vardist}, or \code{dda.resdist}).
#' @param iter Number of bootstrap samples (default: 100).
#' @param progress Logical. Whether to display a progress bar during
#'   resampling (default: \code{TRUE}).
#' @param save_file Character. Optional file path used to save the resulting
#'   object as an R serialized data file (e.g., \code{"results.rds"}).
#' @param alpha Numeric. Significance level used for causal model selection
#'   (default: 0.05).
#' @param data A \code{data.frame} containing ALL raw variables used in the
#'   original DDA model (outcome, predictor, and any covariates).
#' @param agg_stat Character. Specifies the method used for aggregating test
#'   statistics and coefficients across bootstrap samples. Must be one of the
#'   following specifications \code{c("mean", "median", "trimmed",
#'   "winsorized", "midhinge", "tukey")}.
#' @param trim_prob Numeric. Proportion of observations to be trimmed from
#'   each side of the sampling distribution when \code{agg_stat = "trimmed"}
#'   (default: 0.10).
#' @param win_prob Numeric. Proportion of observations to be winsorized from
#'   each side of the sampling distribution when
#'   \code{agg_stat = "winsorized"} (default: 0.10).
#' @param inner_B Optional positive integer. Caps the number of inner
#'   bootstrap resamples (\code{B}) passed to each per-iteration DDA call.
#'   \code{NULL} (default) keeps whatever \code{B} was used in the original
#'   DDA call. Every base DDA function runs its own resampling on each of the
#'   \code{iter} outer iterations, so total cost grows as
#'   \code{iter} \eqn{\times} \code{B}; capping \code{inner_B} is the main
#'   lever for reducing run time without changing how the outer bagging loop
#'   aggregates results. The effect is largest for \code{dda.indep} called
#'   with \code{diff = TRUE}, where each outer iteration runs a full
#'   difference-statistic bootstrap, and for \code{dda.resdist} and
#'   \code{dda.vardist}, which bootstrap their confidence intervals
#'   internally.
#'
#' @details
#' This function uses a fitted DDA output object (obtained from
#' \code{dda.indep}, \code{dda.vardist}, or \code{dda.resdist}) and performs
#' bootstrap aggregation of DDA test statistics. The function computes DDA
#' statistics across \code{iter} bootstrap samples and aggregates the results
#' to evaluate the stability and robustness of DDA model selection. p-values
#' obtained from significance tests are aggregated using the harmonic mean
#' p-value approach (Wilson, 2019). Model selection decisions within each
#' bootstrap sample are obtained with \code{\link{dda.decisions}} (p-value
#' based rules for significance tests, interval-excludes-zero rules for
#' bootstrap difference statistics) and reported as decision proportions
#' across samples.
#'
#' Run time scales with \code{iter} multiplied by the resampling budget of the
#' base DDA call, and the underlying independence statistics are quadratic in
#' the number of observations. Use \code{inner_B} to cap the inner budget.
#' Note also that if the original DDA call used \code{parallelize = TRUE},
#' that setting is inherited by every outer iteration, so a new cluster is
#' started \code{iter} times; for the small inner bootstraps typical of
#' bagging this is usually slower than running the inner calls serially.
#'
#' @references
#' Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis:
#' Foundations and Statistical Methods}. Cambridge, UK: Cambridge University
#' Press.
#'
#' Wilson, D. J. (2019). The harmonic mean p-value for combining dependent
#' tests. \emph{Proceedings of the National Academy of Sciences}, \emph{116}(4),
#' 1195--1200.
#'
#' @return An object of class \code{dda_bagging} (with subclasses
#'   \code{dda_bagging_indep}, \code{dda_bagging_vardist}, or
#'   \code{dda_bagging_resdist}), which contains aggregated and raw results
#'   from tests matching the initial DDA function (\code{dda.indep},
#'   \code{dda.vardist}, or \code{dda.resdist}).
#'
#' @seealso \code{\link{dda.indep}}, \code{\link{dda.vardist}},
#'   \code{\link{dda.resdist}}, \code{\link{dda.decisions}}
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rchisq(n, df = 3) - 3
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' ## --- Fit a base DDA independence model
#'
#' base_model <- dda.indep(y ~ x, pred = "x", data = d, B = 10, hetero = TRUE)
#'
#' ## --- Bootstrap aggregation of the base model
#'
#' bagged <- dda.bagging(base_model, data = d, iter = 5, agg_stat = "mean",
#'   inner_B = 10, progress = FALSE)
#' # Note: n, B and iter are all kept small here to lower computation time.
#' # inner_B caps the resampling budget of each outer iteration.
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
#'   agg_stat = "trimmed", trim_prob = 0.05, progress = TRUE)
#'
#' print(bagged)
#' summary(bagged, show = c("hsic", "dcor", "bp"))
#' print_ols_summary(bagged)
#' }
#'
#' @export
dda.bagging <- function(
    dda_result,
    iter         = 100,
    progress     = TRUE,
    save_file    = NULL,
    alpha        = 0.05,
    data         = NULL,
    agg_stat     = c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey"),
    trim_prob    = 0.10,
    win_prob     = 0.10,
    inner_B      = NULL
) {

  # Capture caller environment immediately. Symbols stored in call_info$all_args
  # (e.g. data = dat, B = my_b) are resolved here, not inside the loop.
  caller_env <- parent.frame()

  agg_stat <- match.arg(agg_stat)

  # --- Input Validation ---
  if (!inherits(dda_result, c("dda.indep", "dda.resdist", "dda.vardist"))) {
    stop("Unsupported DDA object. Must be dda.indep, dda.resdist, or dda.vardist.")
  }
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("Please provide a valid 'data' data.frame.")
  }
  stopifnot(
    is.numeric(iter)     && iter > 0,
    is.numeric(alpha)    && alpha > 0 && alpha < 1,
    is.numeric(trim_prob) && trim_prob >= 0 && trim_prob < 0.5,
    is.numeric(win_prob)  && win_prob  >= 0 && win_prob  < 0.5
  )
  if (!is.null(inner_B)) {
    stopifnot(is.numeric(inner_B) && length(inner_B) == 1 && inner_B > 0)
    inner_B <- as.integer(inner_B)
  }

  # --- Helper: Robust & Finite Aggregation ---
  agg_helper <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x) & !is.nan(x) & is.finite(x)]
    if (length(x) == 0) return(NA_real_)

    switch(agg_stat,
           "mean"       = mean(x),
           "median"     = median(x),
           "trimmed"    = mean(x, trim = trim_prob),
           "winsorized" = {
             q_low  <- quantile(x, probs = win_prob,     na.rm = TRUE, names = FALSE)
             q_high <- quantile(x, probs = 1 - win_prob, na.rm = TRUE, names = FALSE)
             x[x < q_low]  <- q_low
             x[x > q_high] <- q_high
             mean(x)
           },
           "midhinge" = {
             q <- quantile(x, probs = c(0.25, 0.75), names = FALSE, na.rm = TRUE)
             mean(q)
           },
           "tukey" = {
             q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
             (q[1] + 2*q[2] + q[3]) / 4
           }
    )
  }

  # --- Helper: Safe Numeric Extraction ---
  get_numeric <- function(x) {
    if (is.null(x))    return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x[1]))
    if (is.list(x))    return(get_numeric(x[[1]]))
    return(NA_real_)
  }

  # --- Helper: Harmonic Mean P-values ---
  harmonic_p <- function(pvec) {
    pvec <- as.numeric(pvec)
    pvec <- pvec[!is.na(pvec) & !is.nan(pvec)]
    if (length(pvec) == 0) return(NA_real_)
    pvec[pvec <= 0] <- 1e-300

    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
      warning("Package 'harmonicmeanp' not found. Falling back to arithmetic mean.")
      return(mean(pvec))
    }
    harmonicmeanp::p.hmp(pvec, L = length(pvec))
  }

  # --- Helper: Matrix Aggregation that Keeps P-Value Columns Harmonic ---
  # 5-column results (resdist skewness/kurtosis differences) are
  # c(diff, z, p, lower, upper); column 3 is a p-value and is combined with
  # the harmonic mean instead of agg_stat
  agg_mat <- function(mat) {
    if (is.null(mat)) return(NULL)
    res <- apply(mat, 2, agg_helper)
    if (ncol(mat) == 5) res[3] <- harmonic_p(mat[, 3])
    res
  }

  # --- Helper: Decision Proportions ---
  calc_props <- function(dec_vec, levs = c("Target", "Alternative", "Undecided")) {
    dec_vec <- dec_vec[!is.na(dec_vec)]
    tab     <- table(factor(dec_vec, levels = levs))
    sm      <- sum(tab)
    if (sm == 0) {
      empty <- rep(0, length(levs))
      names(empty) <- levs
      return(empty)
    }
    return(tab / sm)
  }

  # --- Extract Core DDA Information ---
  call_info     <- dda_result$call_info
  original_data <- data
  nobs          <- nrow(original_data)
  dda_func      <- get(call_info$function_name)
  obj_type      <- class(dda_result)[1]

  var_names <- dda_result$var.names
  if (is.null(var_names) || length(var_names) != 2) {
    stop("Variable names not found in DDA result.")
  }
  y_name <- var_names[1]
  x_name <- var_names[2]

  # --- Formula Extraction ---
  original_formula <- NULL
  if (!is.null(call_info$formula)) {
    if (inherits(call_info$formula, "formula")) {
      original_formula <- call_info$formula
    } else if (inherits(call_info$formula, "lm")) {
      original_formula <- formula(call_info$formula)
    } else if (!is.null(call_info$all_args$formula)) {
      if (inherits(call_info$all_args$formula, "formula")) {
        original_formula <- call_info$all_args$formula
      } else if (inherits(call_info$all_args$formula, "lm")) {
        original_formula <- formula(call_info$all_args$formula)
      } else {
        original_formula <- tryCatch(as.formula(call_info$all_args$formula), error = function(e) NULL)
      }
    }
  }
  if (is.null(original_formula)) stop("Could not extract formula from DDA result.")

  all_vars       <- all.vars(original_formula)
  cov_names      <- setdiff(all_vars, c(y_name, x_name))
  has_covariates <- length(cov_names) > 0

  if (has_covariates) {
    cov_str          <- paste(cov_names, collapse = " + ")
    full_formula_tar <- as.formula(paste(y_name, "~", x_name, "+", cov_str))
    full_formula_alt <- as.formula(paste(x_name, "~", y_name, "+", cov_str))
  } else {
    full_formula_tar <- as.formula(paste(y_name, "~", x_name))
    full_formula_alt <- as.formula(paste(x_name, "~", y_name))
  }

  # Pre-evaluate all call_info$all_args once in the caller's frame.
  #
  # match.call() stores arguments as unevaluated language objects, so
  # passing them raw via do.call() inside dda.bagging resolves symbols in
  # the wrong environment and silently returns NULL for every iteration.
  # eval(..., envir = caller_env) turns each element into its actual R value.
  # tryCatch keeps the original object if evaluation fails (e.g. a literal
  # formula or a function object that doesn't need eval).
  boot_args_pre <- lapply(call_info$all_args, function(a) {
    tryCatch(eval(a, envir = caller_env), error = function(e) a)
  })
  # Guard: drop first element if it has no name (positional formula arg).
  if (length(boot_args_pre) > 0 && names(boot_args_pre)[1] == "")
    boot_args_pre[[1]] <- NULL

  # Formulas are constant across iterations; build them once. as.formula()
  # parses and evaluates, so rebuilding them inside the loop costs 3 parses
  # per iteration for no benefit.
  inner_formula <- as.formula(paste(y_name, "~", x_name))
  if (has_covariates) {
    cov_formula_y <- as.formula(paste(y_name, "~", paste(cov_names, collapse = " + ")))
    cov_formula_x <- as.formula(paste(x_name, "~", paste(cov_names, collapse = " + ")))
  }

  # --- Bootstrap Execution ---
  bagged_results <- vector("list", iter)
  ols_tar_coefs  <- vector("list", iter)
  ols_alt_coefs  <- vector("list", iter)
  ols_tar_rsq    <- vector("list", iter)
  ols_alt_rsq    <- vector("list", iter)
  ols_tar_pvals  <- vector("list", iter)
  ols_alt_pvals  <- vector("list", iter)

  if (progress) pb <- txtProgressBar(min = 0, max = iter, style = 3)

  for (i in 1:iter) {

    # Draw n observations with replacement from the original data.
    boot_indices <- sample(1:nobs, nobs, replace = TRUE)
    datboot      <- original_data[boot_indices, ]

    # Fit OLS target and alternative models on the bootstrap sample.
    lm_tar <- tryCatch(lm(full_formula_tar, data = datboot), error = function(e) NULL)
    lm_alt <- tryCatch(lm(full_formula_alt, data = datboot), error = function(e) NULL)

    if (!is.null(lm_tar)) {
      s <- summary(lm_tar)
      ols_tar_coefs[[i]] <- coef(lm_tar)
      ols_tar_rsq[[i]]   <- c(s$r.squared, s$adj.r.squared)
      ols_tar_pvals[[i]] <- s$coefficients[, 4]
    }
    if (!is.null(lm_alt)) {
      s <- summary(lm_alt)
      ols_alt_coefs[[i]] <- coef(lm_alt)
      ols_alt_rsq[[i]]   <- c(s$r.squared, s$adj.r.squared)
      ols_alt_pvals[[i]] <- s$coefficients[, 4]
    }

    # Residualize covariates if present, then scale the working variables.
    if (has_covariates) {
      tryCatch({
        ry <- as.vector(scale(resid(lm(cov_formula_y, data = datboot))))
        rx <- as.vector(scale(resid(lm(cov_formula_x, data = datboot))))
      }, error = function(e) stop(paste("Covariate residualization failed:", e$message)))
    } else {
      ry <- as.vector(scale(datboot[[y_name]]))
      rx <- as.vector(scale(datboot[[x_name]]))
    }

    # Build the argument list for this iteration's DDA call, starting from
    # the pre-evaluated original arguments.
    boot_args         <- boot_args_pre
    boot_args$formula <- inner_formula
    boot_args$pred    <- x_name

    boot_processed        <- data.frame(rx, ry)
    names(boot_processed) <- c(x_name, y_name)
    boot_args$data        <- boot_processed

    # If inner_B is specified, override the inner bootstrap B passed to the
    # DDA function on this iteration. This caps the number of resamples used
    # inside dda.resdist / dda.vardist on each outer iteration, which can
    # dramatically reduce total run time without changing how the outer
    # bagging loop aggregates across iter samples.
    if (!is.null(inner_B)) boot_args$B <- inner_B

    bagged_results[[i]] <- tryCatch(
      do.call(dda_func, boot_args),
      error = function(e) NULL
    )

    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)

  # --- Filter Valid Results ---
  is_valid  <- vapply(bagged_results, function(x) !is.null(x) && !all(is.na(x)), logical(1))
  valid_res <- bagged_results[is_valid]
  n_valid   <- length(valid_res)

  if (n_valid == 0) stop("No valid bootstrap iterations succeeded. Check data variance.")

  raw_stats <- list()
  agg       <- list()
  decs      <- list()

  # --- Aggregate OLS Results ---
  raw_stats$ols_tar_coefs <- tryCatch(do.call(rbind, ols_tar_coefs[is_valid]), error = function(e) NULL)
  raw_stats$ols_alt_coefs <- tryCatch(do.call(rbind, ols_alt_coefs[is_valid]), error = function(e) NULL)
  raw_stats$ols_tar_rsq   <- tryCatch(do.call(rbind, ols_tar_rsq[is_valid]),   error = function(e) NULL)
  raw_stats$ols_alt_rsq   <- tryCatch(do.call(rbind, ols_alt_rsq[is_valid]),   error = function(e) NULL)
  raw_stats$ols_tar_pvals <- tryCatch(do.call(rbind, ols_tar_pvals[is_valid]), error = function(e) NULL)
  raw_stats$ols_alt_pvals <- tryCatch(do.call(rbind, ols_alt_pvals[is_valid]), error = function(e) NULL)

  lb <- alpha / 2
  ub <- 1 - alpha / 2

  if (!is.null(raw_stats$ols_tar_coefs) && is.matrix(raw_stats$ols_tar_coefs)) {
    prop_sig_tar   <- apply(raw_stats$ols_tar_pvals, 2, function(x) mean(x < alpha, na.rm = TRUE))
    agg$ols_target <- cbind(
      estimate = apply(raw_stats$ols_tar_coefs, 2, agg_helper),
      apply(raw_stats$ols_tar_coefs, 2, quantile, probs = lb, na.rm = TRUE),
      apply(raw_stats$ols_tar_coefs, 2, quantile, probs = ub, na.rm = TRUE),
      prop_sig_tar
    )
    colnames(agg$ols_target) <- c("estimate", paste0(lb*100, " %"), paste0(ub*100, " %"), paste0("Prop (p<", alpha, ")"))
  }

  if (!is.null(raw_stats$ols_alt_coefs) && is.matrix(raw_stats$ols_alt_coefs)) {
    prop_sig_alt        <- apply(raw_stats$ols_alt_pvals, 2, function(x) mean(x < alpha, na.rm = TRUE))
    agg$ols_alternative <- cbind(
      estimate = apply(raw_stats$ols_alt_coefs, 2, agg_helper),
      apply(raw_stats$ols_alt_coefs, 2, quantile, probs = lb, na.rm = TRUE),
      apply(raw_stats$ols_alt_coefs, 2, quantile, probs = ub, na.rm = TRUE),
      prop_sig_alt
    )
    colnames(agg$ols_alternative) <- c("estimate", paste0(lb*100, " %"), paste0(ub*100, " %"), paste0("Prop(p<", alpha, ")"))
  }

  # ============================================================================
  # indep block
  # ============================================================================
  if (obj_type == "dda.indep") {
    agg$var.names <- var_names

    raw_stats$hsic_yx_stat <- sapply(valid_res, function(x) get_numeric(x$hsic.yx$statistic))
    raw_stats$hsic_xy_stat <- sapply(valid_res, function(x) get_numeric(x$hsic.xy$statistic))
    raw_stats$hsic_yx_pval <- sapply(valid_res, function(x) get_numeric(x$hsic.yx$p.value))
    raw_stats$hsic_xy_pval <- sapply(valid_res, function(x) get_numeric(x$hsic.xy$p.value))

    agg$hsic_yx_stat <- agg_helper(raw_stats$hsic_yx_stat)
    agg$hsic_xy_stat <- agg_helper(raw_stats$hsic_xy_stat)
    agg$hsic_yx_pval <- harmonic_p(raw_stats$hsic_yx_pval)
    agg$hsic_xy_pval <- harmonic_p(raw_stats$hsic_xy_pval)

    if (!is.null(valid_res[[1]]$distance_cor.dcor_yx) || !is.null(valid_res[[1]]$dcor.yx)) {
      dcor_name_yx <- if (!is.null(valid_res[[1]]$distance_cor.dcor_yx)) "distance_cor.dcor_yx" else "dcor.yx"
      dcor_name_xy <- if (!is.null(valid_res[[1]]$distance_cor.dcor_xy)) "distance_cor.dcor_xy" else "dcor.xy"

      raw_stats$dcor_yx_stat <- sapply(valid_res, function(x) get_numeric(x[[dcor_name_yx]]$statistic))
      raw_stats$dcor_xy_stat <- sapply(valid_res, function(x) get_numeric(x[[dcor_name_xy]]$statistic))
      raw_stats$dcor_yx_pval <- sapply(valid_res, function(x) get_numeric(x[[dcor_name_yx]]$p.value))
      raw_stats$dcor_xy_pval <- sapply(valid_res, function(x) get_numeric(x[[dcor_name_xy]]$p.value))

      agg$dcor_yx_stat <- agg_helper(raw_stats$dcor_yx_stat)
      agg$dcor_xy_stat <- agg_helper(raw_stats$dcor_xy_stat)
      agg$dcor_yx_pval <- harmonic_p(raw_stats$dcor_yx_pval)
      agg$dcor_xy_pval <- harmonic_p(raw_stats$dcor_xy_pval)
    }

    if (!is.null(valid_res[[1]]$breusch_pagan)) {
      raw_stats$bp_yx_stat  <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$statistic))
      raw_stats$bp_yx_df    <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$parameter))
      raw_stats$bp_yx_p     <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[1]]$p.value))
      raw_stats$rbp_yx_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$statistic))
      raw_stats$rbp_yx_df   <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$parameter))
      raw_stats$rbp_yx_p    <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[2]]$p.value))
      raw_stats$bp_xy_stat  <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$statistic))
      raw_stats$bp_xy_df    <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$parameter))
      raw_stats$bp_xy_p     <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[3]]$p.value))
      raw_stats$rbp_xy_stat <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$statistic))
      raw_stats$rbp_xy_df   <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$parameter))
      raw_stats$rbp_xy_p    <- sapply(valid_res, function(x) get_numeric(x$breusch_pagan[[4]]$p.value))

      agg$breusch_pagan <- list(
        list(statistic = agg_helper(raw_stats$bp_yx_stat),  parameter = agg_helper(raw_stats$bp_yx_df),  p.value = harmonic_p(raw_stats$bp_yx_p)),
        list(statistic = agg_helper(raw_stats$rbp_yx_stat), parameter = agg_helper(raw_stats$rbp_yx_df), p.value = harmonic_p(raw_stats$rbp_yx_p)),
        list(statistic = agg_helper(raw_stats$bp_xy_stat),  parameter = agg_helper(raw_stats$bp_xy_df),  p.value = harmonic_p(raw_stats$bp_xy_p)),
        list(statistic = agg_helper(raw_stats$rbp_xy_stat), parameter = agg_helper(raw_stats$rbp_xy_df), p.value = harmonic_p(raw_stats$rbp_xy_p))
      )
    }

    if (!is.null(valid_res[[1]]$nlcor.yx)) {
      raw_stats$nlcor_yx_t1 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t1)))
      raw_stats$nlcor_yx_t2 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t2)))
      raw_stats$nlcor_yx_t3 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.yx$t3)))
      raw_stats$nlcor_xy_t1 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t1)))
      raw_stats$nlcor_xy_t2 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t2)))
      raw_stats$nlcor_xy_t3 <- do.call(rbind, lapply(valid_res, function(x) as.numeric(x$nlcor.xy$t3)))

      agg$nlcor.yx <- list(
        t1   = c(agg_helper(raw_stats$nlcor_yx_t1[,1]), agg_helper(raw_stats$nlcor_yx_t1[,2]), agg_helper(raw_stats$nlcor_yx_t1[,3]), harmonic_p(raw_stats$nlcor_yx_t1[,4])),
        t2   = c(agg_helper(raw_stats$nlcor_yx_t2[,1]), agg_helper(raw_stats$nlcor_yx_t2[,2]), agg_helper(raw_stats$nlcor_yx_t2[,3]), harmonic_p(raw_stats$nlcor_yx_t2[,4])),
        t3   = c(agg_helper(raw_stats$nlcor_yx_t3[,1]), agg_helper(raw_stats$nlcor_yx_t3[,2]), agg_helper(raw_stats$nlcor_yx_t3[,3]), harmonic_p(raw_stats$nlcor_yx_t3[,4])),
        func = valid_res[[1]]$nlcor.yx$func
      )
      agg$nlcor.xy <- list(
        t1   = c(agg_helper(raw_stats$nlcor_xy_t1[,1]), agg_helper(raw_stats$nlcor_xy_t1[,2]), agg_helper(raw_stats$nlcor_xy_t1[,3]), harmonic_p(raw_stats$nlcor_xy_t1[,4])),
        t2   = c(agg_helper(raw_stats$nlcor_xy_t2[,1]), agg_helper(raw_stats$nlcor_xy_t2[,2]), agg_helper(raw_stats$nlcor_xy_t2[,3]), harmonic_p(raw_stats$nlcor_xy_t2[,4])),
        t3   = c(agg_helper(raw_stats$nlcor_xy_t3[,1]), agg_helper(raw_stats$nlcor_xy_t3[,2]), agg_helper(raw_stats$nlcor_xy_t3[,3]), harmonic_p(raw_stats$nlcor_xy_t3[,4])),
        func = valid_res[[1]]$nlcor.xy$func
      )

    }

    if (!is.null(valid_res[[1]]$out.diff)) {
      raw_stats$diff_arr <- simplify2array(lapply(valid_res, function(x) as.matrix(x$out.diff)))
      agg$diff_matrix    <- apply(raw_stats$diff_arr, c(1, 2), agg_helper)
    }
  }

  # ============================================================================
  # resdist block
  # ============================================================================
  if (obj_type == "dda.resdist") {
    agg$var.names <- var_names

    # Store the prob.trans flag so print methods can retrieve it for footnotes.
    prob_trans_flag <- isTRUE(if (!is.null(dda_result$probtrans)) dda_result$probtrans else FALSE)
    agg$probtrans   <- prob_trans_flag

    raw_stats$agost_tar_stat  <- sapply(valid_res, function(x) get_numeric(x$agostino$target$statistic[1]))
    raw_stats$agost_tar_z     <- sapply(valid_res, function(x) get_numeric(x$agostino$target$statistic[2]))
    raw_stats$agost_tar_pval  <- sapply(valid_res, function(x) get_numeric(x$agostino$target$p.value))
    raw_stats$agost_alt_stat  <- sapply(valid_res, function(x) get_numeric(x$agostino$alternative$statistic[1]))
    raw_stats$agost_alt_z     <- sapply(valid_res, function(x) get_numeric(x$agostino$alternative$statistic[2]))
    raw_stats$agost_alt_pval  <- sapply(valid_res, function(x) get_numeric(x$agostino$alternative$p.value))

    raw_stats$anscom_tar_stat <- sapply(valid_res, function(x) get_numeric(x$anscombe$target$statistic[1]))
    raw_stats$anscom_tar_z    <- sapply(valid_res, function(x) get_numeric(x$anscombe$target$statistic[2]))
    raw_stats$anscom_tar_pval <- sapply(valid_res, function(x) get_numeric(x$anscombe$target$p.value))
    raw_stats$anscom_alt_stat <- sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$statistic[1]))
    raw_stats$anscom_alt_z    <- sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$statistic[2]))
    raw_stats$anscom_alt_pval <- sapply(valid_res, function(x) get_numeric(x$anscombe$alternative$p.value))

    agg$agostino.target.statistic      <- agg_helper(raw_stats$agost_tar_stat)
    agg$agostino.target.z              <- agg_helper(raw_stats$agost_tar_z)
    agg$agostino.target.p.value        <- harmonic_p(raw_stats$agost_tar_pval)
    agg$agostino.alternative.statistic <- agg_helper(raw_stats$agost_alt_stat)
    agg$agostino.alternative.z         <- agg_helper(raw_stats$agost_alt_z)
    agg$agostino.alternative.p.value   <- harmonic_p(raw_stats$agost_alt_pval)

    agg$anscombe.target.statistic      <- agg_helper(raw_stats$anscom_tar_stat)
    agg$anscombe.target.z              <- agg_helper(raw_stats$anscom_tar_z)
    agg$anscombe.target.p.value        <- harmonic_p(raw_stats$anscom_tar_pval)
    agg$anscombe.alternative.statistic <- agg_helper(raw_stats$anscom_alt_stat)
    agg$anscombe.alternative.z         <- agg_helper(raw_stats$anscom_alt_z)
    agg$anscombe.alternative.p.value   <- harmonic_p(raw_stats$anscom_alt_pval)

    for (k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if (!is.null(valid_res[[1]][[k]])) {
        mat_k <- tryCatch(do.call(rbind, lapply(valid_res, function(x) as.numeric(x[[k]]))), error = function(e) NULL)
        if (!is.null(mat_k) && ncol(mat_k) >= 2) {
          raw_stats[[k]] <- mat_k
          agg[[k]]       <- agg_mat(mat_k)
        }
      }
    }
  }

  # ============================================================================
  # vardist block
  # ============================================================================
  if (obj_type == "dda.vardist") {
    agg$var.names <- var_names

    raw_stats$agost_pre_stat  <- sapply(valid_res, function(x) get_numeric(x$agostino$predictor$statistic[1]))
    raw_stats$agost_pre_z     <- sapply(valid_res, function(x) get_numeric(x$agostino$predictor$statistic[2]))
    raw_stats$agost_pre_pval  <- sapply(valid_res, function(x) get_numeric(x$agostino$predictor$p.value))
    raw_stats$agost_out_stat  <- sapply(valid_res, function(x) get_numeric(x$agostino$outcome$statistic[1]))
    raw_stats$agost_out_z     <- sapply(valid_res, function(x) get_numeric(x$agostino$outcome$statistic[2]))
    raw_stats$agost_out_pval  <- sapply(valid_res, function(x) get_numeric(x$agostino$outcome$p.value))

    raw_stats$anscom_pre_stat <- sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$statistic[1]))
    raw_stats$anscom_pre_z    <- sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$statistic[2]))
    raw_stats$anscom_pre_pval <- sapply(valid_res, function(x) get_numeric(x$anscombe$predictor$p.value))
    raw_stats$anscom_out_stat <- sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$statistic[1]))
    raw_stats$anscom_out_z    <- sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$statistic[2]))
    raw_stats$anscom_out_pval <- sapply(valid_res, function(x) get_numeric(x$anscombe$outcome$p.value))

    agg$agostino.predictor.statistic.skew <- agg_helper(raw_stats$agost_pre_stat)
    agg$agostino.predictor.statistic.z    <- agg_helper(raw_stats$agost_pre_z)
    agg$agostino.predictor.p.value        <- harmonic_p(raw_stats$agost_pre_pval)
    agg$agostino.outcome.statistic.skew   <- agg_helper(raw_stats$agost_out_stat)
    agg$agostino.outcome.statistic.z      <- agg_helper(raw_stats$agost_out_z)
    agg$agostino.outcome.p.value          <- harmonic_p(raw_stats$agost_out_pval)

    agg$anscombe.predictor.statistic.kurt <- agg_helper(raw_stats$anscom_pre_stat)
    agg$anscombe.predictor.statistic.z    <- agg_helper(raw_stats$anscom_pre_z)
    agg$anscombe.predictor.p.value        <- harmonic_p(raw_stats$anscom_pre_pval)
    agg$anscombe.outcome.statistic.kurt   <- agg_helper(raw_stats$anscom_out_stat)
    agg$anscombe.outcome.statistic.z      <- agg_helper(raw_stats$anscom_out_z)
    agg$anscombe.outcome.p.value          <- harmonic_p(raw_stats$anscom_out_pval)

    for (k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) {
      if (!is.null(valid_res[[1]][[k]])) {
        mat_k <- tryCatch(do.call(rbind, lapply(valid_res, function(x) as.numeric(x[[k]]))), error = function(e) NULL)
        if (!is.null(mat_k) && ncol(mat_k) >= 2) {
          raw_stats[[k]] <- mat_k
          agg[[k]]       <- agg_mat(mat_k)
        }
      }
    }
  }

  # ============================================================================
  # model selection decisions
  # ============================================================================
  # One decision per test per bootstrap sample, computed by dda.decisions()
  # so that dda.bagging and simulation code share the same rules; decs holds
  # the proportion of each decision across the valid samples.
  dec_list <- lapply(valid_res, function(x) {
    tryCatch(dda.decisions(x, alpha = alpha), error = function(e) NULL)
  })
  dec_keys <- unique(unlist(lapply(dec_list, names)))
  dec_levs <- if (obj_type == "dda.indep") {
    c("Target", "Alternative", "Confounding", "Undecided")
  } else {
    c("Target", "Alternative", "Undecided")
  }
  for (k in dec_keys) {
    dvec <- vapply(dec_list, function(d) {
      if (!is.null(d) && k %in% names(d)) d[[k]] else NA_character_
    }, character(1))
    decs[[k]] <- calc_props(dvec, levs = dec_levs)
  }

  # --- Compile & Return ---
  if (!is.null(save_file)) {
    saveRDS(
      list(bagged_results       = bagged_results,
           raw_stats            = raw_stats,
           aggregated_stats     = agg,
           decision_percentages = decs,
           n_valid_iterations   = n_valid,
           agg_stat_used        = agg_stat),
      file = save_file
    )
  }

  out <- list(
    bagged_results       = bagged_results,
    raw_stats            = raw_stats,
    aggregated_stats     = agg,
    decision_percentages = decs,
    n_valid_iterations   = n_valid,
    agg_stat_used        = agg_stat
  )
  class(out) <- c(paste0("dda_bagging_", gsub("dda.", "", obj_type)), "dda_bagging")
  return(out)
}
