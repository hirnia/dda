# Internal helper to dynamically re-aggregate statistics from raw vectors
#' @noRd
reaggregate_bagging <- function(object, agg_stat = NULL, trim_prob = 0.10, win_prob = 0.10) {
  if (is.null(agg_stat)) return(object)

  agg_stat <- match.arg(agg_stat, c("mean", "median", "trimmed", "winsorized", "midhinge", "tukey"))
  if (agg_stat == object$agg_stat_used) return(object)

  agg_helper <- function(x) dda_agg(x, agg_stat, trim_prob, win_prob)

  raw <- object$raw_stats
  agg <- object$aggregated_stats
  obj_type <- class(object)[1]

  # Re-aggregate OLS
  if (!is.null(raw$ols_tar_coefs)) agg$ols_target[,"estimate"] <- apply(raw$ols_tar_coefs, 2, agg_helper)
  if (!is.null(raw$ols_alt_coefs)) agg$ols_alternative[,"estimate"] <- apply(raw$ols_alt_coefs, 2, agg_helper)

  # Re-aggregate specific DDA statistics (p-values are not re-aggregated here)
  if (obj_type == "dda_bagging_indep") {
    agg$hsic_yx_stat <- agg_helper(raw$hsic_yx_stat)
    agg$hsic_xy_stat <- agg_helper(raw$hsic_xy_stat)

    if (!is.null(raw$dcor_yx_stat)) {
      agg$dcor_yx_stat <- agg_helper(raw$dcor_yx_stat)
      agg$dcor_xy_stat <- agg_helper(raw$dcor_xy_stat)
    }
    if (!is.null(raw$bp_yx_stat)) {
      agg$breusch_pagan[[1]]$statistic <- agg_helper(raw$bp_yx_stat)
      agg$breusch_pagan[[1]]$parameter <- agg_helper(raw$bp_yx_df)
      agg$breusch_pagan[[2]]$statistic <- agg_helper(raw$rbp_yx_stat)
      agg$breusch_pagan[[2]]$parameter <- agg_helper(raw$rbp_yx_df)
      agg$breusch_pagan[[3]]$statistic <- agg_helper(raw$bp_xy_stat)
      agg$breusch_pagan[[3]]$parameter <- agg_helper(raw$bp_xy_df)
      agg$breusch_pagan[[4]]$statistic <- agg_helper(raw$rbp_xy_stat)
      agg$breusch_pagan[[4]]$parameter <- agg_helper(raw$rbp_xy_df)
    }
    if (!is.null(raw$nlcor_yx_t1)) {
      agg$nlcor.yx$t1[1:3] <- c(agg_helper(raw$nlcor_yx_t1[,1]), agg_helper(raw$nlcor_yx_t1[,2]), agg_helper(raw$nlcor_yx_t1[,3]))
      agg$nlcor.yx$t2[1:3] <- c(agg_helper(raw$nlcor_yx_t2[,1]), agg_helper(raw$nlcor_yx_t2[,2]), agg_helper(raw$nlcor_yx_t2[,3]))
      agg$nlcor.yx$t3[1:3] <- c(agg_helper(raw$nlcor_yx_t3[,1]), agg_helper(raw$nlcor_yx_t3[,2]), agg_helper(raw$nlcor_yx_t3[,3]))

      agg$nlcor.xy$t1[1:3] <- c(agg_helper(raw$nlcor_xy_t1[,1]), agg_helper(raw$nlcor_xy_t1[,2]), agg_helper(raw$nlcor_xy_t1[,3]))
      agg$nlcor.xy$t2[1:3] <- c(agg_helper(raw$nlcor_xy_t2[,1]), agg_helper(raw$nlcor_xy_t2[,2]), agg_helper(raw$nlcor_xy_t2[,3]))
      agg$nlcor.xy$t3[1:3] <- c(agg_helper(raw$nlcor_xy_t3[,1]), agg_helper(raw$nlcor_xy_t3[,2]), agg_helper(raw$nlcor_xy_t3[,3]))
    }
    if (!is.null(raw$diff_arr)) {
      agg$diff_matrix <- apply(raw$diff_arr, c(1,2), agg_helper)
    }

  } else if (obj_type == "dda_bagging_resdist") {
    agg$agostino.target.statistic <- agg_helper(raw$agost_tar_stat)
    agg$agostino.target.z <- agg_helper(raw$agost_tar_z)
    agg$agostino.alternative.statistic <- agg_helper(raw$agost_alt_stat)
    agg$agostino.alternative.z <- agg_helper(raw$agost_alt_z)
    agg$anscombe.target.statistic <- agg_helper(raw$anscom_tar_stat)
    agg$anscombe.target.z <- agg_helper(raw$anscom_tar_z)
    agg$anscombe.alternative.statistic <- agg_helper(raw$anscom_alt_stat)
    agg$anscombe.alternative.z <- agg_helper(raw$anscom_alt_z)

    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS3", "RCC", "RHS4")) {
      if (!is.null(raw[[k]])) {
        new_agg <- apply(raw[[k]], 2, agg_helper)
        # keep the aggregated p-value computed in dda.bagging (column 3
        # of the 5-column skewness/kurtosis difference results)
        if (ncol(raw[[k]]) == 5 && length(agg[[k]]) >= 3) new_agg[3] <- agg[[k]][3]
        agg[[k]] <- new_agg
      }
    }

  } else if (obj_type == "dda_bagging_vardist") {
    agg$agostino.predictor.statistic.skew <- agg_helper(raw$agost_pre_stat)
    agg$agostino.predictor.statistic.z <- agg_helper(raw$agost_pre_z)
    agg$agostino.outcome.statistic.skew <- agg_helper(raw$agost_out_stat)
    agg$agostino.outcome.statistic.z <- agg_helper(raw$agost_out_z)
    agg$anscombe.predictor.statistic.kurt <- agg_helper(raw$anscom_pre_stat)
    agg$anscombe.predictor.statistic.z <- agg_helper(raw$anscom_pre_z)
    agg$anscombe.outcome.statistic.kurt <- agg_helper(raw$anscom_out_stat)
    agg$anscombe.outcome.statistic.z <- agg_helper(raw$anscom_out_z)

    for(k in c("skewdiff", "kurtdiff", "cor12diff", "cor13diff", "RHS", "RCC", "Rtanh")) {
      if (!is.null(raw[[k]])) {
        new_agg <- apply(raw[[k]], 2, agg_helper)
        # keep the aggregated p-value if a 5-column result is ever present
        if (ncol(raw[[k]]) == 5 && length(agg[[k]]) >= 3) new_agg[3] <- agg[[k]][3]
        agg[[k]] <- new_agg
      }
    }
  }

  object$aggregated_stats <- agg
  object$agg_stat_used <- agg_stat
  return(object)
}

#' @title Print Methods for Bootstrap Aggregated DDA Objects
#'
#' @description \code{print} returns aggregated test statistics and summaries
#' of bootstrap aggregated Direction Dependence Analysis (DDA) objects. The
#' function supports independence properties (obtained from
#' \code{dda.indep}), residual distributions (obtained from
#' \code{dda.resdist}), variable distributions (obtained from
#' \code{dda.vardist}).
#'
#' @param x An object of class \code{dda_bagging_indep},
#'   \code{dda_bagging_resdist}, or \code{dda_bagging_vardist}.
#' @param agg_stat Character. Specifies the method used for aggregating test
#'   statistics and coefficients across bootstrap samples. Must be one of the
#'   following specifications \code{c("mean", "median", "trimmed",
#'   "winsorized", "midhinge", "tukey")}. If \code{NULL}, the method
#'   established in \code{dda.bagging()} is used.
#' @param trim_prob Numeric. Proportion of observations to be trimmed on each
#'   side of the sampling distribution when \code{agg.stat = "trimmed"}
#'   (default: 0.10).
#' @param win_prob Numeric. Proportion of observations to be winsorized on
#'   each side of the sampling distribution when
#'   \code{agg.stat = "winsorized"} (default: 0.10).
#' @param digits Integer. Number of digits used for rounding (default: 4).
#' @param alpha Numeric. Significance level used for causal model selection
#'   (default: 0.05).
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Returns a summary of the bootstrap aggregated DDA object.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' base_model <- dda.indep(y ~ x, pred = "x", data = d, B = 20,
#'                          hetero = TRUE, nlfun = 2, diff = TRUE)
#' bagged <- dda.bagging(base_model, data = d, iter = 5,
#'                        agg_stat = "mean", progress = FALSE)
#'
#' # Default print
#' print(bagged)
#'
#' \donttest{
#' # Override aggregation method at print time
#' print(bagged, agg_stat = "trimmed", trim_prob = 0.05)
#'
#' # Print aggregated OLS summaries
#' summary_ols(bagged)
#' }
#' @export
#' @rdname print.dda_bagging
#' @method print dda_bagging_indep

print.dda_bagging_indep <- function(x,
                                    agg_stat = NULL,
                                    trim_prob = 0.10,
                                    win_prob = 0.10,
                                    digits = 4,
                                    alpha = 0.05,
                                    ...) {
  # Rename internal variable to 'x' for the print generic
  object <- reaggregate_bagging(x, agg_stat, trim_prob, win_prob)
  stats <- object$aggregated_stats

  # Try to find variable names (default y, x if missing)
  varnames <- NULL
  if (!is.null(stats$var.names) && length(stats$var.names) == 2) {
    varnames <- stats$var.names
  } else {
    varnames <- c("y", "x")
  }

  cat("\nBOOTSTRAP AGGREGATED DDA: Independence Properties\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n")
  cat("Aggregation method:", object$agg_stat_used, "\n\n")

  # Skip values that are NULL, NA, or NaN
  should_print <- function(val) {
    !is.null(val) && !is.na(val) && !is.nan(val)
  }

  # -------------------- Target Model (yx) --------------------
  cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n\n")

  # Omnibus Tests
  print_hsic <- should_print(stats$hsic_yx_stat)
  print_dcor <- should_print(stats$dcor_yx_stat)

  if (print_hsic || print_dcor) {
    cat("Omnibus Independence Tests:\n")
    if (print_hsic) {
      cat(paste0("HSIC = ", signif(as.numeric(stats$hsic_yx_stat), digits),
                 ", agg. p-value = ", signif(as.numeric(stats$hsic_yx_pval), digits)), "\n")
    }
    if (print_dcor) {
      cat(paste0("dCor = ", signif(as.numeric(stats$dcor_yx_stat), digits),
                 ", agg. p-value = ", signif(as.numeric(stats$dcor_yx_pval), digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Target)
  if (!is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 2) {
    cat("Homoscedasticity Tests:\n")
    cat(paste0("Standard Breusch-Pagan test: BP = ", signif(as.numeric(stats$breusch_pagan[[1]]$statistic), digits),
               ", df = ", signif(as.numeric(stats$breusch_pagan[[1]]$parameter), digits),
               ", agg. p-value = ", signif(as.numeric(stats$breusch_pagan[[1]]$p.value), digits)), "\n")
    cat(paste0("Robust Breusch-Pagan test:   BP = ", signif(as.numeric(stats$breusch_pagan[[2]]$statistic), digits),
               ", df = ", signif(as.numeric(stats$breusch_pagan[[2]]$parameter), digits),
               ", agg. p-value = ", signif(as.numeric(stats$breusch_pagan[[2]]$p.value), digits)), "\n")
    cat("\n")
  }

  # Non-linear Correlation (Target)
  if (!is.null(stats$nlcor.yx)) {
    nl <- stats$nlcor.yx
    func <- nl$func
    # Construct labels based on variable names and function
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation\n"))
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", varnames[2], "]"),
        paste0("Cor[", "r_", varnames[1], ", ", func, "(", varnames[2], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[1], "), ", func, "(", varnames[2], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func, "\n"))
      labs <- c(
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "]"),
        paste0("Cor[r_", varnames[1], ", ", varnames[2], "^", func, "]"),
        paste0("Cor[r_", varnames[1], "^", func, ", ", varnames[2], "^", func, "]")
      )
    }

    tab <- rbind(as.numeric(nl$t1), as.numeric(nl$t2), as.numeric(nl$t3))
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # -------------------- Alternative Model (xy) --------------------
  cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n\n")

  # Omnibus Tests
  print_hsic_alt <- should_print(stats$hsic_xy_stat)
  print_dcor_alt <- should_print(stats$dcor_xy_stat)

  if (print_hsic_alt || print_dcor_alt) {
    cat("Omnibus Independence Tests:\n")
    if (print_hsic_alt) {
      cat(paste0("HSIC = ", signif(as.numeric(stats$hsic_xy_stat), digits),
                 ", agg. p-value = ", signif(as.numeric(stats$hsic_xy_pval), digits)), "\n")
    }
    if (print_dcor_alt) {
      cat(paste0("dCor = ", signif(as.numeric(stats$dcor_xy_stat), digits),
                 ", agg. p-value = ", signif(as.numeric(stats$dcor_xy_pval), digits)), "\n")
    }
    cat("\n")
  }

  # Homoscedasticity Tests (Alternative)
  if (!is.null(stats$breusch_pagan) && length(stats$breusch_pagan) >= 4) {
    cat("Homoscedasticity Tests:\n")
    cat(paste0("Standard Breusch-Pagan test: BP = ", signif(as.numeric(stats$breusch_pagan[[3]]$statistic), digits),
               ", df = ", signif(as.numeric(stats$breusch_pagan[[3]]$parameter), digits),
               ", agg. p-value = ", signif(as.numeric(stats$breusch_pagan[[3]]$p.value), digits)), "\n")
    cat(paste0("Robust Breusch-Pagan test:   BP = ", signif(as.numeric(stats$breusch_pagan[[4]]$statistic), digits),
               ", df = ", signif(as.numeric(stats$breusch_pagan[[4]]$parameter), digits),
               ", agg. p-value = ", signif(as.numeric(stats$breusch_pagan[[4]]$p.value), digits)), "\n")
    cat("\n")
  }

  # Non-linear Correlation (Alternative)
  if (!is.null(stats$nlcor.xy)) {
    nl <- stats$nlcor.xy
    func <- nl$func
    if (is.na(suppressWarnings(as.numeric(func)))) {
      cat(paste("Non-linear Correlation Tests:", func, "Transformation\n"))
      labs <- c(
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", varnames[1], "]"),
        paste0("Cor[", "r_", varnames[2], ", ", func, "(", varnames[1], ")]"),
        paste0("Cor[", func, "(", "r_", varnames[2], "), ", func, "(", varnames[1], ")]")
      )
    } else {
      cat(paste("Non-linear Correlation Tests: Power Transformation using", func, "\n"))
      labs <- c(
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "]"),
        paste0("Cor[r_", varnames[2], ", ", varnames[1], "^", func, "]"),
        paste0("Cor[r_", varnames[2], "^", func, ", ", varnames[1], "^", func, "]")
      )
    }

    tab <- rbind(as.numeric(nl$t1), as.numeric(nl$t2), as.numeric(nl$t3))
    colnames(tab) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    rownames(tab) <- labs
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  # Difference Statistics CIs
  if (!is.null(stats$diff_matrix)) {
    boot_type <- "Percentile"
    if (!is.null(object$bagged_results[[1]]$boot.args)) {
      type_arg <- object$bagged_results[[1]]$boot.args[1]
      if (!is.null(type_arg)) {
        if(type_arg == "bca") boot_type <- "BCa"
        if(type_arg == "perc") boot_type <- "Percentile"
      }
    }
    conf_lev <- "95%"
    if (!is.null(object$bagged_results[[1]]$boot.args)) {
      lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
      if(!is.na(lev)) conf_lev <- paste0(lev*100, "%")
    }

    cat(paste0(conf_lev, " ", boot_type, " Bootstrap CIs for Difference Statistics (", object$n_valid_iterations, " samples):\n"))

    dmat <- stats$diff_matrix
    # ensure it is strictly numeric and un-named
    dmat_clean <- matrix(as.numeric(dmat), nrow=nrow(dmat), ncol=ncol(dmat))
    rownames(dmat_clean) <- rownames(dmat)
    colnames(dmat_clean) <- colnames(dmat)

    print.default(format(dmat_clean, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  cat("---\n")
  cat(paste("Note: Difference statistics > 0 suggest", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}

#' @export
#' @rdname print.dda_bagging
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' base_rd <- dda.resdist(y ~ x, pred = "x", data = d, B = 20)
#' bagged_rd <- dda.bagging(base_rd, data = d, iter = 5,
#'                           progress = FALSE)
#' print(bagged_rd)
#' \donttest{
#' print(bagged_rd, agg_stat = "median")
#' }
#' @method print dda_bagging_resdist
print.dda_bagging_resdist <- function(x, agg_stat = NULL, trim_prob = 0.10, win_prob = 0.10, digits = 4, ...) {
  object <- reaggregate_bagging(x, agg_stat, trim_prob, win_prob)
  stats <- object$aggregated_stats
  varnames <- if (!is.null(stats$var.names)) stats$var.names else c("target", "alternative")

  cat("\nBOOTSTRAP AGGREGATED DDA: Residual Distributions\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n")
  cat("Aggregation method:", object$agg_stat_used, "\n\n")

  cat("Skewness and kurtosis tests:\n")

  skew_row <- as.numeric(c(stats$agostino.target.statistic, stats$agostino.target.z, stats$agostino.target.p.value,
                           stats$agostino.alternative.statistic, stats$agostino.alternative.z, stats$agostino.alternative.p.value))
  kurt_row <- as.numeric(c(stats$anscombe.target.statistic, stats$anscombe.target.z, stats$anscombe.target.p.value,
                           stats$anscombe.alternative.statistic, stats$anscombe.alternative.z, stats$anscombe.alternative.p.value))

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", " z-value", " Pr(>|z|)", "alternative", " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = digits, scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  cat("\n")

  # Determine bootstrap info
  boot_type <- "Percentile"
  conf_lev <- "95%"
  if (!is.null(object$bagged_results[[1]]$boot.args)) {
    type_arg <- object$bagged_results[[1]]$boot.args[1]
    if (!is.null(type_arg)) {
      if(type_arg == "bca") boot_type <- "BCa"
      if(type_arg == "perc") boot_type <- "Percentile"
    }
    lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
    if(!is.na(lev)) conf_lev <- paste0(lev*100, "%")
  }

  cat(paste0("Skewness and kurtosis difference tests and ", conf_lev, " ", boot_type, " bootstrap CIs:\n\n"))

  citests <- rbind(as.numeric(stats$skewdiff), as.numeric(stats$kurtdiff))
  rownames(citests) <- c("Skewness", "Kurtosis")
  if (ncol(citests) == 5) colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
  if (ncol(citests) == 3) colnames(citests) <- c("diff", "lower", "upper")
  print.default(format(citests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for joint higher moment differences:\n"))
  joint_mat <- rbind(as.numeric(stats$cor12diff), as.numeric(stats$RHS3), as.numeric(stats$cor13diff), as.numeric(stats$RHS4), as.numeric(stats$RCC))
  rownames(joint_mat) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
  colnames(joint_mat) <- c("estimate", "lower", "upper")
  print.default(format(joint_mat, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\n---\n")
  probtrans <- if(!is.null(stats$probtrans)) stats$probtrans else FALSE
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], "\n"))
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], "\n"))
  if(isFALSE(probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], "\n"))
  } else {
    cat(paste("      Under prob.trans = TRUE, skewness and kurtosis differences < 0 and\n",
              "      co-skewness and co-kurtosis differences > 0 suggest", varnames[2], "->", varnames[1], "\n"))
  }

  invisible(object)
}

#' @export
#' @rdname print.dda_bagging
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' base_vd <- dda.vardist(y ~ x, pred = "x", data = d, B = 20)
#' bagged_vd <- dda.bagging(base_vd, data = d, iter = 5,
#'                           progress = FALSE)
#' print(bagged_vd)
#' \donttest{
#' print(bagged_vd, agg_stat = "winsorized", win_prob = 0.10)
#' }
#' @method print dda_bagging_vardist
print.dda_bagging_vardist <- function(x, agg_stat = NULL, trim_prob = 0.10, win_prob = 0.10, digits = 4, ...) {
  object <- reaggregate_bagging(x, agg_stat, trim_prob, win_prob)
  stats <- object$aggregated_stats
  varnames <- if (!is.null(stats$var.names)) stats$var.names else c("Outcome", "Predictor")

  cat("\nBOOTSTRAP AGGREGATED DDA: Variable Distributions\n")
  cat("Number of Bootstrap Samples:", object$n_valid_iterations, "\n")
  cat("Aggregation method:", object$agg_stat_used, "\n\n")

  cat("Skewness and kurtosis tests:\n")

  skew_row <- as.numeric(c(stats$agostino.outcome.statistic.skew, stats$agostino.outcome.statistic.z, stats$agostino.outcome.p.value,
                           stats$agostino.predictor.statistic.skew, stats$agostino.predictor.statistic.z, stats$agostino.predictor.p.value))
  kurt_row <- as.numeric(c(stats$anscombe.outcome.statistic.kurt, stats$anscombe.outcome.statistic.z, stats$anscombe.outcome.p.value,
                           stats$anscombe.predictor.statistic.kurt, stats$anscombe.predictor.statistic.z, stats$anscombe.predictor.p.value))

  sigtests <- rbind(skew_row, kurt_row)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
  print.default(format(sigtests, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\n")

  # Determine bootstrap info
  boot_type <- "Percentile"
  conf_lev <- "95%"
  if (!is.null(object$bagged_results[[1]]$boot.args)) {
    type_arg <- object$bagged_results[[1]]$boot.args[1]
    if (!is.null(type_arg)) {
      if(type_arg == "bca") boot_type <- "BCa"
      if(type_arg == "perc") boot_type <- "Percentile"
    }
    lev <- as.numeric(object$bagged_results[[1]]$boot.args[2])
    if(!is.na(lev)) conf_lev <- paste0(lev*100, "%")
  }

  cat(paste0(conf_lev, " ", boot_type, " bootstrap CIs for higher moment differences:\n"))
  citests <- rbind(as.numeric(stats$skewdiff), as.numeric(stats$kurtdiff))
  rownames(citests) <- c("Skewness", "Kurtosis")
  colnames(citests) <- c("diff", "lower", "upper")
  print.default(format(citests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for differences in higher-order correlations:\n"))
  hoctests <- rbind(as.numeric(stats$cor12diff), as.numeric(stats$cor13diff))
  rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]")
  colnames(hoctests) <- c("estimate", "lower", "upper")
  print.default(format(hoctests, digits = digits), print.gap = 2L, quote = FALSE)

  cat(paste0("\n", conf_lev, " ", boot_type, " bootstrap CIs for Likelihood Ratio approximations:\n"))
  LRtests <- rbind(as.numeric(stats$RHS), as.numeric(stats$Rtanh), as.numeric(stats$RCC))
  rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
  colnames(LRtests) <- c("estimate", "lower", "upper")
  print.default(format(LRtests, digits = digits), print.gap = 2L, quote = FALSE)

  cat("\n---\n")
  cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], "\n"))

  invisible(object)
}
