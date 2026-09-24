# Helper: Largest Remainder Method for rounding proportions to sum exactly to 1
#' @noRd
round_preserve_sum <- function(x, digits = 2) {
  if (sum(x, na.rm = TRUE) == 0) return(x)

  multiplier <- 10^digits
  scaled <- x * multiplier
  rounded <- floor(scaled)
  remainder <- scaled - rounded

  diff <- round(multiplier - sum(rounded))

  if (diff > 0 && diff <= length(x)) {
    # Add 1 to the elements with the largest decimal remainders
    idx <- order(remainder, decreasing = TRUE)[1:diff]
    rounded[idx] <- rounded[idx] + 1
  }

  return(rounded / multiplier)
}

# Internal helper to print decision summary
#' @noRd
print_bagging_decisions <- function(object, show = NULL, moment = NULL, type = "generic", digits = 2) {

  # --- Header ---
  header_type <- switch(type,
                        "indep" = "Independence Properties",
                        "resdist" = "Residual Distributions",
                        "vardist" = "Variable Distributions",
                        "Generic")

  cat("\n")
  cat(paste0("BOOTSTRAP AGGREGATED DDA: ", header_type), "\n")
  cat(paste0("Number of Bootstrap Samples: ", object$n_valid_iterations), "\n\n")

  decisions <- object$decision_percentages
  if (length(decisions) == 0) {
    cat("No decision proportions available.\n")
    return()
  }

  # --- Mappings ---
  alias_map <- list(
    "skew"   = c("dec_agost", "dec_skewdiff"),
    "kurt"   = c("dec_anscom", "dec_kurtdiff"),
    "coskew" = c("dec_cor12diff", "dec_RHS", "dec_RHS3"),
    "cokurt" = c("dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh"),
    "hsic"   = c("hsic", "diff_hsic"),
    "dcor"   = c("dcor", "diff_dcor"),
    "mi"     = c("diff_mi"),
    "bp"     = c("dec_bp"),
    "nlcor"  = c("dec_nl.min")
  )

  keys_m3 <- c("dec_agost", "dec_skewdiff", "dec_cor12diff", "dec_RHS", "dec_RHS3")
  keys_m4 <- c("dec_anscom", "dec_kurtdiff", "dec_cor13diff", "dec_RCC", "dec_RHS4", "dec_Rtanh")

  # Updated labels to match Docx track changes
  label_map <- list(
    "hsic"          = "HSIC",
    "dcor"          = "dCor",
    "diff_hsic"     = "HSIC Difference",
    "diff_dcor"     = "dCor Difference",
    "diff_mi"       = "MI Difference",
    "dec_bp"        = "Breusch-Pagan",
    "dec_nl.min"    = "Non-linear Correlation",
    "dec_agost"     = "Separate D'Agostino Tests",
    "dec_anscom"    = "Separate Anscombe-Glynn Tests",
    "dec_skewdiff"  = "Skewness Difference",
    "dec_kurtdiff"  = "Kurtosis Difference",
    "dec_cor12diff" = "Co-Skewness Difference",
    "dec_cor13diff" = "Co-Kurtosis Difference",
    "dec_RHS"       = "Hyvarinen-Smith Co-Skewness Difference",
    "dec_RHS3"      = "Hyvarinen-Smith Co-Skewness Difference",
    "dec_RHS4"      = "Hyvarinen-Smith Co-Kurtosis Difference",
    "dec_RCC"       = "Chen-Chan Co-Kurtosis Difference",
    "dec_Rtanh"     = "Hyvarinen-Smith tanh Difference"
  )

  all_keys <- names(decisions)
  keys_to_show <- c()

  # --- Logic for 'show' (including "all" option) and 'moment' ---
  if (!is.null(show)) {
    if ("all" %in% tolower(show)) {
      keys_to_show <- all_keys
    } else {
      expanded_show <- c()
      for (s in show) {
        if (s %in% names(alias_map)) {
          expanded_show <- c(expanded_show, alias_map[[s]])
        } else {
          expanded_show <- c(expanded_show, s)
        }
      }
      keys_to_show <- intersect(expanded_show, all_keys)
    }
  }

  moment_keys <- c()
  if (!is.null(moment)) {
    if (3 %in% moment) moment_keys <- c(moment_keys, keys_m3)
    if (4 %in% moment) moment_keys <- c(moment_keys, keys_m4)
    moment_keys <- intersect(moment_keys, all_keys)
  }

  if (is.null(show) && is.null(moment)) {
    keys_to_show <- all_keys
  } else if (is.null(show) && !is.null(moment)) {
    keys_to_show <- moment_keys
  } else if (!is.null(show) && !is.null(moment)) {
    keys_to_show <- unique(c(keys_to_show, moment_keys))
  }

  keys_to_show <- intersect(all_keys, keys_to_show)

  if (length(keys_to_show) == 0) {
    if(!is.null(show) || !is.null(moment)) {
      cat("No statistics matched the specified filters.\n")
    }
    return()
  }

  # Format string for dynamic decimal places
  fmt <- paste0("%.", digits, "f")

  # --- Print Decision Tables ---
  for (dname in keys_to_show) {
    prop <- decisions[[dname]]

    if (all(is.nan(prop)) || all(is.na(prop))) next

    display_title <- if (!is.null(label_map[[dname]])) label_map[[dname]] else dname
    cat(display_title, "\n")

    # Extract exact values; Confounding is a real category for the separate
    # independence tests (both directions significant), not a relabeling of
    # Undecided
    t_val <- if ("Target" %in% names(prop)) prop["Target"] else 0
    a_val <- if ("Alternative" %in% names(prop)) prop["Alternative"] else 0
    c_val <- if ("Confounding" %in% names(prop)) prop["Confounding"] else 0
    u_val <- if ("Undecided" %in% names(prop)) prop["Undecided"] else 0

    # Apply Largest Remainder Method rounding dynamically based on user digits
    if (type == "indep") {
      vals_rounded <- round_preserve_sum(c(t_val, a_val, c_val, u_val), digits = digits)
      df_print <- data.frame(
        Target      = sprintf(fmt, vals_rounded[1]),
        Alternative = sprintf(fmt, vals_rounded[2]),
        Confounding = sprintf(fmt, vals_rounded[3]),
        Undecided   = sprintf(fmt, vals_rounded[4]),
        check.names = FALSE
      )
    } else {
      vals_rounded <- round_preserve_sum(c(t_val, a_val, u_val), digits = digits)
      df_print <- data.frame(
        Target      = sprintf(fmt, vals_rounded[1]),
        Alternative = sprintf(fmt, vals_rounded[2]),
        Undecided   = sprintf(fmt, vals_rounded[3]),
        check.names = FALSE
      )
    }

    print(df_print, row.names = FALSE)
    cat("\n")
  }

  # --- Footnote with Actual Variable Names ---
  stats <- object$aggregated_stats
  varnames <- NULL

  # Try to get variable names from aggregated stats
  if (!is.null(stats$var.names) && length(stats$var.names) == 2) {
    varnames <- stats$var.names
  } else {
    # Fallback defaults
    if (type == "indep") {
      varnames <- c("y", "x")
    } else if (type == "resdist") {
      varnames <- c("target", "alternative")
    } else {
      varnames <- c("Outcome", "Predictor")
    }
  }

  # Print footnote using actual variable names
  cat("---\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1]), "\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2]), "\n")
}

#' @title Summary Methods for Bootstrap Aggregated DDA Objects
#'
#' @description \code{summary} returns the proportion of model selection
#' decisions (i.e., target model, alternative model, or presence of hidden
#' confounding) across the bootstrap samples for bootstrap aggregated
#' Direction Dependence Analysis (DDA) objects.
#'
#' @param object An object of class \code{dda_bagging_indep}, \code{dda_bagging_vardist},
#'   or \code{dda_bagging_resdist}.
#' @param show Character vector specifying the DDA statistics used to
#'   perform model selection. The function accepts the following
#'   specifications \code{c("all", "hsic", "dcor", "mi", "bp", "nlcor",
#'   "skew", "kurt", "coskew", "cokurt")}. If \code{NULL},
#'   shows default statistics based on the object type.
#' @param moment Numeric vector indicating which moments to include in the summary
#'   (e.g., \code{c(3, 4)}). Applicable only to \code{dda_bagging_vardist} and \code{dda_bagging_resdist}.
#' @param digits Integer. Number of decimal places to print for proportions (default: 2).
#' @param ... Additional arguments passed to \code{summary}.
#'
#' @return Invisibly returns the original object.
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
#' bagged <- dda.bagging(base_model, iter = 10, data = d, progress = FALSE)
#'
#' ## proportion of model selection decisions across bootstrap samples
#' summary(bagged)
#'
#' \donttest{
#' ## restrict the table to selected statistics and moments
#' summary(bagged, show = c("skew", "coskew"), moment = 3)
#' }
#'
#' @export
#' @rdname summary.dda_bagging
#' @method summary dda_bagging_indep
summary.dda_bagging_indep <- function(object, show = NULL, digits = 2, ...) {
  print_bagging_decisions(object, show = show, type = "indep", digits = digits)
  invisible(object)
}

#' @export
#' @rdname summary.dda_bagging
#' @method summary dda_bagging_vardist
summary.dda_bagging_vardist <- function(object, show = NULL, moment = NULL, digits = 2, ...) {
  print_bagging_decisions(object, show = show, moment = moment, type = "vardist", digits = digits)
  invisible(object)
}

#' @export
#' @rdname summary.dda_bagging
#' @method summary dda_bagging_resdist
summary.dda_bagging_resdist <- function(object, show = NULL, moment = NULL, digits = 2, ...) {
  print_bagging_decisions(object, show = show, moment = moment, type = "resdist", digits = digits)
  invisible(object)
}
