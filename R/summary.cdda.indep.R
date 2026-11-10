#' @title Summary of \code{cdda.indep} Class Objects
#'
#' @description \code{summary} returns test statistics from the
#' \code{cdda.indep} class object.
#'
#' @name summary.cdda.indep
#'
#' @param object An object of class \code{cdda.indep} when using
#'   \code{summary}.
#' @param hetero A logical value indicating whether separate homoscedasticity
#'   tests should be returned when using \code{summary}, default is
#'   \code{FALSE}.
#' @param hsic A logical value indicating whether separate HSIC tests should
#'   be returned when using \code{summary}, default is \code{TRUE}.
#' @param hsic.diff A logical value indicating whether HSIC difference
#'   statistics should be returned when using \code{summary}, default is
#'   \code{FALSE}.
#' @param dcor A logical value indicating whether separate Distance
#'   Correlation (dCor) tests should be returned when using \code{summary},
#'   default is \code{TRUE}.
#' @param dcor.diff A logical value indicating whether dCor difference
#'   statistics should be returned when using \code{summary}, default is
#'   \code{FALSE}.
#' @param mi.diff A logical value indicating whether Mutual Information (MI)
#'   difference statistics should be returned when using \code{summary},
#'   default is \code{FALSE}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @examples
#' summary(result, hetero = FALSE)
#'
#' @export
#' @rdname cdda.indep
#' @method summary cdda.indep
summary.cdda.indep <- function(object,
                               nlfun = FALSE,
                               hetero = FALSE,
                               hsic = TRUE,
                               hsic.diff = FALSE,
                               dcor = TRUE,
                               dcor.diff = FALSE,
                               mi.diff = FALSE,
                               ...
                               ){

  cdda.output <- object

  varnames <- cdda.output[[1]][[1]]$var.names
  mod.names <- names(cdda.output[[1]])

  varnames <- cdda.output[[1]][[1]]$var.names
  mod.names <- names(cdda.output[[1]])
  n.resamples <- cdda.output[[1]][[1]]$boot.args[3]
  ci.level <- as.numeric(cdda.output[[1]][[1]]$boot.args[2]) * 100

  if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "bca") {
    boot.print <- paste("\n", ci.level, "% BCa bootstrap CIs (", n.resamples, " resamples)", "\n\n", sep = "")
  } else if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "perc") {
    boot.print <- paste("\n", ci.level, "% percentile bootstrap CIs (", n.resamples, " resamples)", "\n\n", sep = "")
  } else {
    boot.print <- NULL
  }
  cat("\n")

  ### Indep Summary ### ----------------------------------------------------

  if (hetero == TRUE) {
    if(is.null(cdda.output[[1]][[1]]$breusch.pagan)){
      stop("Homoscedasticity tests not found. Use 'hetero = TRUE'.")
    }

    ### BP ### ------------------------------------------------------------
    bptests <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {

      tar.bp <- unlist(cdda.output[[1]][[i]]$breusch.pagan[[1]][c(1,2,4)])
      alt.bp <- unlist(cdda.output[[2]][[i]]$breusch.pagan[[1]][c(1,2,4)])
      bptests[i, ] <- c(tar.bp, alt.bp)
    }

    rownames(bptests) <- mod.names
    colnames(bptests) <- rep(c("X2-value", "df", "p-value"), 2)
    bptests <- round(bptests, 3)

    cat(paste("Homoscedasticity test Breusch-Pagan", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model        Alternative Model ", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    print(bptests)
    cat("---", "\n", "\n")

    ### Robust BP ### --------------------------------------------------------
    rbptests <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {

      tar.rbp <- unlist(cdda.output[[1]][[i]]$breusch.pagan[[2]][c(1,2,4)])
      alt.rbp <- unlist(cdda.output[[2]][[i]]$breusch.pagan[[2]][c(1,2,4)])
      rbptests[i, ] <- c(tar.rbp, alt.rbp)
    }

    rownames(rbptests) <- mod.names
    colnames(rbptests) <- rep(c("X2-value", "df", "p-value"), 2)
    rbptests <- round(rbptests, 3)

    cat(paste("Homoscedasticity test Robust Breusch-Pagan", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model        Alternative Model ", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    print(rbptests)
    cat("---", "\n", "\n")
  }

  if(nlfun == TRUE) {

    ### Non-linear Correlation Tests ### -------------------------------------

    if(is.null(cdda.output[[1]][[1]]$nlfun)) stop("Non-linear function is missing.")

    nlsigtests <- matrix(NA, length(mod.names), 8)

    for (i in seq_along(mod.names)) {

      tar.nl <- rbind(cdda.output[[1]][[i]]$nlcor.yx$t1,
                      cdda.output[[1]][[i]]$nlcor.yx$t2,
                      cdda.output[[1]][[i]]$nlcor.yx$t3)
      tar.nl <- tar.nl[which.min(tar.nl[,4]),]


      alt.nl <- rbind(cdda.output[[2]][[i]]$nlcor.yx$t1,
                      cdda.output[[2]][[i]]$nlcor.yx$t2,
                      cdda.output[[2]][[i]]$nlcor.yx$t3)
      alt.nl <- alt.nl[which.min(alt.nl[,4]),]

      nlsigtests[i, ] <- c(tar.nl, alt.nl)
    }

    rownames(nlsigtests) <- mod.names
    colnames(nlsigtests) <- rep(c("statistic", "t-value", "df", "p-value"), 2)
    nlsigtests <- round(nlsigtests, 3)

    cat(paste("Non-linear Correlation Tests", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model                  Alternative Model ", "\n"))
    cat(paste0("-----------------------------------------------------------------------", "\n"))
    print(nlsigtests)
    cat("---", "\n", "\n")
  }

  if (hsic == TRUE) {

    ### HSIC ### ------------------------------------------------------------
    if(is.null(cdda.output[[1]][[1]]$hsic.yx)) stop("Difference tests not found, set 'diff = TRUE'.")

    hsictests <- matrix(NA, length(mod.names), 4)

    for (i in seq_along(mod.names)) {

      tar.hsic <- c(cdda.output[[1]][[i]]$hsic.yx$statistic, cdda.output[[1]][[i]]$hsic.yx$p.value)
      alt.hsic <- c(cdda.output[[2]][[i]]$hsic.yx$statistic, cdda.output[[2]][[i]]$hsic.yx$p.value)
      hsictests[i, ] <- c(tar.hsic, alt.hsic)
    }

    rownames(hsictests) <- mod.names
    colnames(hsictests) <- rep(c("HSIC", "p-value"), 2)
    hsictests <- round(hsictests, 3)

    cat(paste("Hilbert-Schmidt Independence Criterion", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model                 Alternative Model ", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    print(hsictests)
    cat("---", "\n", "\n")
  }

  if (hsic.diff == TRUE){

    ### HSIC Diff ### ----------------------------------------------------------
    if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")

    hsicdiff <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {

      tar.hsicd <- unlist(cdda.output[[1]][[i]]$out.diff[1, ])
      alt.hsicd <- unlist(cdda.output[[2]][[i]]$out.diff[1, ])
      hsicdiff[i, ] <- c(tar.hsicd, alt.hsicd)
    }

    rownames(hsicdiff) <- mod.names
    colnames(hsicdiff) <- rep(c("diff", "lower", "upper"), 2)
    hsicdiff <- round(hsicdiff, 3)

    cat(boot.print)

    cat(paste("Hilbert-Schmidt Independence Criterion Differences", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model                 Alternative Model ", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    print(hsicdiff)
    cat("---", "\n", "\n")
  }

  ### Distance Correlation ### -------------------------------------------------

  ## dcor ## -------------------------------------------------------------------
  if (dcor == TRUE){

    dcor <- matrix(NA, length(mod.names), 4)

    for (i in seq_along(mod.names)) {
      tar.dcor <- unlist(c(cdda.output[[1]][[i]]$distance.cor.dcor.yx$statistic,
                           cdda.output[[1]][[i]]$distance.cor.dcor.yx$p.value))
      alt.dcor <- unlist(c(cdda.output[[2]][[i]]$distance.cor.dcor.yx$statistic,
                           cdda.output[[2]][[i]]$distance.cor.dcor.yx$p.value))
      dcor[i, ] <- c(tar.dcor, alt.dcor)
    }

    rownames(dcor) <- mod.names
    colnames(dcor) <- rep(c("dCor", "p-value"), 2)
    dcor <- round(dcor, 3)

    cat(paste("Distance Correlation", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model           Alternative Model ", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))

    print(dcor)
    cat("---", "\n", "\n")
  }

  ## dcor.diff ## --------------------------------------------------------------
  if (dcor.diff == TRUE){

    if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")

    dcortests <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {
      tar.dcor.diff <- unlist(cdda.output[[1]][[i]]$out.diff[2, ])
      alt.dcor.diff <- unlist(cdda.output[[2]][[i]]$out.diff[2, ])
      dcortests[i, ] <- c(tar.dcor.diff, alt.dcor.diff)
    }

    rownames(dcortests) <- mod.names
    colnames(dcortests) <- rep(c("diff", "lower", "upper"), 2)
    dcortests <- round(dcortests, 3)

    cat(boot.print)

    cat(paste("Distance Correlation Differences", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model           Alternative Model ", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))

    print(dcortests)
    cat("---", "\n", "\n")
  }
  ### Mutual Information ### ---------------------------------------------------
  if (mi.diff == TRUE){

    if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")

    mitests <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {
      tar.mi <- unlist(cdda.output[[1]][[i]]$out.diff[3, ])
      alt.mi <- unlist(cdda.output[[2]][[i]]$out.diff[3, ])
      mitests[i, ] <- c(tar.mi, alt.mi)
    }

    rownames(mitests) <- mod.names
    colnames(mitests) <- rep(c("diff", "lower", "upper"), 2)
    mitests <- round(mitests, 3)

    cat(boot.print)

    cat(paste("Mutual Information Difference", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    cat(paste0("            Target Model         Alternative Model ", "\n"))
    cat(paste0("---------------------------------------------------------------------", "\n"))
    print(mitests)
    cat("---", "\n", "\n")

  }
}
