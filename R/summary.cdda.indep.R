#' @title Summary of \code{cdda.indep} Class Objects
#'
#' @description \code{summary} returns test statistics from the
#' \code{cdda.indep} class object.
#'
#' @name summary.cdda.indep
#'
#' @param object An object of class \code{cdda.indep} when using
#'   \code{summary}.
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
  ## OLD
  ## n.resamples <- cdda.output[[1]][[1]]$boot.args[3]
  ## ci.level <- as.numeric(cdda.output[[1]][[1]]$boot.args[2]) * 100
  ##
  ## if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "bca") {
  ##   boot.print <- paste("\n", ci.level, "% BCa bootstrap CIs (", n.resamples, " resamples)", "\n\n", sep = "")
  ## } else if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "perc") {
  ##   boot.print <- paste("\n", ci.level, "% percentile bootstrap CIs (", n.resamples, " resamples)", "\n\n", sep = "")
  ## } else {
  ##   boot.print <- NULL
  ## }
  ## NEW 0.2.0
  boot.args <- cdda.output$cdda.diff$boot.args
  n.resamples <- boot.args[3]
  ci.level <- as.numeric(boot.args[2]) * 100

  if (!is.null(boot.args[1]) && boot.args[1] == "bca") {
    boot.print <- paste("\n", ci.level, "% BCa bootstrap CIs (", n.resamples, " resamples)", "\n\n", sep = "")
  } else if (!is.null(boot.args[1]) && boot.args[1] == "perc") {
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
    cdda.print.table(bptests)
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
    cdda.print.table(rbptests)
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
    cdda.print.table(nlsigtests)
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
    cdda.print.table(hsictests)
    cat("---", "\n", "\n")
  }

  ## OLD
  ## if (hsic.diff == TRUE){
  ##
  ##   ### HSIC Diff ### ----------------------------------------------------------
  ##   if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")
  ##
  ##   hsicdiff <- matrix(NA, length(mod.names), 6)
  ##
  ##   for (i in seq_along(mod.names)) {
  ##
  ##     tar.hsicd <- unlist(cdda.output[[1]][[i]]$out.diff[1, ])
  ##     alt.hsicd <- unlist(cdda.output[[2]][[i]]$out.diff[1, ])
  ##     hsicdiff[i, ] <- c(tar.hsicd, alt.hsicd)
  ##   }
  ##
  ##   rownames(hsicdiff) <- mod.names
  ##   colnames(hsicdiff) <- rep(c("diff", "lower", "upper"), 2)
  ##   hsicdiff <- round(hsicdiff, 3)
  ##
  ##   cat(boot.print)
  ##
  ##   cat(paste("Hilbert-Schmidt Independence Criterion Differences", "\n"))
  ##   cdda.print.table(hsicdiff)
  ##   cat("---", "\n", "\n")
  ## }
  ## NEW 0.2.0
  show.diff <- function(row, label, note = NULL){
    d <- cdda.output$cdda.diff$diff
    if (length(d) == 0) stop("Difference tests not found, set 'diff = TRUE'.")
    tab <- round(t(sapply(d, function(m) m[row, ])), 3)
    rownames(tab) <- mod.names
    cat(boot.print)
    cat(label, "\n")
    cat("Values > 0 suggest", cdda.output$df.original$pred.name, "->",
        cdda.output$df.original$response.name, "\n")
    print(tab, print.gap = 2L)
    if (!is.null(note)) cat(note, "\n")
    cat("---", "\n", "\n")
  }

  if (hsic.diff == TRUE) show.diff("HSIC", "Hilbert-Schmidt Independence Criterion Differences")

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
    cdda.print.table(dcor)
    cat("---", "\n", "\n")
  }

  ## OLD
  ## ## dcor.diff ## --------------------------------------------------------------
  ## if (dcor.diff == TRUE){
  ##
  ##   if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")
  ##
  ##   dcortests <- matrix(NA, length(mod.names), 6)
  ##
  ##   for (i in seq_along(mod.names)) {
  ##     tar.dcor.diff <- unlist(cdda.output[[1]][[i]]$out.diff[2, ])
  ##     alt.dcor.diff <- unlist(cdda.output[[2]][[i]]$out.diff[2, ])
  ##     dcortests[i, ] <- c(tar.dcor.diff, alt.dcor.diff)
  ##   }
  ##
  ##   rownames(dcortests) <- mod.names
  ##   colnames(dcortests) <- rep(c("diff", "lower", "upper"), 2)
  ##   dcortests <- round(dcortests, 3)
  ##
  ##   cat(boot.print)
  ##
  ##   cat(paste("Distance Correlation Differences", "\n"))
  ##   cdda.print.table(dcortests)
  ##   cat("---", "\n", "\n")
  ## }
  ## NEW 0.2.0
  if (dcor.diff == TRUE) show.diff("dCor", "Distance Correlation Differences")
  ## OLD
  ## ### Mutual Information ### ---------------------------------------------------
  ## if (mi.diff == TRUE){
  ##
  ##   if(is.null(cdda.output[[1]][[1]]$out.diff)) stop("Difference tests not found, set 'diff = TRUE'.")
  ##
  ##   mitests <- matrix(NA, length(mod.names), 6)
  ##
  ##   for (i in seq_along(mod.names)) {
  ##     tar.mi <- unlist(cdda.output[[1]][[i]]$out.diff[3, ])
  ##     alt.mi <- unlist(cdda.output[[2]][[i]]$out.diff[3, ])
  ##     mitests[i, ] <- c(tar.mi, alt.mi)
  ##   }
  ##
  ##   rownames(mitests) <- mod.names
  ##   colnames(mitests) <- rep(c("diff", "lower", "upper"), 2)
  ##   mitests <- round(mitests, 3)
  ##
  ##   cat(boot.print)
  ##
  ##   cat(paste("Mutual Information Difference", "\n"))
  ##   cdda.print.table(mitests)
  ##   cat("---", "\n", "\n")
  ##
  ## }
  ## NEW 0.2.0
  if (mi.diff == TRUE) show.diff("MI", "Mutual Information Differences")
}


#' @title Print a Target and Alternative Model Table
#'
#' @description Prints a matrix with target model columns followed by
#'   alternative model columns and places the model labels over their columns.
#'   Used by the summary methods of cdda.indep and cdda.vardist.
#'
#' @keywords internal
#' @noRd
cdda.print.table <- function(tab, gap = 2L, digits = 3L){

  k <- ncol(tab) / 2
  heads <- colnames(tab)

  # numbers with a fixed number of decimals, degrees of freedom as integers
  cells <- matrix("", nrow(tab), ncol(tab))
  for (j in seq_len(ncol(tab))) {
    if (heads[j] == "df") cells[, j] <- format(round(tab[, j]))
    else cells[, j] <- formatC(tab[, j], format = "f", digits = digits)
  }

  # each model label must fit over its columns
  w <- pmax(nchar(heads), apply(nchar(cells), 2, max))
  labs <- c("Target Model", "Alternative Model")
  for (m in 1:2) {
    cols <- if (m == 1) 1:k else (k + 1):(2 * k)
    span <- sum(w[cols]) + gap * (k - 1)
    if (span < nchar(labs[m])) w[min(cols)] <- w[min(cols)] + nchar(labs[m]) - span
  }
  span.tar <- sum(w[1:k]) + gap * (k - 1)

  rn <- formatC(c("", rownames(tab)), width = -max(nchar(rownames(tab))))
  sep <- strrep(" ", gap)
  join <- function(v) {
    v <- sprintf("%*s", w, v)
    paste0(paste(v[1:k], collapse = sep), strrep(" ", 2 * gap),
           paste(v[(k + 1):(2 * k)], collapse = sep))
  }
  lines <- paste0(rn, sep, c(join(heads), apply(cells, 1, join)))
  labels <- paste0(rn[1], sep, formatC(labs[1], width = -span.tar),
                   strrep(" ", 2 * gap), labs[2])
  rule <- strrep("-", max(nchar(c(lines, labels))))

  cat(rule, labels, rule, lines, sep = "\n")
}
