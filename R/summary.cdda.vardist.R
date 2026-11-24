#' @title Summary of \code{cdda.vardist} Class Objects
#'
#' @description \code{summary} returns test statistics from the
#' \code{cdda.vardist} class object.
#'
#' @name summary.cdda.vardist
#'
#' @param object An object of class \code{cdda.vardist} when using
#'   \code{summary}.
#' @param skew A logical value indicating whether skewness differences and
#'   separate D'Agostino skewness tests should be returned when using
#'   \code{summary}, default is \code{TRUE}.
#' @param coskew A logical value indicating whether co-skewness differences
#'   (CoSk and HS.CoSk) should be returned when using \code{summary}, default is \code{FALSE}.
#' @param kurt A logical value indicating whether excess kurtosis differences
#'   and Anscombe-Glynn kurtosis tests should be returned when using
#'   \code{summary}, default is \code{TRUE}.
#' @param cokurt A logical value indicating whether co-kurtosis differences
#'   (CoKurt, CC.CoKurt, and Rtanh) should be returned when using \code{summary}, default is \code{FALSE}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @export
#' @rdname cdda.vardist
#' @method summary cdda.vardist


summary.cdda.vardist <- function(object,
                                 skew = TRUE,
                                 coskew = FALSE,
                                 kurt = TRUE,
                                 cokurt = FALSE,
                                 ...
                                 ){
  cdda.output <- object

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

  ## NEW 0.2.0: reading of the separate normality tests
  sep.note <- paste0("Note: A model is supported when its predictor deviates from normality\n",
                     "      and its outcome does not. The predictor is ", cdda.output$df.original$pred.name,
                     " in the target model\n",
                     "      and ", cdda.output$df.original$response.name, " in the alternative model.\n")

    if (skew == TRUE){

      ### ----------- Boot CI Loop Setup ------------ ###
      momdif.stats.skew <- matrix(NA, length(mod.names), 6)

      rownames(momdif.stats.skew) <- mod.names
      colnames(momdif.stats.skew) <-  rep(c("diff", "lower", "upper"), 2)

      ## Boot CI Loop ##
      for (i in seq_along(mod.names)) {

        momdif.stats.skew[i, ]<- c(cdda.output[[1]][[i]]$skewdiff,
                                   cdda.output[[2]][[i]]$skewdiff
                                  )
      }

      ### Agostino ### ---------------------------------------------------------
      ## OLD
      ## agos <- matrix(NA, length(mod.names), 6)
      ## rownames(agos) <- mod.names
      ## colnames(agos) <- rep(c("statistic", "z-value", "p-value"), 2)
      ##
      ## for (i in seq_along(mod.names)) {
      ##
      ##   agos[i, ]<- c(unlist(cdda.output[[1]][[i]]$agostino$predictor),
      ##                 unlist(cdda.output[[2]][[i]]$agostino$predictor)
      ##   )
      ##   agos[i, ] <- round(agos[i, ], 3)
      ## }
      ## NEW 0.2.0: predictor and outcome tests of both panels
      agos <- matrix(NA, 2 * length(mod.names), 6)
      rownames(agos) <- paste(rep(mod.names, each = 2), c("predictor", "outcome"))
      colnames(agos) <- rep(c("statistic", "z-value", "p-value"), 2)

      for (i in seq_along(mod.names)) {
        for (v in c("predictor", "outcome")) {
          agos[paste(mod.names[i], v), ] <- c(unlist(cdda.output[[1]][[i]]$agostino[[v]]),
                                              unlist(cdda.output[[2]][[i]]$agostino[[v]]))
        }
      }
      agos <- round(agos, 3)

      cat(paste("D'Agostino Skewness Tests", "\n"))
      cdda.print.table(agos)
      cat(sep.note)

      cat(boot.print)
      cat(paste("Skewness differences"), "\n")
      cdda.print.table(round(momdif.stats.skew, 3))

      cat("---", "\n")
    }

  if (coskew == TRUE){

    if(is.null(cdda.output[[1]][[1]]$cor12diff)) { stop("CoSk differences not found.") }
    hoctests.skew <-  matrix(NA, length(mod.names), 6)

    for(i in seq_along(mod.names)){
      hoctests.skew[i, ] <- c(cdda.output[[1]][[i]]$cor12diff,
                              cdda.output[[2]][[i]]$cor12diff)
      hoctests.skew[i, ] <- round(hoctests.skew[i, ], 3)
    }

    cat(boot.print)

    cat("CoSk differences", "\n")

    rownames(hoctests.skew) <- mod.names
    colnames(hoctests.skew) <- rep(c("diff", "lower", "upper"), 2)

    cdda.print.table(hoctests.skew)
    cat("\n")

    if(cdda.output[[1]][[1]]$boot.warning) { cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n", "        CoKurt should also be computed for the model", varnames[1], "->", varnames[2], "\n") }

    ### HS.CoSk ### ------------------------------------------------------------
    rhs <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {
      rhs[i, ] <- c(cdda.output[[1]][[i]]$RHS, cdda.output[[2]][[i]]$RHS)
      rhs[i, ] <- round(rhs[i, ], 3)
    }

    rownames(rhs) <- mod.names
    colnames(rhs) <- rep(c("diff", "lower", "upper"), 2)

    cat("HS.CoSk differences", "\n")
    cdda.print.table(rhs)

    cat("---", "\n")
  }

    if (kurt == TRUE){

      ### Anscombe ### ---------------------------------------------------------
      ## OLD
      ## ansc <- matrix(NA, length(mod.names), 6)
      ## rownames(ansc) <- mod.names
      ## colnames(ansc) <- rep(c("statistic", "z-value", "p-value"), 2)
      ##
      ## for (i in seq_along(mod.names)) {
      ##
      ##   ansc[i, ]<- c(unlist(cdda.output[[1]][[i]]$anscombe$predictor),
      ##                 unlist(cdda.output[[2]][[i]]$anscombe$predictor)
      ##                )
      ##   ansc[i, ] <- round(ansc[i, ], 3)
      ##  }
      ## NEW 0.2.0: predictor and outcome tests of both panels
      ansc <- matrix(NA, 2 * length(mod.names), 6)
      rownames(ansc) <- paste(rep(mod.names, each = 2), c("predictor", "outcome"))
      colnames(ansc) <- rep(c("statistic", "z-value", "p-value"), 2)

      for (i in seq_along(mod.names)) {
        for (v in c("predictor", "outcome")) {
          ansc[paste(mod.names[i], v), ] <- c(unlist(cdda.output[[1]][[i]]$anscombe[[v]]),
                                              unlist(cdda.output[[2]][[i]]$anscombe[[v]]))
        }
      }
      ansc <- round(ansc, 3)

      cat(paste("Anscombe-Glynn Kurtosis Tests", "\n"))
      cdda.print.table(ansc)
      cat(sep.note)

      cat("---", "\n")

    ### ----------- Boot CI Loop Setup ------------ ###

    momdif.stats.kurt <- matrix(NA, length(mod.names), 6)

    rownames(momdif.stats.kurt) <- mod.names
    colnames(momdif.stats.kurt) <- rep(c("diff", "lower", "upper"), 2)

    ## Boot CI Loop ##
    for (i in seq_along(mod.names)) {

      momdif.stats.kurt[i, ]<- c(cdda.output[[1]][[i]]$kurtdiff, cdda.output[[2]][[i]]$kurtdiff)
     }

    cat(boot.print)

    cat(paste("Excess Kurtosis differences"), "\n")
    cdda.print.table(round(momdif.stats.kurt, 3))

    cat("---", "\n")
  }

    if (cokurt == TRUE){

    ### CoKurt ### -------------------------------------------------------------

    hoctests.kurt <-  matrix(NA, length(mod.names), 6)

    for(i in seq_along(mod.names)){
      hoctests.kurt[i, ] <- c(cdda.output[[1]][[i]]$cor13diff, cdda.output[[2]][[i]]$cor13diff)
      hoctests.kurt[i, ] <- round(hoctests.kurt[i, ], 3)
    }

    colnames(hoctests.kurt) <- rep(c("diff", "lower", "upper"), 2)
    rownames(hoctests.kurt) <- mod.names

    cat(boot.print)

    cat("CoKurt differences", "\n")

    cdda.print.table(hoctests.kurt)
    cat("---", "\n")

    if(cdda.output[[1]][[1]]$boot.warning) { cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n", "        CoKurt should also be computed for the model", varnames[1], "->", varnames[2], "\n") }

    ### CC.CoKurt ### ----------------------------------------------------------

    rcc <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {
      rcc[i, ] <- c(cdda.output[[1]][[i]]$RCC, cdda.output[[2]][[i]]$RCC)
      rcc[i, ] <- round(rcc[i, ], 3)
    }

    rownames(rcc) <- mod.names
    colnames(rcc) <- rep(c("diff", "lower", "upper"), 2)

    cat(boot.print)

    cat("CC.CoKurt differences", "\n")
    cdda.print.table(rcc)

    cat("---", "\n")

    ### Rtanh ### --------------------------------------------------------------

    rtanh <- matrix(NA, length(mod.names), 6)

    for (i in seq_along(mod.names)) {

      rtanh[i, ] <- c(cdda.output[[1]][[i]]$Rtanh, cdda.output[[2]][[i]]$Rtanh)
      rtanh[i, ] <- round(rtanh[i, ], 3)
    }

    rownames(rtanh) <- mod.names
    colnames(rtanh) <- rep(c("diff", "lower", "upper"), 2)

    cat(boot.print)

    cat("Rtanh differences", "\n")
    cdda.print.table(rtanh)

    cat("---", "\n")
    }
}



