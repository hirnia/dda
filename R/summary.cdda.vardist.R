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
#'   should be returned when using \code{summary}, default is \code{FALSE}.
#' @param kurt A logical value indicating whether excess kurtosis differences
#'   and Anscombe-Glynn kurtosis tests should be returned when using
#'   \code{summary}, default is \code{TRUE}.
#' @param cokurt A logical value indicating whether co-kurtosis differences
#'   should be returned when using \code{summary}, default is \code{FALSE}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @examples
#' summary(result, skew = FALSE, kurt = FALSE, coskew = TRUE)
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
  mod_names <- names(cdda.output[[1]])
  n_resamples <- cdda.output[[1]][[1]]$boot.args[3]
  ci.level <- as.numeric(cdda.output[[1]][[1]]$boot.args[2]) * 100

  if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "bca") {
    boot_print <- paste("\n", ci.level, "% BCa bootstrap CIs (", n_resamples, " resamples)", "\n\n", sep = "")
  } else if (!is.null(cdda.output[[1]][[1]]$boot.args[1]) && cdda.output[[1]][[1]]$boot.args[1] == "perc") {
    boot_print <- paste("\n", ci.level, "% percentile bootstrap CIs (", n_resamples, " resamples)", "\n\n", sep = "")
  } else {
    boot_print <- NULL
  }
  cat("\n")

    if (skew == TRUE){

      ### ----------- Boot CI Loop Setup ------------ ###
      momdif.stats.skew <- matrix(NA, length(mod_names), 6)

      rownames(momdif.stats.skew) <- mod_names
      colnames(momdif.stats.skew) <-  rep(c("diff", "lower", "upper"), 2)

      ## Boot CI Loop ##
      for (i in seq_along(mod_names)) {

        momdif.stats.skew[i, ]<- c(cdda.output[[1]][[i]]$skewdiff,
                                   cdda.output[[2]][[i]]$skewdiff
                                  )
      }

      ### Agostino ### ---------------------------------------------------------
      agos <- matrix(NA, length(mod_names), 6)
      rownames(agos) <- mod_names
      colnames(agos) <- rep(c("statistic", "z-value", "p-value"), 2)

      for (i in seq_along(mod_names)) {

        agos[i, ]<- c(unlist(cdda.output[[1]][[i]]$agostino$predictor),
                      unlist(cdda.output[[2]][[i]]$agostino$predictor)
        )
        agos[i, ] <- round(agos[i, ], 3)
      }

      cat(paste("D'Agostino Skewness Tests", "\n"))
      cat(paste("----------------------------------------------------------------------", "\n"))
      cat(paste("            Target Model                 Alternative Model ", "\n"))
      cat(paste("----------------------------------------------------------------------", "\n"))
      print(agos)

      cat("\n")
      cat(paste("Skewness differences"), "\n")
      cat(paste("----------------------------------------------------------------------", "\n"))
      cat(paste("            Target Model       Alternative Model"), "\n")
      cat(paste("----------------------------------------------------------------------", "\n"))
      print(round(momdif.stats.skew, 3))

      cat("---", "\n")
    }

  if (coskew == TRUE){

    if(is.null(cdda.output[[1]][[1]]$cor12diff)) { stop("Co-Skewness differences not found. Specify coskew = TRUE.") }
    hoctests.skew <-  matrix(NA, length(mod_names), 6)

    for(i in seq_along(mod_names)){
      hoctests.skew[i, ] <- c(cdda.output[[1]][[i]]$cor12diff,
                              cdda.output[[2]][[i]]$cor12diff)
      hoctests.skew[i, ] <- round(hoctests.skew[i, ], 3)
    }

    cat(boot_print)

    cat(paste("Co-Skewness differences"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("             Target Model            Alternative Model"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))

    rownames(hoctests.skew) <- mod_names
    colnames(hoctests.skew) <- rep(c("diff", "lower", "upper"), 2)

    print.default(format(hoctests.skew, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")

    if(cdda.output[[1]][[1]]$boot.warning) { cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n", "        Cor^2[3,1] - Cor^2[1,3] should also be computed for the model", varnames[1], "->", varnames[2], "\n") }

    ### Hyvarinen-Smith (co-skewness) ### --------------------------------------
    rhs <- matrix(NA, length(mod_names), 6)

    for (i in seq_along(mod_names)) {
      rhs[i, ] <- c(cdda.output[[1]][[i]]$RHS, cdda.output[[2]][[i]]$RHS)
      rhs[i, ] <- round(rhs[i, ], 3)
    }

    rownames(rhs) <- mod_names
    colnames(rhs) <- rep(c("diff", "lower", "upper"), 2)

    cat(paste("Hyvarinen-Smith Co-Skewness differences", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("            Target Model        Alternative Model ", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    print(rhs)

    cat("---", "\n")
  }

    if (kurt == TRUE){

      ### Anscombe ### ---------------------------------------------------------
      ansc <- matrix(NA, length(mod_names), 6)
      rownames(ansc) <- mod_names
      colnames(ansc) <- rep(c("statistic", "z-value", "p-value"), 2)

      for (i in seq_along(mod_names)) {

        ansc[i, ]<- c(unlist(cdda.output[[1]][[i]]$anscombe$predictor),
                      unlist(cdda.output[[2]][[i]]$anscombe$predictor)
                     )
        ansc[i, ] <- round(ansc[i, ], 3)
       }

      cat(boot_print)

      cat(paste("Anscombe-Glynn Kurtosis Tests", "\n"))
      cat(paste("----------------------------------------------------------------------", "\n"))
      cat(paste("            Target Model                 Alternative Model ", "\n"))
      cat(paste("----------------------------------------------------------------------", "\n"))
      print(ansc)

      cat("---", "\n")

    ### ----------- Boot CI Loop Setup ------------ ###

    momdif.stats.kurt <- matrix(NA, length(mod_names), 6)

    rownames(momdif.stats.kurt) <- mod_names
    colnames(momdif.stats.kurt) <- rep(c("diff", "lower", "upper"), 2)

    ## Boot CI Loop ##
    for (i in seq_along(mod_names)) {

      momdif.stats.kurt[i, ]<- c(cdda.output[[1]][[i]]$kurtdiff, cdda.output[[2]][[i]]$kurtdiff)
     }

    cat(boot_print)

    cat(paste("Excess Kurtosis differences"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("           Target Model            Alternative Model"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))
    print(round(momdif.stats.kurt, 3))

    cat("---", "\n")
  }

    if (cokurt == TRUE){

    ### CoKurt Bootstrap CIS: Corr ### ----------------------------------------

    hoctests.kurt <-  matrix(NA, length(mod_names), 6)

    for(i in seq_along(mod_names)){
      hoctests.kurt[i, ] <- c(cdda.output[[1]][[i]]$cor13diff, cdda.output[[2]][[i]]$cor13diff)
      hoctests.kurt[i, ] <- round(hoctests.kurt[i, ], 3)
    }

    colnames(hoctests.kurt) <- rep(c("diff", "lower", "upper"), 2)
    rownames(hoctests.kurt) <- mod_names

    cat(boot_print)

    cat(paste("Co-Kurtosis differences"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("           Target Model              Alternative Model"), "\n")
    cat(paste("----------------------------------------------------------------------", "\n"))

    print.default(format(hoctests.kurt, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("---", "\n")

    #cat("\n", paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], sep = " "), "\n", "\n") #used orig pred & outcome varz
    if(cdda.output[[1]][[1]]$boot.warning) { cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n", "        Cor^2[3,1] - Cor^2[1,3] should also be computed for the model", varnames[1], "->", varnames[2], "\n") }

    #### Chen-Chan (co-kurtosis) ### -------------------------------------------

    rcc <- matrix(NA, length(mod_names), 6)

    for (i in seq_along(mod_names)) {
      rcc[i, ] <- c(cdda.output[[1]][[i]]$RCC, cdda.output[[2]][[i]]$RCC)
      rcc[i, ] <- round(rcc[i, ], 3)
    }

    rownames(rcc) <- mod_names
    colnames(rcc) <- rep(c("diff", "lower", "upper"), 2)

    cat(boot_print)

    cat(paste("Chen-Chan Co-Kurtosis differences", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("             Target Model         Alternative Model ", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    print(rcc)

    cat("---", "\n")

    ### Hyperbolic Tangent (Rtanh) ### -----------------------------------------

    rtanh <- matrix(NA, length(mod_names), 6)

    for (i in seq_along(mod_names)) {

      rtanh[i, ] <- c(cdda.output[[1]][[i]]$Rtanh, cdda.output[[2]][[i]]$Rtanh)
      rtanh[i, ] <- round(rtanh[i, ], 3)
    }

    rownames(rtanh) <- mod_names
    colnames(rtanh) <- rep(c("diff", "lower", "upper"), 2)

    cat(boot_print)

    cat(paste("Hyperbolic Tangent", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    cat(paste("            Target Model          Alternative Model ", "\n"))
    cat(paste("----------------------------------------------------------------------", "\n"))
    print(rtanh)

    cat("---", "\n")
    }
}



