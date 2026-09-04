#' @title Plots of \code{cdda.vardist} Class Objects
#' @description \code{plot} returns graphs of \code{cdda.vardist} results.
#'
#' @name plot.cdda.vardist
#'
#' @param x      An object of class \code{cdda.vardist} when using \code{print} or \code{plot}.
#' @param stat   A character indicating the statistic to be plotted, default is
#'               \code{"rhs"}, with options \code{c("coskew", "cokurt", "rhs", "rcc", "rtanh")}.
#' @param ylim   A numeric vector of length 2 indicating the y-axis limits. If \code{NULL}, the function will set the limits automatically.
#' @param ...    Additional arguments to be passed to the function.
#'
#' @examples plot(result, stat = "rtanh", ylim = c(-0.05, 0.05))
#'
#' @export
#' @rdname cdda.vardist
#' @method plot cdda.vardist
plot.cdda.vardist <- function(x,
                              stat = NULL,
                              ylim =  NULL,
                              ...
                              ){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(is.null(stat)){
    stop("stat argument must be specified. as 'rhs', 'coskew', 'cokurt', 'rcc', or 'rtanh'")
  }
  else if(stat != "rhs" & stat != "cokurt" & stat != "rcc" &
          stat != "rtanh" & stat != "coskew"){
    stop("stat must be one of 'rhs', 'coskew', 'cokurt', 'rcc', or 'rtanh'.")
  }

  obj <- x

  ## Plot Label and Assignment ##
  mod.vals <- obj[[4]][["mod_data"]]
  mod.levels <- x.axis.labels <- obj[[4]][["mod_levels"]] #modval levels for pick-a-point
  y.title <- obj[[3]][[1]] #Test statistic CI header
  ci.level <- as.numeric(obj[[1]][[1]]$boot.args[2]) * 100

  ### Limits for X-axis ####

  if(is.numeric(mod.levels)){
    x.range <- c(min(mod.levels), max(mod.levels))
    plot.axis <- mod.levels
  }
  else {x.range <- c(1, length(mod.levels))
  plot.axis <- seq_along(mod.levels)
  }

  out <- matrix(NA, length(plot.axis), 6)

  if (stat == "rhs"){
    y.title <- paste0("Hyvarinen-Smith Co-Skewness Differences (", ci.level, "% CI)")
    if(is.null(obj[[1]][[1]]$RHS)){
      stop( "Hyvarinen-Smith Co-Skewness not found. Specify coskew = TRUE." )
      y.title <- "RHS"
    }

    for (i in seq_along(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$RHS, obj[[2]][[i]]$RHS)
    }
  }

  if (stat == "cokurt"){

    y.title <- paste0("Co-Kurtosis Differences (", ci.level, "% CI)")
    if(is.null(obj[[1]][[1]]$cor13diff)){
      stop( "mitests not found. Specify cokurt = TRUE." )
      y.title <- "Co-Kurtosis"
    }

    for(i in seq_along(plot.axis)){
      out[i, ] <- c(obj[[1]][[i]]$cor13diff, obj[[2]][[i]]$cor13diff)
    }
  }

  if (stat == "coskew"){

    y.title <- paste0("Co-Skewness Differences (", ci.level, "% CI)")


    if(is.null(obj[[1]][[1]]$cor12diff)) {
      stop("Co-Skewness differences not found. Specify coskew = TRUE.")
      }
    hoctests.skew <-  matrix(NA, length(mod.levels), 6)

    for(i in seq_along(plot.axis)){
      out[i, ] <- c(obj[[1]][[i]]$cor12diff, obj[[2]][[i]]$cor12diff)
    }
  }

  if (stat == "rcc"){
    y.title <- paste0("Chen-Chan Co-Kurtosis Differences (", ci.level, "% CI)")
    if(is.null(obj[[1]][[1]]$RCC)){
      stop( "Chen-Chan Co-Kurtosis Differences not found. Specify cokurt = TRUE." )
    }

    for (i in seq_along(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$RCC, obj[[2]][[i]]$RCC)
    }
  }

  if (stat == "rtanh"){
    y.title <- paste0("Hyperbolic Tangent Differences (", ci.level, "% CI)")
    if(is.null(obj[[1]][[1]]$Rtanh)){
      stop( "Hyperbolic Tangent Differences not found. Specify cokurt = TRUE." )
    }

    for (i in seq_along(plot.axis)) {
      out[i, ] <- c(obj[[1]][[i]]$Rtanh, obj[[2]][[i]]$Rtanh)
    }
  }

  ### Vardist Plotting ################################################

  ### Matrix Creation ###

  out.tar <- data.frame(condition = plot.axis, out.mean = out[,1],
                        out.low = out[,2], out.upp = out[,3])
  out.alt <- data.frame(condition = plot.axis, out.mean = out[,4],
                        out.low = out[,5], out.upp = out[,6])

  ### Set the y-axis ###

  if(is.null(ylim)) { y.range <- c(min(c(out.tar$out.low, out.alt$out.low)) - 0.1,
                                   max(c(out.tar$out.upp, out.alt$out.upp) + 0.1) )
  } else{y.range <- ylim}

  ### Start two plots ###

  par(mfrow = c(1, 2))

  plot(plot.axis, out.tar$out.mean, type = "n",
       ylim = y.range, xlim = x.range, xaxt = "n",
       xlab = "Moderator Values", ylab = y.title, main = "Target Model")
  axis(1, at = seq_along(out.tar$condition), labels = x.axis.labels)
  polygon(x = c(out.tar$condition, rev(out.tar$condition)),
          y = c(out.tar$out.low, rev(out.tar$out.upp)),
          border = FALSE, col = "lightgrey")
  points(out.tar$condition, out.tar$out.mean, type = "l")
  abline(h = 0, lty = "dashed")


  plot(plot.axis, out.alt$out.mean, type = "n",
       ylim = y.range, xlim = x.range, xaxt = "n",
       xlab = "Moderator Values", ylab = y.title, main = "Alternative Model")
  axis(1, at = seq_along(out.tar$condition), labels = x.axis.labels)
  polygon(x = c(out.alt$condition, rev(out.alt$condition)),
          y = c(out.alt$out.low, rev(out.alt$out.upp)),
          border = FALSE, col = "lightgrey")
  points(out.alt$condition, out.alt$out.mean, type = "l")
  abline(h = 0, lty = "dashed")
}

par(mfrow = c(1, 1))


