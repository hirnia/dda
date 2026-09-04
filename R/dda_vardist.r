#' @title Direction Dependence Analysis: Variable Distributions
#'
#' @description \code{dda.vardist} evaluates patterns of asymmetry of
#' variable distributions for causally competing models
#' (\code{y ~ x} vs. \code{x ~ y}).
#' \code{print} returns DDA test statistics associated with
#' \code{dda.vardist} objects.
#'
#' @param formula Symbolic formula of the model to be tested or a
#'   \code{lm} object.
#' @param pred Variable name of the predictor which serves as the outcome in
#'   the alternative model.
#' @param data An optional data frame containing the variables in the model
#'   (by default variables are taken from the environment which
#'   \code{dda.vardist} is called from).
#' @param B Number of bootstrap samples.
#' @param boot.type A character indicating the type of bootstrap confidence
#'   intervals. Must be one of the two specifications
#'   \code{c("perc", "bca")}. \code{boot.type = "perc"} is the default.
#' @param conf.level Confidence level for bootstrap confidence intervals.
#'
#' @return An object of class \code{dda.vardist} containing the results of
#'   direction dependence tests of variable distributions.
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction
#'   Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK:
#'   Cambridge University Press.
#'
#' @seealso \code{\link{cdda.vardist}} for a conditional version.
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rchisq(n, df = 4) - 4
#' e <- rchisq(n, df = 3) - 3
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' result <- dda.vardist(y ~ x, pred = "x", data = d, B = 50)
#' # Note: Only 50 bootstrap samples are used here to lower computation time
#'
#' \dontrun{
#' ## --- Realistic settings; run time is substantial
#'
#' result <- dda.vardist(y ~ x, pred = "x", data = d, B = 500,
#'   boot.type = "bca", conf.level = 0.95)
#'
#' print(result)
#' }
#'
#' @export
#' @rdname dda.vardist
dda.vardist <- function(formula,
                        pred = NULL,
                        data = list(),
                        B = 200,
                        boot.type = "perc",
                        conf.level = 0.95
                        ){

   ### --- helper functions for bootstrap CIs

   mysd <- function(x){sqrt(sum((x-mean(x))^2)/length(x))}

   cor.ij <- function(x,y, i=1, j=1){
       n <- length(x)
       mx <- mean(x)
       my <- mean(y)
       Cov <- sum((x - mx)^i * (y - my)^j)/n
       Cov/(mysd(x)^i * mysd(y)^j)
    }

    boot.diff <- function(dat, g){
              dat <- dat[g, ]
              x <- dat[, 1]  # "purified" predictor
              y <- dat[, 2]  # "purified" outcome

              x <- as.vector(scale(x))
			        y <- as.vector(scale(y))

              skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
              kurt.diff <- (moments::kurtosis(x)-3)^2 - (moments::kurtosis(y)-3)^2
              cor12.diff <- (cor.ij(x, y, i = 2, j = 1)^2) - (cor.ij(x, y, i = 1, j = 2)^2)
			        cor13.diff <- ((cor.ij(x, y, i = 3, j = 1)^2) - (cor.ij(x, y, i = 1, j = 3)^2)) * sign(moments::kurtosis(x)-3)

              Rtanh <- cor(x, y) * mean(x * tanh(y) - tanh(x) * y)

			        Cxy <- mean(x^3 * y) - 3*cor(x,y)*var(x)
			        Cyx <- mean(x * y^3) - 3*cor(x,y)*var(y)
              RCC <- (Cxy + Cyx) * (Cxy - Cyx)

              xx <- sign(moments::skewness(x)) * x
			        yy <- sign(moments::skewness(y)) * y
              RHS <- cor(xx, yy) * mean( (xx^2 * yy) - (xx * yy^2) )

              result <- c(skew.diff, kurt.diff, cor12.diff, cor13.diff, RHS, RCC, Rtanh)
              names(result) <- c("skew.diff", "kurt.diff", "cor21.diff", "cor13.diff", "RHS", "RCC", "Rtanh")
              return(result)
    }

  if(is.null(pred)) stop( "Tentative predictor is missing." )
	if(B <= 0) stop( "Number of resamples 'B' must be positive." )
	if(conf.level < 0 || conf.level > 1) stop("'conf.level' must be between 0 and 1")
	if( !boot.type %in% c("bca", "perc") ) stop( "Unknown argument in boot.type." )

	### --- prepare outcome, predictor, and model matrix for covariates

    if (!inherits(formula, "formula") ) {
           X <- if (is.matrix(formula$x) ) formula$x
           else model.matrix(terms(formula), model.frame(formula) )
           y <- if (is.vector(formula$y) ) formula$y
           else model.response(model.frame(formula))

           delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
           if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[, delete.pred]   # tentative predictor
		   X <- X[ , -delete.pred] # model matrix with covariates
		   if ( !is.matrix(X) ) X <- as.matrix(X)
	}
	else {
       mf <- model.frame(formula, data = data)
		   y <- model.response(mf)   # tentative outcome
		   X <- model.matrix(formula, data = data)

		   delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
           if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

		   x <- X[, delete.pred]   # tentative predictor
		   X <- X[ , -delete.pred] # model matrix with covariates
		   if (!is.matrix(X)) X <- as.matrix(X)
	}


  ry <- lm.fit(X, y)$residuals
	rx <- lm.fit(X, x)$residuals

	ry <- as.vector(scale(ry))
	rx <- as.vector(scale(rx))




    dat <- data.frame(predictor = rx, outcome = ry)

	### --- run separate normality tests

	agostino.out <- apply(dat, 2, moments::agostino.test)
	agostino.out <- lapply(agostino.out, unclass)
	agostino.out$predictor[3:5] <- NULL
  agostino.out$outcome[3:5] <- NULL

	anscombe.out <- apply(dat, 2, moments::anscombe.test)
	anscombe.out <- lapply(anscombe.out, unclass)
	anscombe.out$predictor[3:5] <- NULL
  anscombe.out$outcome[3:5] <- NULL

	output <- list(agostino.out, anscombe.out)
	names(output) <- c("agostino", "anscombe")

  output$anscombe$outcome$statistic[1] <- output$anscombe$outcome$statistic[1] - 3     # change kurtosis to ex-kurtosis
	output$anscombe$predictor$statistic[1] <- output$anscombe$predictor$statistic[1] - 3 # change kurtosis to ex-kurtosis

    ### --- run bootstrap confidence intervals

   if(B > 0){
     boot.res <- boot::boot(dat, boot.diff, R = B)

    #if( boot.type == "bca" && any( is.na( empinf( boot.res ) ) ) ) stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'")
     if( boot.type == "bca"){
       boot.out <- suppressWarnings( try(
         lapply(as.list(1:7), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type="bca", t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res), silent = TRUE))
       if (inherits(boot.out, "try-error")) { stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'") }
     }

     if( boot.type == "perc"){
       suppressWarnings(boot.out <- lapply(as.list(1:7), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type="perc", t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res))

     }

			 names(boot.out) <- c("skew.diff", "kurt.diff", "cor12.diff", "cor13.diff", "RHS", "RCC", "Rtanh")

       ci.skewdiff  <- unclass(boot.out$skew.diff)[[4]][4:5] ; names(ci.skewdiff) <- c("lower", "upper")
			 ci.kurtdiff  <- unclass(boot.out$kurt.diff)[[4]][4:5] ; names(ci.kurtdiff) <- c("lower", "upper")
			 ci.cor12diff <- unclass(boot.out$cor12.diff)[[4]][4:5] ; names(ci.cor12diff) <- c("lower", "upper")
			 ci.cor13diff <- unclass(boot.out$cor13.diff)[[4]][4:5] ; names(ci.cor13diff) <- c("lower", "upper")
			 ci.RHS       <- unclass(boot.out$RHS)[[4]][4:5] ; names(ci.RHS) <- c("lower", "upper")
			 ci.RCC       <- unclass(boot.out$RCC)[[4]][4:5] ; names(ci.RCC) <- c("lower", "upper")
			 ci.Rtanh     <- unclass(boot.out$Rtanh)[[4]][4:5] ; names(ci.Rtanh) <- c("lower", "upper")

			 output <- c(output,
			       list(skewdiff = c(boot.res$t0[1], ci.skewdiff)),
						 list(kurtdiff = c(boot.res$t0[2], ci.kurtdiff)),
						 list(cor12diff = c(boot.res$t0[3], ci.cor12diff)),
						 list(cor13diff = c(boot.res$t0[4], ci.cor13diff)),
						 list(RHS = c(boot.res$t0[5], ci.RHS)),
						 list(RCC = c(boot.res$t0[6], ci.RCC)),
						 list(Rtanh = c(boot.res$t0[7], ci.Rtanh)),
						 list(boot.args = c(boot.type, conf.level, B)),
						 list(boot.warning = FALSE)
						)
			if(sign(output$anscombe$predictor$statistic[1]) != sign(output$anscombe$outcome$statistic[1])) { output$boot.warning <- TRUE }
	}

	response.name <- all.vars(formula(formula))[1]  # get name of response variable
	output <- c(output, list(var.names = c(response.name, pred)))

	call_info <- list(
	  "function_call" = match.call(),
	  "function_name" = "dda.vardist",  # or deparse(substitute(sys.function()))
	  "all_args" = as.list(match.call())[-1],
	  "formula" = formula,
	  "data_name" = deparse(substitute(data))#, RM original_data from all ddas
	 # "original_data" = if(missing(data) || is.null(data)) NULL else data
	)

	output <- c(output,
	            list(call_info = call_info)
	           )

	class(output) <- "dda.vardist"
	return(output)

}

#' @title Print Method for \code{dda.vardist} Objects
#'
#' @param x An object of class \code{dda.vardist} when using \code{print}.
#' @param ... Additional arguments to be passed to the function.
#'
#' @examples
#' print(result)
#'
#' @export
#' @rdname dda.vardist
#' @method print dda.vardist
print.dda.vardist <- function(x, ...){
   varnames <- x$var.names

	 cat("\n")
     cat("DIRECTION DEPENDENCE ANALYSIS: Variable Distributions", "\n", "\n")
     cat("Skewness and kurtosis tests:", "\n")

	     sigtests <- rbind( c(x[[1]]$outcome$statistic, x[[1]]$outcome$p.value, x[[1]]$predictor$statistic, x[[1]]$predictor$p.value),
	                        c(x[[2]]$outcome$statistic, x[[2]]$outcome$p.value, x[[2]]$predictor$statistic, x[[2]]$predictor$p.value)
	                      )
         sigtests <- round(sigtests, 4)
         rownames(sigtests) <- c("Skewness", "Kurtosis")
		     colnames(sigtests) <- c(varnames[1], " z-value", " Pr(>|z|)", varnames[2], " z-value", " Pr(>|z|)")
         print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

	 if(!is.null(x$boot.args)){
	     ci.level <- as.numeric(x$boot.args[2]) * 100
		   cat("\n")

	       # Print Skewness and Kurtosis based measures

		   if(x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for higher moment differences:", "\n", sep = "")
       if(x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for higher moment differences:", "\n", sep = "")

	       citests <- rbind(x$skewdiff, x$kurtdiff)
		     citests <- round(citests, 4)
	       rownames(citests) <- c("Skewness", "Kurtosis")
	       colnames(citests) <- c("diff", "lower", "upper")
	       print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
		     cat("\n")

           # Print Co-Skewness and Co-Kurtosis based measures

   	       if(x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for differences in higher-order correlations:", "\n", sep = "")
           if(x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for differences in higher-order correlations:", "\n", sep = "")

		   hoctests <- rbind(x$cor12diff, x$cor13diff)
		   hoctests <- round(hoctests, 4)
		   rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]" )
		   colnames(hoctests) <- c("estimate", "lower", "upper")
		   print.default(format( hoctests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
		   cat("\n")

           # Print LR approximative measures

   	       if(x$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")
           if(x$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")

		   LRtests <- rbind(x$RHS, x$Rtanh, x$RCC )
		   LRtests <- round(LRtests, 4)
		   rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (tanh)", "Chen-Chan (co-kurtosis)")
		   colnames(LRtests) <- c("estimate", "lower", "upper")
		   print.default(format( LRtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

		     cat("\n")
	       cat(paste("Number of resamples:", x$boot.args[3]))
	       cat("\n")
	       cat("---")
	       cat("\n")
	       cat(paste("Note: (Cor^2[i,j] - Cor^2[j,i]) > 0 suggests the model", varnames[2], "->", varnames[1], sep = " "))
	       cat("\n")
		   if(x$boot.warning) { cat("Warning: Excess-kurtosis values of", varnames[2], "and", varnames[1], "have unequal signs", "\n", "        Cor^2[3,1] - Cor^2[1,3] should also be computed for the model", varnames[1], "->", varnames[2], "\n") }
      }
}


