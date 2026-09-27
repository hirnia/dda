#' @title Direction Dependence Analysis: Residual Distributions
#'
#' @description \code{dda.resdist} evaluates patterns of asymmetry of error
#' distributions of causally competing models (\code{y ~ x} vs. \code{x ~ y}).
#' \code{print} returns DDA test statistics associated with \code{dda.resdist}
#' objects.
#'
#' @param formula Symbolic formula of the target model to be tested or an
#'   \code{lm} object.
#' @param pred Variable name of the predictor which serves as the outcome in
#'   the alternative model.
#' @param data An optional data frame containing the variables in the model
#'   (by default variables are taken from the environment which
#'   \code{dda.resdist} is called from).
#' @param B Number of bootstrap samples.
#' @param boot.type A vector of character strings representing the type of
#'   bootstrap confidence intervals. Must be one of the two values
#'   \code{c("perc", "bca")}; \code{boot.type = "perc"} is the default.
#' @param prob.trans A logical value indicating whether a probability integral
#'   transformation should be performed prior computation of skewness and
#'   kurtosis tests.
#' @param conf.level Confidence level for bootstrap confidence intervals.
#' @param robust A logical value indicating whether Siegel's (1982) repeated
#'   median estimator should be used for model estimation. If
#'   \code{robust = TRUE} repeated median estimation is applied for the
#'   causally competing models, otherwise ordinary least squares (OLS)
#'   estimation is used.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return An object of class \code{dda.resdist} containing the results of
#'   direction dependence tests of error distributions.
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction
#'   Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK:
#'   Cambridge University Press.
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rchisq(n, df = 4) - 4
#' e <- rnorm(n, sd = sqrt(6))
#' y <- 0.5 * x + e
#' d <- data.frame(x, y)
#'
#' result <- dda.resdist(y ~ x, pred = "x", data = d,
#'             B = 50, conf.level = 0.90, prob.trans = FALSE, robust = FALSE)
#'
#' print(result)
#'
#' ## --- prob.trans = TRUE maps the predictor onto the marginal of the outcome,
#' ## --- so it is informative when the outcome is itself non-normal.
#'
#' d.skew <- data.frame(x, y = 0.5 * x + (rchisq(n, df = 3) - 3))
#'
#' dda.resdist(y ~ x, pred = "x", data = d.skew,
#'             B = 50, conf.level = 0.90, prob.trans = TRUE)
#'
#' \dontrun{
#'
#' result <- dda.resdist(y ~ x, pred = "x", data = d, B = 500,
#'   boot.type = "bca", conf.level = 0.95, prob.trans = FALSE)
#'
#' print(result)
#'
#' ## --- Siegel (1982) repeated median estimation
#'
#' dda.resdist(y ~ x, pred = "x", data = d, B = 500, robust = TRUE)
#' }
#'
#' @export
#' @rdname dda.resdist
dda.resdist <- function(formula,
                        pred = NULL,
                        data = list(),
                        B = 200,
                        boot.type = "perc",
                        prob.trans = FALSE,
                        conf.level = 0.95,
                        robust = FALSE
){

   ### --- helper function for bootstrap CIs

   mysd <- function(x){sqrt(sum((x-mean(x))^2)/length(x))}
   cor.ij <- function(x,y, i=1, j=1){
       n <- length(x)
       mx <- mean(x)
       my <- mean(y)
       Cov <- sum((x - mx)^i * (y - my)^j)/n
       Cov/(mysd(x)^i * mysd(y)^j)
    }


	prob.int <- function(x, y) {

		fval <- ecdf(y)(y)
		xval <- y
		quant_y <- approxfun(fval, xval, ties = mean) # 0.2.0 ties = mean for tied values in resampled data
		cum_x <- ecdf(x)
		transform <- function(x) quant_y(cum_x(x))
		x <- transform(x)
		x[is.na(x)] <- min(y)
		list(x = x, y = y)
	}


    boot.diff <- function(dat, g, prob.trans, robust){
              dat <- dat[g, ]

			  x <- dat[, 1]  # "purified" alternative error
              y <- dat[, 2]  # "purified" target error

              x <- as.vector(scale(x))
			  y <- as.vector(scale(y))

              if( isTRUE(prob.trans) ){

				   xboot <- dat[, 3] # tentative (cov-adjusted) predictor
                   yboot <- dat[, 4] # tentative (cov-adjusted) outcome

				   trans.boot <- prob.int(xboot, yboot)

				   xboot.trans <- trans.boot$x
                   yboot.trans <- trans.boot$y

		           if( isTRUE(robust) ){ # 0.2.0 robust estimation
		             tar.boot.trans <- RobustLinearReg::siegel_regression(yboot.trans ~ xboot.trans)
		             alt.boot.trans <- RobustLinearReg::siegel_regression(xboot.trans ~ yboot.trans)
		           } else {
		             tar.boot.trans <- lm(yboot.trans ~ xboot.trans)
		             alt.boot.trans <- lm(xboot.trans ~ yboot.trans)
		           }

				   xtrans <- as.vector(scale(resid(alt.boot.trans)))
				   ytrans <- as.vector(scale(resid(tar.boot.trans)))

				   # 0.2.0 target minus alternative under prob.trans, so that > 0 speaks for the target model
			       skew.diff <- (moments::skewness(ytrans)^2) - (moments::skewness(xtrans)^2)
                   kurt.diff <- (moments::kurtosis(ytrans)-3)^2 - (moments::kurtosis(xtrans)-3)^2

			  } else{

			       skew.diff <- (moments::skewness(x)^2) - (moments::skewness(y)^2)
                   kurt.diff <- (moments::kurtosis(x)-3)^2 - (moments::kurtosis(y)-3)^2
              }


			  cor12.diff <- (cor.ij(y, x, i = 2, j = 1)^2) - (cor.ij(y, x, i = 1, j = 2)^2)  # rescaled to ensure > 0 under target model
			  cor13.diff <- ((cor.ij(y, x, i = 3, j = 1)^2) - (cor.ij(y, x, i = 1, j = 3)^2)) * sign(moments::kurtosis(y)-3)  # rescaled to ensure > 0 under target model

			  xx <- sign(moments::skewness(x)) * x
			  yy <- sign(moments::skewness(y)) * y
			  RHS3 <- cor(xx, yy) * mean( (yy^2 * xx) - (yy * xx^2) ) # rescaled to ensure > 0 under target model

              RHS4 <- cor(x, y) * mean( (y^3 * x) - (y * x^3) ) * sign(moments::kurtosis(y)-3)

			  C1 <- mean(y^3 * x) - 3*cor(x,y)*var(y) # with x = alternative error and y = target error
			  C2 <- mean(y * x^3) - 3*cor(x,y)*var(x)
              RCC <- (C1 + C2) * (C1 - C2)	# scaled to ensure that > 0 under target model

              result <- c(skew.diff, kurt.diff, cor12.diff, cor13.diff, RHS3, RCC, RHS4)
              names(result) <- c("skew.diff", "kurt.diff", "cor21.diff", "cor13.diff", "RHS3", "RCC", "RHS4")
              return(result)
    }

    skew.diff.test <- function(x, y){

           agostino.zvalue <- function(x){
                     n  <- length(x)
                     s3 <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
                     y <- s3 * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
                     b2 <- 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)/((n - 2) * (n + 5) * (n + 7) * (n + 9))
                     w <- -1 + sqrt(2 * (b2 - 1))
                     d <- 1/sqrt(log(sqrt(w)))
                     a <- sqrt(2/(w - 1))
                     z <- d * log(y/a + sqrt((y/a)^2 + 1))
                     return(z)
            }

		zval <- (agostino.zvalue(x) - agostino.zvalue(y))/sqrt(2 - 2*cor(x,y)^3)
		pval <- (1 - pnorm(abs(zval))) * 2 # two-sided pvalue
		return(list(z.value = zval, p.value = pval))
	}

	kurt.diff.test <- function(x, y){

           anscombe.zvalue <- function(x){
                     n   <- length(x)
                     b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
                     eb2 <- 3 * (n - 1)/(n + 1)
                     vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
                     m3  <- (6 * (n^2 - 5 * n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
                     a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
                     xx  <- (b - eb2)/sqrt(vb2)

					 # modification when cube rooted value is negative: use sign(x) * abs(x)^(1/3) instead of (-x)^(1/3)
					 cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
					 z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))

                     return(z)
            }

        zval <- (anscombe.zvalue(x) - anscombe.zvalue(y))/sqrt(2 - 2*cor(x,y)^4)
		pval <- (1 - pnorm(abs(zval))) * 2 # two-sided pvalue
		return(list(z.value = zval, p.value = pval))
	}


	myanscombe.test <- function(x, alternative = c("two.sided", "less", "greater")){

	# modified anscombe.test() function to handle cube roots of negative numbers

	          n   <- length(x)
			  DNAME <- deparse(substitute(x))
              x <- sort(x[complete.cases(x)])
              n <- length(x)
              s <- match.arg(alternative)
              alter <- switch(s, two.sided = 0, less = 1, greater = 2)

              b   <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
              eb2 <- 3 * (n - 1)/(n + 1)
              vb2 <- 24 * n * (n - 2) * (n - 3)/((n + 1)^2 * (n + 3) * (n + 5))
              m3  <- (6 * (n^2 - 5 * n + 2)/((n + 7) * (n + 9))) * sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
              a   <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2))
              xx  <- (b - eb2)/sqrt(vb2)

			 # modification when cube rooted value is negative: use sign(x) * abs(x)^(1/3) instead of (-x)^(1/3)
				cr <- sign((1 - 2/a)/(1 + xx * sqrt(2/(a - 4)))) * abs((1 - 2/a)/(1 + xx * sqrt(2/(a - 4))))^(1/3)
				z   <- (1 - 2/(9 * a) - cr ) / sqrt(2/(9 * a))

		pval <- pnorm(z, lower.tail = FALSE)
     if (alter == 0) {
        pval <- 2 * pval
        if (pval > 1)
            pval <- 2 - pval
        alt <- "kurtosis is not equal to 3"
      }
     else if (alter == 1) {
        alt <- "kurtosis is greater than 3"
     }
     else {
        pval <- 1 - pval
        alt <- "kurtosis is lower than 3"
     }
      RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval,
        alternative = alt, method = "Anscombe-Glynn kurtosis test",
        data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
	}


    if(is.null(pred)) stop( "Tentative predictor is missing." )
	if(B < 0) stop( "Number of resamples 'B' must be positive." )
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

	### --- estimate competing models

	ry <- as.vector(scale(ry))
	rx <- as.vector(scale(rx))

	if( isTRUE(robust) ){ # 0.2.0 robust estimation
	  tar <- RobustLinearReg::siegel_regression(ry ~ rx)
	  alt <- RobustLinearReg::siegel_regression(rx ~ ry)
	} else {
	  tar <- lm(ry ~ rx)
	  alt <- lm(rx ~ ry)
	}

	dat <- data.frame(alternative = as.vector(scale(resid(alt))),
	                  target = as.vector(scale(resid(tar))),
                      pred.adj = rx,
					  out.adj = ry
					  )  # note: residuals are standardized prior DDA

    if( isTRUE(prob.trans) ){  # add prob transformed if needed

		trans <- prob.int(rx, ry)
		rx.trans <- trans$x
        ry.trans <- trans$y

		if( isTRUE(robust) ){ # 0.2.0 robust estimation
		  tar.trans <- RobustLinearReg::siegel_regression(ry.trans ~ rx.trans)
		  alt.trans <- RobustLinearReg::siegel_regression(rx.trans ~ ry.trans)
		} else {
		  tar.trans <- lm(ry.trans ~ rx.trans)
		  alt.trans <- lm(rx.trans ~ ry.trans)
		}

		dat$alternative.trans <- as.vector(scale(resid(alt.trans)))
		dat$target.trans <- as.vector(scale(resid(tar.trans)))

	}

	### --- run separate normality tests


	if( isTRUE(prob.trans) ){

	  agostino.out <- apply(dat[c("alternative.trans", "target.trans")], 2, moments::agostino.test)
	  names(agostino.out) <- c("alternative", "target")
	  agostino.out <- lapply(agostino.out, unclass)

	  anscombe.out <- apply(dat[c("alternative.trans", "target.trans")], 2, myanscombe.test)
	  names(anscombe.out) <- c("alternative", "target")
	  anscombe.out <- lapply(anscombe.out, unclass)

	} else {

	  agostino.out <- apply(dat[c("alternative", "target")], 2, moments::agostino.test)
	  agostino.out <- lapply(agostino.out, unclass)

	  anscombe.out <- apply(dat[c("alternative", "target")], 2, myanscombe.test)
	  anscombe.out <- lapply(anscombe.out, unclass)

	}
	  agostino.out$alternative[3:5] <- NULL
      agostino.out$target[3:5] <- NULL
	  anscombe.out$alternative[3:5] <- NULL
      anscombe.out$target[3:5] <- NULL

	output <- list(agostino.out, anscombe.out)
	names(output) <- c("agostino", "anscombe")

	output$anscombe$target$statistic[1] <- output$anscombe$target$statistic[1] - 3 # change kurtosis to ex-kurtosis
	output$anscombe$alternative$statistic[1] <- output$anscombe$alternative$statistic[1] - 3 # change kurtosis to ex-kurtosis

	### --- run asymptotic difference tests

	if(isTRUE(prob.trans)){

	 # 0.2.0 target minus alternative under prob.trans, so that > 0 speaks for the target model
	 output <- c(output,
	            list(skewdiff = unlist(skew.diff.test(dat$target.trans, dat$alternative.trans))),
				list(kurtdiff = unlist(kurt.diff.test(dat$target.trans, dat$alternative.trans)))
	  )

	 output$skewdiff <- c( moments::skewness(dat$target.trans)^2 - moments::skewness(dat$alternative.trans)^2, output$skewdiff)  # add point estimtes of skew and kurt differences to vector
	 output$kurtdiff <- c((moments::kurtosis(dat$target.trans)-3)^2 - (moments::kurtosis(dat$alternative.trans)-3)^2, output$kurtdiff)


	} else {

	  output <- c(output,
	            list(skewdiff = unlist(skew.diff.test(dat$alternative, dat$target))),
				list(kurtdiff = unlist(kurt.diff.test(dat$alternative, dat$target)))
	  )

	output$skewdiff <- c( moments::skewness(dat$alternative)^2 - moments::skewness(dat$target)^2, output$skewdiff)  # add point estimtes of skew and kurt differences to vector
	output$kurtdiff <- c((moments::kurtosis(dat$alternative)-3)^2 - (moments::kurtosis(dat$target)-3)^2, output$kurtdiff)

	}

    ### --- run bootstrap confidence intervals

   if(B > 0){
             suppressWarnings(boot.res <- boot::boot(dat, boot.diff, R = B, prob.trans = prob.trans, robust = robust))

			 if( boot.type == "bca"){
			 boot.out <- suppressWarnings( try(
			         lapply(as.list(1:7), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type="bca", t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res), silent = TRUE))
			 if (inherits(boot.out, "try-error")) { stop("Acceleration constant cannot be calculated. Increase the number of resamples or use boot.type = 'perc'") }
			 }

			 if( boot.type == "perc"){
			 suppressWarnings(boot.out <- lapply(as.list(1:7), function(i, boot.res) boot::boot.ci(boot.res, conf=conf.level, type="perc", t0=boot.res$t0[i], t=boot.res$t[,i]), boot.res=boot.res))
			 }

			 names(boot.out) <- c("skew.diff", "kurt.diff", "cor12.diff", "cor13.diff", "RHS3", "RCC", "RHS4")

             ci.skewdiff  <- unclass(boot.out$skew.diff)[[4]][4:5] ; names(ci.skewdiff) <- c("lower", "upper")
			 ci.kurtdiff  <- unclass(boot.out$kurt.diff)[[4]][4:5] ; names(ci.kurtdiff) <- c("lower", "upper")
			 ci.cor12diff <- unclass(boot.out$cor12.diff)[[4]][4:5] ; names(ci.cor12diff) <- c("lower", "upper")
			 ci.cor13diff <- unclass(boot.out$cor13.diff)[[4]][4:5] ; names(ci.cor13diff) <- c("lower", "upper")
			 ci.RHS3      <- unclass(boot.out$RHS3)[[4]][4:5] ; names(ci.RHS3) <- c("lower", "upper")
			 ci.RCC       <- unclass(boot.out$RCC)[[4]][4:5] ; names(ci.RCC) <- c("lower", "upper")
			 ci.RHS4      <- unclass(boot.out$RHS4)[[4]][4:5] ; names(ci.RHS4) <- c("lower", "upper")

			 output$skewdiff <- c(output$skewdiff, ci.skewdiff)
			 output$kurtdiff <-	c(output$kurtdiff, ci.kurtdiff)

			 output <- c(output,
						 list(cor12diff = c(boot.res$t0[3], ci.cor12diff)),
						 list(cor13diff = c(boot.res$t0[4], ci.cor13diff)),
						 list(RHS3 = c(boot.res$t0[5], ci.RHS3)),
						 list(RCC = c(boot.res$t0[6], ci.RCC)),
						 list(RHS4 = c(boot.res$t0[7], ci.RHS4)),
						 list(boot.args = c(boot.type, conf.level, B)),
						 list(boot.warning = FALSE)
						)
			if(sign(output$anscombe$alternative$statistic[1]) != sign(output$anscombe$target$statistic[1])) { output$boot.warning <- TRUE }
	}

	response.name <- all.vars(formula(formula))[1]  # get name of response variable
	output <- c(output, list(var.names = c(response.name, pred), probtrans = prob.trans))

	output$call_info <- list(function_call = match.call(), formula = formula) # 0.2.0 call information for dda.bagging

	class(output) <- "dda.resdist"
	return(output)

}

#' @title Print Method for \code{dda.resdist} Objects
#'
#' @param x An object of class \code{dda.resdist} when using \code{print}.
#' @param ... Additional arguments to be passed to the method.
#'
#' @examples
#' \dontrun{
#' print(result)
#' }
#'
#' @export
#' @rdname dda.resdist
#' @method print dda.resdist
print.dda.resdist <- function(x, ...){
  object <- x

  varnames <- object$var.names
  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n")
  cat("Skewness and kurtosis tests:", "\n")

  sigtests <- rbind( c(object[[1]]$target$statistic, object[[1]]$target$p.value, object[[1]]$alternative$statistic, object[[1]]$alternative$p.value),
                     c(object[[2]]$target$statistic, object[[2]]$target$p.value, object[[2]]$alternative$statistic, object[[2]]$alternative$p.value)
  )
  sigtests <- round(sigtests, 4)
  rownames(sigtests) <- c("Skewness", "Kurtosis")
  colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
  print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)

  if(is.null(object$boot.args)){
    cat("\n")
    cat("Skewness and kurtosis difference tests:", "\n")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")
  }

  if(!is.null(object$boot.args)){
    ci.level <- as.numeric(object$boot.args[2]) * 100
    cat("\n")

    # Print skewness and kurtosis results

    if(object$boot.args[1] == "bca")  cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "BCa bootstrap CIs:", "\n", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat("Skewness and kurtosis difference tests and ", ci.level, "% ", "Percentile bootstrap CIs:", "\n", "\n", sep = "")

    citests <- rbind(object$skewdiff, object$kurtdiff)
    citests <- round(citests, 4)
    rownames(citests) <- c("Skewness", "Kurtosis")
    colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
    print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    cat("\n")

    # Print co-skewness and co-kurtosis results

    if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for joint higher moment differences:", "\n", sep = "")
    if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for joint higher moment differences:", "\n", sep = "")

    jointtests <- rbind(object$cor12diff, object$RHS3, object$cor13diff, object$RHS4, object$RCC)
    jointtests <- round(jointtests, 4)
    rownames(jointtests) <- c("Co-Skewness", "Hyvarinen-Smith (Co-Skewness)", "Co-Kurtosis", "Hyvarinen-Smith (Co-Kurtosis)", "Chen-Chan (Co-Kurtosis)")
    colnames(jointtests) <- c("estimate", "lower", "upper")
    print.default(format( jointtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    # hoctests <- rbind(object$cor12diff, object$cor13diff)
    # hoctests <- round(hoctests, 4)
    # rownames(hoctests) <- c("Cor^2[2,1] - Cor^2[1,2]", "Cor^2[3,1] - Cor^2[1,3]" )
    # colnames(hoctests) <- c("estimate", "lower", "upper")
    # print.default(format( hoctests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    # cat("\n")

    # # Print LR approximative measures

    # if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")
    # if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs for Likelihood Ratio approximations:", "\n", sep = "")

    # LRtests <- rbind(object$RHS3, object$RHS4, object$RCC )
    # LRtests <- round(LRtests, 4)
    # rownames(LRtests) <- c("Hyvarinen-Smith (co-skewness)", "Hyvarinen-Smith (co-kurtosis)", "Chen-Chan (co-kurtosis)")
    # colnames(LRtests) <- c("estimate", "lower", "upper")
    # print.default(format( LRtests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)

    cat("\n")
    cat(paste("Number of resamples:", object$boot.args[3]))
    cat("\n")
  }

  cat("---")
  cat("\n")
  cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
  cat("\n")
  cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
  cat("\n")
  if(isFALSE(object$probtrans)){
    cat(paste("      Difference statistics > 0 suggest the model", varnames[2], "->", varnames[1], sep = " "))
  } else {
    cat(paste("      As of package version 0.2.0, difference statistics > 0 the Target model when using prob.trans = TRUE.", varnames[2], "->", varnames[1], sep = " "))
  }
  cat("\n")
  if(isTRUE(object$boot.warning)) { # 0.2.0 isTRUE for B = 0
    cat("Warning: Excess-kurtosis values of residuals have unequal signs", "\n", "        Also compute Co-Kurtosis and Hyvarinen-Smith Co-Kurtosis for", varnames[1], "->", varnames[2], "\n") }
	cat("\n")
  }


