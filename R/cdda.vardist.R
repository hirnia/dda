## OLD (none)
## NEW 0.2.0
## Edited in the roxygen block below: new @details describing the coskew,
## cokurt, rhs, rcc, and rtanh statistics, and Chen & Chan (2013) and
## Hyvarinen & Smith (2013) in @references.

#' @title Conditional Direction Dependence Analysis: Variable Distributions
#' @description \code{cdda.vardist} evaluates variable distributions of competing
#'              conditional models (\code{y ~ x * m} vs. \code{x ~ y * m} with \code{m}
#'              being a continuous or categorical moderator).
#' @name cdda.vardist
#'
#' @param formula     Symbolic formula of the model to be tested or an \code{lm} object.
#' @param pred        A character indicating the variable name of the predictor
#'                    which serves as the outcome in the alternative model.
#' @param mod         A character indicating the variable name of the moderator.
#' @param modval      Characters or a numeric sequence specifying the moderator
#'                    values used in post-hoc probing. Possible characters
#'                    include \code{c("mean", "median", "JN")}. \code{modval = "mean"}
#'                    tests the interaction effect at the moderator values
#'                    \code{M - 1SD}, \code{M}, and \code{M + 1SD};
#'                    \code{modval = "median"} uses \code{Q1}, \code{Md},
#'                    and \code{Q3}. The Johnson-Neyman approach is applied
#'                    when \code{modval = "JN"} with conditional effects being
#'                    evaluated at the boundary values of the significance
#'                    regions. When a numeric sequence is specified, the
#'                    pick-a-point approach is used for the selected numeric values.
#' @param data        A required data frame containing the variables in the model.
#' @param conf.level  Confidence level for bootstrap confidence intervals.
#' @param B           Number of bootstrap samples.
#' @param boot.type   A character indicating the type of bootstrap confidence intervals. Must be one of the two values \code{c("perc", "bca")}. \code{boot.type = "perc"} is the default.
#'
#' @returns An object of class \code{cdda.vardist} containing the results of
#'          conditional direction dependence tests of variable distributions.
#'
#' @details At each moderator value \code{dda.vardist} is applied twice. In the
#'   target panel \code{x} and \code{y} are the predictor and outcome purified
#'   against the moderator and the product of moderator and predictor; values
#'   > 0 suggest \code{x -> y}. In the alternative panel the roles of the two
#'   variables switch (purified against the moderator and the product of
#'   moderator and outcome); values > 0 suggest \code{y -> x}.
#'
#'   The statistics returned by \code{summary} and drawn by \code{plot} are
#'   computed from the standardized purified variables, with \eqn{\rho}{rho}
#'   denoting their Pearson correlation:
#'   \itemize{
#'     \item \code{coskew}: co-skewness difference
#'       \eqn{cor_{21}(x, y)^2 - cor_{12}(x, y)^2}{cor21(x, y)^2 - cor12(x, y)^2},
#'       where \eqn{cor_{ij}(x, y) = E[(x - \bar{x})^i (y - \bar{y})^j] / (s_x^i s_y^j)}{corij(x, y) = E[(x - mean(x))^i (y - mean(y))^j] / (sd(x)^i sd(y)^j)}
#'       is the standardized co-moment of order \eqn{(i, j)}{(i, j)}.
#'     \item \code{cokurt}: co-kurtosis difference
#'       \eqn{cor_{31}(x, y)^2 - cor_{13}(x, y)^2}{cor31(x, y)^2 - cor13(x, y)^2},
#'       multiplied by the sign of the excess kurtosis of \code{x}.
#'     \item \code{rhs}: Hyvarinen-Smith co-skewness, a skewness-based
#'       approximation of the likelihood ratio of the two models (Hyvarinen &
#'       Smith, 2013), \eqn{\rho E(x^2 y - x y^2)}{rho E(x^2 y - x y^2)}. Before
#'       it is computed \code{x} and \code{y} are multiplied by the sign of
#'       their skewness so that both are right-skewed. It is informative when
#'       the variables are skewed.
#'     \item \code{rcc}: Chen-Chan co-kurtosis, a kurtosis-based measure built
#'       from the fourth-order cross-cumulants
#'       \eqn{C_{xy} = E(x^3 y) - 3 \rho var(x)}{Cxy = E(x^3 y) - 3 rho var(x)} and
#'       \eqn{C_{yx} = E(x y^3) - 3 \rho var(y)}{Cyx = E(x y^3) - 3 rho var(y)}
#'       (Chen & Chan, 2013), \eqn{(C_{xy} + C_{yx})(C_{xy} - C_{yx}) = C_{xy}^2 - C_{yx}^2}{(Cxy + Cyx)(Cxy - Cyx) = Cxy^2 - Cyx^2}.
#'       It is informative when the variables are symmetric but have nonzero
#'       excess kurtosis.
#'     \item \code{rtanh}: Hyvarinen-Smith hyperbolic tangent, an
#'       approximation of the likelihood ratio that uses \eqn{\tanh}{tanh} as
#'       nonlinearity (Hyvarinen & Smith, 2013),
#'       \eqn{\rho E(x \tanh(y) - \tanh(x) y)}{rho E(x tanh(y) - tanh(x) y)}.
#'       It is designed for symmetric variables with positive excess kurtosis.
#'   }
#'
#' @examples
#' set.seed(321)
#' n <- 700
#'
#' ## --- generate moderator
#'
#' z <- sort(rnorm(n))
#' z1 <- z[z <= 0]
#' z2 <- z[z > 0]
#'
#' ## --- x -> y when z <= 0
#'
#' x1 <- rchisq(length(z1), df = 4) - 4
#' e1 <- rnorm(length(z1), sd = sqrt(6))
#' y1 <- 0.5 * x1 + e1
#'
#' ## --- y -> x when z > 0
#'
#' y2 <- rchisq(length(z2), df = 4) - 4
#' e2 <- rnorm(length(z2), sd = sqrt(6))
#' x2 <- 0.25 * y2 + e2
#'
#' y <- c(y1, y2); x <- c(x1, x2)
#'
#' d <- data.frame(x, y, z)
#'
#' m <- lm(y ~ x * z, data = d)
#'
#' result <- cdda.vardist(m, pred = "x", mod = "z", B = 50,
#'                       modval = c(-1, 1), data = d)
#' # Note: Only 50 bootstrap samples are generated here to lower computation time
#'
#' print(result)
#' summary(result, skew = FALSE, kurt = FALSE, coskew = TRUE)
#' plot(result, stat = "rhs")
#'
#' \dontrun{
#'
#' result <- cdda.vardist(m, pred = "x", mod = "z", B = 500,
#'                        modval = c(-1, 1), data = d)
#'
#' summary(result, coskew = TRUE, cokurt = TRUE)
#' plot(result, stat = "rtanh")
#' }
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK: Cambridge University Press.
#'
#' Chen, Z., & Chan, L. (2013). Causality in linear nongaussian acyclic models
#' in the presence of latent Gaussian confounders. \emph{Neural Computation,
#' 25}(6), 1605-1641. \doi{10.1162/NECO_a_00444}
#'
#' Hyvarinen, A., & Smith, S. M. (2013). Pairwise likelihood ratios for
#' estimation of non-Gaussian structural equation models. \emph{Journal of
#' Machine Learning Research, 14}, 111-152.
#'
#' @seealso \code{\link{dda.vardist}} for an unconditional version.
#'
#' @export
#' @rdname cdda.vardist
cdda.vardist <- function(formula,
                         pred = NULL,
                         mod = NULL,
                         data = list(),
                         modval = NULL,
                         B = 200,
                         boot.type = "perc",
                         conf.level = 0.95
                         ## Optional robust estimation (not active): add the argument
                         ## , robust = FALSE
                        ){

  if (length(data) == 0) stop("Please specify a data frame.")

  if (!inherits(formula, "formula") ) {
    X <- if (is.matrix(formula$x) ) formula$x
    else model.matrix(terms(formula), model.frame(formula) )
    y <- if (is.vector(formula$y) ) formula$y
    else model.response(model.frame(formula))

    delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
    if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

    if(is.factor(data[,mod])){
      delete.mod <- which(colnames(X) %in% paste(mod, levels(data[,mod]), sep ="")) # get position of moderator
    }
    else if(is.numeric(data[,mod])){
      delete.mod <- which(colnames(X) == mod)
    }

    else { stop( "Specified moderator not found in the target model." ) }
    ## OLD (none)
    ## NEW 0.2.0
    if ( length(delete.mod) == 0 ) stop( "Specified moderator not found in the target model." )

    delete.interaction <- grep(":", colnames(X))
    if ( length(delete.interaction) == 0 ) stop( "No interaction term specified in the target model." )

    x   <- X[, delete.pred] # tentative predictor
    moderator <- X[, delete.mod]  # moderator

    X <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates
    if (!is.matrix(X)) X <- as.matrix(X)

  }

  else {
    mf <- model.frame(formula, data = data)
    y  <- model.response(mf)   # tentative outcome
    X  <- model.matrix(formula, data = data)

    delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
    if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

    ## OLD
    ## if(is.factor(data[,mod])){
    ##   delete.mod <- which(colnames(X) %in% paste(mod, levels(data[,mod]), sep ="")) # get position of moderator
    ## }
    ## if(is.numeric(data[,mod])){
    ##   delete.mod <- which(colnames(X) == mod)
    ## }
    ## if (  is.null(delete.mod) == TRUE ) stop( "Specified moderator not found in the target model." )
    ## NEW 0.2.0
    delete.mod <- integer(0)
    if(is.factor(data[,mod])){
      delete.mod <- which(colnames(X) %in% paste(mod, levels(data[,mod]), sep ="")) # get position of moderator
    }
    if(is.numeric(data[,mod])){
      delete.mod <- which(colnames(X) == mod)
    }
    if ( length(delete.mod) == 0 ) stop( "Specified moderator not found in the target model." )

    delete.interaction <- grep(":", colnames(X))
    if ( length(delete.interaction) == 0 ) stop( "No interaction term specified in the target model." )

    x   <- X[,  delete.pred] # tentative predictor
    moderator <- X[, delete.mod]  # moderator

    X   <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates

    if (!is.matrix(X)) X <- as.matrix(X)
  }

  # compute competing models

  target.model <- round(summary(lm(y ~ x + moderator*x + X))$coefficients, 4)

  alternate.model <- round(summary(lm(x ~ y + moderator*y + X))$coefficients, 4)

  response.name <- all.vars(formula(formula))[1]  # get name of response variable

  nam.tar <- rownames(target.model)
  nam.tar <- gsub("moderator", mod, rownames(target.model))
  nam.tar <- gsub("x", pred, nam.tar)
  rownames(target.model) <- sub("^X", "", nam.tar)

  nam.alt <- gsub("moderator", mod, rownames(alternate.model))
  nam.alt <- gsub("y", response.name, nam.alt)
  rownames(alternate.model) <- sub("^X", "", nam.alt)

  if (is.factor(data[,mod])) {

    ## OLD
    ## moderator.levels <- unique(data[,mod])
    ## NEW 0.2.0
    moderator.levels <- levels(data[,mod])

    indep.temp.yx <- list()
    indep.temp.xy <- list()

    vardist.temp.yx <- list()
    vardist.temp.xy <- list()

    for (i in (levels(data[,mod]))) {

      moderator.rlvl <- relevel(data[,mod], ref = i)

      #Target model: x -> y
      design.tar <- model.matrix(~ X + moderator.rlvl + moderator.rlvl:x + x)

      design.tar <- as.data.frame(design.tar[,-1])
      design.tar$x <- NULL

      ## Optional robust estimation, option 2 (not active): robust partialling of
      ## the moderator, interaction and covariate terms; the same swap applies to
      ## aux.yx.alt and aux.xy.alt below and needs MASS in Imports
      ## fit.aux <- function(f) if (isTRUE(robust)) MASS::rlm(f, maxit = 50) else lm(f)
      ## aux.yx.tar <- fit.aux(y ~ -1 + as.matrix(design.tar))
      ## aux.xy.tar <- fit.aux(x ~ -1 + as.matrix(design.tar))
      aux.yx.tar <- lm(y ~ -1 + as.matrix(design.tar))
      aux.xy.tar <- lm(x ~ -1 + as.matrix(design.tar))

      ry.aux.tar <- resid(aux.yx.tar)
      rx.aux.tar <- resid(aux.xy.tar)

      #Alternative model: y -> x
      design.alt<- model.matrix(~ X + moderator.rlvl + moderator.rlvl:y + y)
      design.alt <- as.data.frame(design.alt[,-1])
      design.alt$y <- NULL

      aux.yx.alt <- lm(y ~ -1 + as.matrix(design.alt))
      aux.xy.alt <- lm(x ~ -1 + as.matrix(design.alt))

      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)

      #DDA Implementation
      ## Optional robust estimation, option 1 (not active): add robust = robust to both dda.vardist calls below
      indep.temp.yx[[i]] <- unclass(dda.vardist(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                                B = B, boot.type = boot.type, conf.level = conf.level
      ))
      indep.temp.xy[[i]] <- unclass(dda.vardist(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                                 B = B, boot.type = boot.type, conf.level = conf.level
      ))
    }

    ## OLD
    ## names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", unique(data[,mod]), sep = ' ')
    ## NEW 0.2.0
    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", levels(data[,mod]), sep = ' ')
  }

  else {

    ## OLD
    ## if (is.numeric(modval)){
    ##
    ##   #if( is.null(modval) ) stop("No moderation values identified for pick-a-point approach.")
    ## NEW 0.2.0
    if (is.null(modval)) stop("Please specify modval: 'mean', 'median', 'JN', or numeric moderator values.")

    if (is.numeric(modval)){

      modmat <- data.frame( matrix(NA, length(moderator), length(modval)) )

      for(i in seq_along(modval)){ modmat[,i] <- moderator - modval[i] }  # compute transformed moderator for each value
      colnames(modmat) <- modval
    }

    else if (modval == "mean"){

      modmat <- data.frame( (moderator - (mean(moderator) - sd(moderator)) ),
                            (moderator - mean(moderator) ),
                            (moderator - (mean(moderator) + sd(moderator)) )
      )
      colnames(modmat) <- c("M-1SD", "M", "M+1SD")
    }

    else if (modval == "median"){
      modmat <- data.frame( "Q1" = (moderator - quantile(moderator, probs = c(0.25)) ),
                            "Md" = (moderator - quantile(moderator, probs = c(0.50)) ),
                            "Q3" = (moderator - quantile(moderator, probs = c(0.75)) )
      )
    }


    else if (modval == "JN"){

      mdat <- data.frame(x, y, X, moderator)

      unwrpjn <- unclass(interactions::johnson_neyman( lm(y ~ x + moderator + moderator:x + X, data = mdat), pred = x,
                                         modx = moderator, control.fdr = FALSE ) )
      twobounds <- unwrpjn[["bounds"]]

      ## OLD (none)
      ## NEW 0.2.0
      values <- NA
      if(twobounds[1] < min(moderator, na.rm = TRUE) &
         twobounds[2] > max(moderator, na.rm = TRUE)) values <- NA

      if(twobounds[1] > min(moderator, na.rm = TRUE) &
         twobounds[2] > max(moderator, na.rm = TRUE)) values <- c(min(moderator), twobounds[1])

      if(twobounds[1] < min(moderator, na.rm = TRUE) &
         twobounds[2] < max(moderator, na.rm = TRUE)) values <- c(twobounds[2], max(moderator))

      if(twobounds[1] > min(moderator, na.rm = TRUE) &
         twobounds[2] < max(moderator, na.rm = TRUE)) values <- c(min(moderator), twobounds[1], twobounds[2], max(moderator))

      ## OLD (none)
      ## NEW 0.2.0
      if (all(is.na(values))) stop("Johnson-Neyman bounds lie outside the observed moderator range. ",
                                   "Use modval = 'mean', 'median', or numeric values.")
      modmat <- data.frame( matrix(NA, length(moderator), length(values)) )
      for( i in seq_along(values) ){ modmat[,i] <- moderator - values[i] }  # compute transformed moderator for each value
      colnames(modmat) <- paste0(round(values, 2))
    }

    else stop("Invalid or no moderator value specified.")

    moderator.levels <- colnames(modmat)

    indep.temp.yx <- list()
    indep.temp.xy <- list()

    vardist.temp.yx <- list()
    vardist.temp.xy <- list()

    for ( i in 1:(ncol(modmat)) ) {

      #Target model: x -> y
      ## Optional robust estimation, option 2 (not active): robust partialling of
      ## the moderator, interaction and covariate terms with -1 for the intercept
      ## already in X; the same swap applies to aux.yx.alt and aux.xy.alt below
      ## and needs MASS in Imports
      ## fit.aux <- function(f) if (isTRUE(robust)) MASS::rlm(f, maxit = 50) else lm(f)
      ## aux.yx.tar <- fit.aux(y ~ -1 + X + modmat[,i] + modmat[,i]:x)
      ## aux.xy.tar <- fit.aux(x ~ -1 + X + modmat[,i] + modmat[,i]:x)
      aux.yx.tar <- lm(y ~ X + modmat[,i] + modmat[,i]:x)
      aux.xy.tar <- lm(x ~ X + modmat[,i] + modmat[,i]:x)

      ry.aux.tar <- resid(aux.yx.tar)
      rx.aux.tar <- resid(aux.xy.tar)

      #Alternative model: y -> x
      aux.yx.alt <- lm(y ~ X +  modmat[,i] +  modmat[,i]:y)
      aux.xy.alt <- lm(x ~ X +  modmat[,i] +  modmat[,i]:y)

      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)

      #DDA Implementation
      ## Optional robust estimation, option 1 (not active): add robust = robust to both dda.vardist calls below
      indep.temp.yx[[i]] <- unclass(dda.vardist(ry.aux.tar ~ rx.aux.tar, B = B,
                                                pred = "rx.aux.tar", boot.type = boot.type, conf.level = conf.level))
      indep.temp.xy[[i]] <- unclass(dda.vardist(rx.aux.alt ~ ry.aux.alt, B = B,
                                                pred = "ry.aux.alt", boot.type = boot.type, conf.level = conf.level))

    }

    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", colnames(modmat), sep = ' ')
  }

  cdda.var.output <- list()
  cdda.var.output[[1]] <- indep.temp.yx
  cdda.var.output[[2]] <- indep.temp.xy
  cdda.var.output[[3]] <- list("target.model" = target.model, "alternate.model" = alternate.model)
  cdda.var.output[[4]] <- list("response.name" = response.name,
                           "pred.name"= pred,
                           "mod.name" = mod,
                           "mod.levels" = moderator.levels,
                           "mod.data" = data[,mod]
                           )
  cdda.var.output[[5]] <- list(
    "function.call" = match.call(),
    "function.name" = as.character(match.call()[[1]]),
    "all.args" = as.list(match.call())[-1],  # Remove function name
    "environment" = parent.frame(),
    "formula" = if(exists("formula")) formula else NULL,
    "data.name" = deparse(substitute(data))
  )
  names(cdda.var.output) <- c("cdda.target", "cdda.alternative", "models",
                              "df.original", "call.info")
  class(cdda.var.output) <- "cdda.vardist"
  return(cdda.var.output)

}

#' @name print.cdda.vardist
#' @title       Print Method for \code{cdda.vardist} Objects
#' @description \code{print} returns the output of standard linear model coefficients
#'              for causally competing target and alternative models.
#' @param x     An object of class \code{cdda.vardist} when using \code{print} or \code{plot}.
#' @param ...   Additional arguments to be passed to the function.
#'
#' @export
#' @rdname cdda.vardist
#' @method print cdda.vardist
print.cdda.vardist <- function(x, ...){

  cdda.var.output <- x

  cat(paste("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Target Model", "\n")
  cat("\n")
  print(cdda.var.output[[3]]$target.model)

  cat(paste("-----------------------------------------------------"))
  cat("\n", "\n")

  cat("OLS Summary: Alternative Model", "\n")
  cat("\n")
  print(cdda.var.output[[3]]$alternate.model)
}
