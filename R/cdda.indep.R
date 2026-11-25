## OLD
## #' @param diff        A logical value indicating whether differences in HSIC, dCor, and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.
## #' @return A list of class \code{cdda.indep} containing the results of
## #'   CDDA independence tests for pre-specified moderator values.
## NEW 0.2.0
## Edited in the roxygen block below: @param diff, @return, new @details.

#' @title Conditional Direction Dependence Analysis: Independence Properties
#' @description \code{cdda.indep} evaluates asymmetries of predictor-error independence of competing
#'              conditional models (\code{y ~ x * m} vs. \code{x ~ y * m}
#'              with \code{m} being a continuous or categorical moderator).
#'
#' @name cdda.indep
#'
#' @param formula     Symbolic formula of the model to be tested or an \code{lm} object.
#' @param pred        A character indicating the variable name of the predictor which serves as the outcome in the alternative model.
#' @param mod         A character indicating the variable name of the moderator.
#' @param modval      Characters or a numeric sequence specifying the moderator
#'                    values used in post-hoc probing. Possible characters include
#'                    \code{c("mean", "median", "JN")}. \code{modval = "mean"}
#'                    tests the interaction effect at the moderator values
#'                    \code{M - 1SD}, \code{M}, and \code{M + 1SD};
#'                    \code{modval = "median"} uses \code{Q1}, \code{Md},
#'                    and \code{Q3}. The Johnson-Neyman approach is applied
#'                    when \code{modval = "JN"} with conditional effects being
#'                    evaluated at the boundary values of the significance
#'                    regions. When a numeric sequence is specified, the
#'                    pick-a-point approach is used for the selected numeric values.
#' @param data        A required data frame containing the variables in the model.
#' @param hetero Determines handling of separate homoscedasticity tests:
#'   \itemize{
#'     \item \strong{For} \code{cdda.indep}: A logical value indicating
#'       whether separate homoscedasticity tests (i.e., standard and robust
#'       Breusch-Pagan tests) should be computed. Default is \code{FALSE}.
#'     \item \strong{For} \code{summary}: A logical value indicating
#'       whether separate homoscedasticity tests should be returned in the
#'       output. Default is \code{FALSE}.
#'   }
#' @param diff        A logical value indicating whether differences in HSIC, dCor, and MI values should be computed. Bootstrap confidence intervals are computed using B bootstrap samples.
#' @param nlfun Determines handling of non-linear correlation tests
#'   depending on the function used:
#'   \itemize{
#'     \item \strong{For} \code{cdda.indep}: Either a numeric value or
#'       a function of .Primitive type used for non-linear correlation
#'       tests. When numeric, the value is used in a power transformation.
#'     \item \strong{For} \code{summary}: A logical value indicating
#'       whether non-linear correlation tests should be returned in the
#'       output. Default is \code{FALSE}.
#'   }
#' @param hsic.method A character indicating the inference method for the Hilbert-Schmidt Independence Criterion. Must be one of the four specifications \code{c("gamma", "eigenvalue", "bootstrap", "permutation")}. \code{hsic.method = "gamma"} is the default.
#' @param B           Number of permutations for separate dCor tests and number of resamples when \code{hsic.method = c("bootstrap", "permutation")} or \code{diff = TRUE}.
#' @param boot.type   A character indicating the type of bootstrap confidence intervals. Must be one of the two specifications \code{c("perc", "bca")}. \code{boot.type = "perc"} is the default.
#' @param conf.level  Confidence level for bootstrap confidence intervals.
#' @param parallelize A logical value indicating whether bootstrapping is performed on multiple cores. Only used if \code{diff = TRUE}.
#' @param cores       A numeric value indicating the number of cores. Only used if \code{parallelize = TRUE}.
#' @param robust      A logical value indicating whether Siegel's (1982) repeated median estimator should be used for model estimation. If \code{robust = TRUE} repeated median estimation is applied for the causally competing models, otherwise ordinary least squares (OLS) estimation is used. Moderator, interaction, and covariate terms are partialled out with OLS in both cases.
#'
#' @return A list of class \code{cdda.indep} containing the results of
#'   CDDA independence tests for pre-specified moderator values. When
#'   \code{diff = TRUE} the element \code{cdda.diff} holds one matrix of HSIC,
#'   dCor and MI differences (estimate and CI limits) per moderator value.
#'
#' @details At each moderator value the target model is estimated from
#'   variables purified against the moderator and the product of moderator and
#'   predictor. The alternative model is estimated from variables purified
#'   against the moderator and the product of moderator and outcome.
#'   Difference statistics compare these two models: independence of the
#'   alternative model minus independence of the target model, so values > 0
#'   suggest the target model.
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
#' result <- cdda.indep(m,
#'                      pred = "x",
#'                      mod = "z",
#'                      modval = c(-1, 1),
#'                      data = d,
#'                      hetero = TRUE,
#'                      diff = TRUE,
#'                      parallelize = TRUE,
#'                      cores = 2,
#'                      nlfun = 2,
#'                      B = 2)
#' # Note: Only 2 bootstrap samples are generated here to lower computation time.
#' # Permutation p-values cannot fall below 1 / (B + 1), so the dCor p-values
#' # below are bounded at 0.33; raise B for statistical inference.
#'
#' print(result)
#' summary(result, hetero = TRUE)
#' plot(result, stat = "dcor.diff")
#'
#' \dontrun{
#'
#' result <- cdda.indep(m,
#'                      pred = "x",
#'                      mod = "z",
#'                      modval = c(-1, 1),
#'                      data = d,
#'                      hetero = TRUE,
#'                      diff = TRUE,
#'                      parallelize = TRUE,
#'                      cores = 2,
#'                      nlfun = 2,
#'                      B = 500)
#'
#' summary(result, hsic.diff = TRUE, dcor.diff = TRUE, mi.diff = TRUE)
#' plot(result, stat = "hsic.diff")
#' }
#'
#' @references Wiedermann, W., & von Eye, A. (2025). \emph{Direction Dependence Analysis: Foundations and Statistical Methods}. Cambridge, UK: Cambridge University Press.
#' @seealso \code{\link{dda.indep}} for an unconditional version.
#'
#' @export
#' @rdname cdda.indep
cdda.indep <- function(
    formula = NULL,
    pred = NULL,
    mod = NULL,
    modval = NULL,
    data = list(),
    hetero = FALSE,
    diff = FALSE,
    nlfun = NULL,
    hsic.method = "gamma",
    B = 200,
    boot.type = "perc",
    conf.level = 0.95,
    parallelize = FALSE,
    cores = 1,
    robust = FALSE,
    ...
){

  if (length(data) == 0) stop( "Please specify a data frame." )

  if (!inherits(formula, "formula") ) {
    X <- if (is.matrix(formula$x) ) formula$x
    else model.matrix(terms(formula), model.frame(formula) )
    y <- if (is.vector(formula$y) ) formula$y
    else model.response(model.frame(formula))

    delete.pred <- which( colnames(X) == pred ) # get position of tentative predictor
    if ( length(delete.pred) == 0 ) stop( "Specified predictor not found in the target model." )

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
    moderator <- X[, delete.mod] # moderator

    X   <- X[, -c(delete.pred, delete.mod, delete.interaction)] # model matrix with only covariates

    if (!is.matrix(X)) X <- as.matrix(X)
  }

  #compute competing models
  target.model <- round(summary(lm(y ~ x + moderator * x + X))$coefficients, 4)
  alternate.model <- round(summary(lm(x ~ y + moderator * y + X))$coefficients, 4)

  response.name <- all.vars(formula(formula))[1]  # get name of response variable

  nam.tar <- rownames(target.model)
  nam.tar <- gsub("moderator", mod, rownames(target.model))
  nam.tar <- gsub("x", pred, nam.tar)
  rownames(target.model) <- sub("^X", "", nam.tar)

  nam.alt <- gsub("moderator", mod, rownames(alternate.model))
  nam.alt <- gsub("y", response.name, nam.alt)
  rownames(alternate.model) <- sub("^X", "", nam.alt)

  ## OLD (none)
  ## NEW 0.2.0
  # difference statistics: independence of the alternative model
  # minus independence of the target model at the same moderator value
  fit.cm <- function(f) {
    if (isTRUE(robust)) RobustLinearReg::siegel_regression(f) else lm(f)
  }

  boot.cdda <- function(dat, g){
    dat <- dat[g, ]
    ## NEW 0.2.0
    ## mi.dat <- dat[!duplicated(g), ]  # distinct rows; repeated rows sit at distance zero
    hsic <- function(a, b) dHSIC::dhsic.test(a, b, method = "gamma")$statistic
    c(hsic(dat$ry.alt, dat$e.alt) - hsic(dat$rx.tar, dat$e.tar),
      dccpp::dcor(dat$ry.alt, dat$e.alt) - dccpp::dcor(dat$rx.tar, dat$e.tar),
      # OLD: MI difference from marginal entropies (Hyvarinen & Smith, 2013). The
      # joint-entropy term of the two purified pairs is not included.
      # NEW 0.2.0: MI estimated directly for each model (Kraskov et al., 2004).
      ## OLD
      (max.entropy(dat$ry.alt) + max.entropy(dat$e.alt)) -
        (max.entropy(dat$rx.tar) + max.entropy(dat$e.tar)))
      ## NEW 0.2.0
      ## mi.ksg(mi.dat$ry.alt, mi.dat$e.alt) - mi.ksg(mi.dat$rx.tar, mi.dat$e.tar))
  }

  max.entropy <- function(x){
    sdx <- sd(x)
    x  <- as.vector(scale(x))
    k1 <- 79.047
    k2 <- 7.412889
    g  <- 0.37457
    GaussE <- log(2*pi)/2+1/2
    NegE <- k1 * (mean(log(cosh(x))) - g)^2 + k2 * mean(x * exp(-x^2/2))^2
    GaussE - NegE + log(sdx)
  }

  ## NEW 0.2.0
  ## Kraskov, Stoegbauer, & Grassberger (2004) estimator 1 (eq. 8), maximum norm
  ## mi.ksg <- function(u, e, k = 3){
  ##   u <- as.vector(scale(u))
  ##   e <- as.vector(scale(e))
  ##   n <- length(u)
  ##   nu <- ne <- numeric(n)
  ##   for (i in 1:n) {
  ##     du <- abs(u - u[i])
  ##     de <- abs(e - e[i])
  ##     eps <- sort(pmax(du, de), partial = k + 1)[k + 1]  # distance to the k-th neighbour
  ##     nu[i] <- sum(du < eps) - 1
  ##     ne[i] <- sum(de < eps) - 1
  ##   }
  ##   digamma(k) + digamma(n) - mean(digamma(nu + 1) + digamma(ne + 1))
  ## }

  cdda.diff <- function(rx.tar, ry.tar, rx.alt, ry.alt){
    dat <- data.frame(rx.tar = rx.tar,
                      e.tar  = resid(fit.cm(ry.tar ~ rx.tar)),  # target model
                      ry.alt = ry.alt,
                      e.alt  = resid(fit.cm(rx.alt ~ ry.alt)))  # alternative model


    if (parallelize) {
      boot.res <- boot::boot(dat, boot.cdda, R = B, parallel = "snow", ncpus = cores)
    } else {
      boot.res <- boot::boot(dat, boot.cdda, R = B)
    }

    type <- boot.type
    if (type == "bca" && any(is.na(boot::empinf(boot.res)))) {
      warning("Acceleration constant cannot be calculated. Falling back to percentile bootstrap method. ",
              "Consider increasing the number of resamples (B) for more stable results.", call. = FALSE)
      type <- "perc"
    }

    ci <- suppressWarnings(t(sapply(1:3, function(k)
      boot::boot.ci(boot.res, conf = conf.level, type = type, index = k)[[4]][4:5])))

    out <- cbind(boot.res$t0, ci)
    dimnames(out) <- list(c("HSIC", "dCor", "MI"), c("estimate", "lower", "upper"))
    out
  }

  diff.temp <- list()
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
      ## NEW 0.2.0: diff calculated outside of core dda function, so diff = F
      indep.temp.yx[[i]] <- unclass(dda.indep(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = FALSE, hetero = hetero,
                                              parallelize = parallelize, cores = cores,
                                              robust = robust # 0.2.0 robust estimation of the competing models
      ))
      indep.temp.xy[[i]] <- unclass(dda.indep(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = FALSE, hetero = hetero,
                                              parallelize = parallelize, cores = cores,
                                              robust = robust # 0.2.0 robust estimation of the competing models
      ))

      if (diff) diff.temp[[i]] <- cdda.diff(rx.aux.tar, ry.aux.tar, rx.aux.alt, ry.aux.alt)
    }

    ## OLD
    ## names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", unique(data[,mod]), sep = ' ')
    ## NEW 0.2.0
    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", levels(data[,mod]), sep = ' ')
    if (diff) names(diff.temp) <- names(indep.temp.yx)
  }

  else {

    ## OLD
    ## if (is.numeric(modval)){
    ##
    ##   if( is.null(modval) ) stop("No moderation values identified for pick-a-point approach.")
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

      jnoutput <- unclass(interactions::johnson_neyman( lm(y ~ x + moderator + moderator:x + X, data = mdat), pred = x,
                                          modx = moderator, control.fdr = FALSE ) )
      twobounds <- jnoutput[["bounds"]]

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
      aux.yx.tar <- lm(y ~ X + modmat[,i] + modmat[,i]:x)
      aux.xy.tar <- lm(x ~ X + modmat[,i] + modmat[,i]:x)

      #Alternative model: y -> x
      aux.yx.alt <- lm(y ~ X + modmat[,i] + modmat[,i]:y)
      aux.xy.alt <- lm(x ~ X + modmat[,i] + modmat[,i]:y)

      #Target model residuals
      ry.aux.tar <- resid(aux.yx.tar)
      rx.aux.tar <- resid(aux.xy.tar)

      #Alternative model residuals
      ry.aux.alt <- resid(aux.yx.alt)
      rx.aux.alt <- resid(aux.xy.alt)

      #DDA Implementation
      ## NEW 0.2.0: diff = FALSE, as it's calculated outside of the package
      indep.temp.yx[[i]] <- unclass(dda.indep(ry.aux.tar ~ rx.aux.tar, pred = "rx.aux.tar",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = FALSE,
                                              hetero = hetero, parallelize = parallelize, cores = cores,
                                              robust = robust # 0.2.0 robust estimation of the competing models
      ))
      indep.temp.xy[[i]] <- unclass(dda.indep(rx.aux.alt ~ ry.aux.alt, pred = "ry.aux.alt",
                                              hsic.method = hsic.method, nlfun = nlfun,
                                              B = B, boot.type = boot.type,
                                              conf.level = conf.level, diff = FALSE,
                                              hetero = hetero, parallelize = parallelize, cores = cores,
                                              robust = robust # 0.2.0 robust estimation of the competing models
      ))

      if (diff) diff.temp[[i]] <- cdda.diff(rx.aux.tar, ry.aux.tar, rx.aux.alt, ry.aux.alt)

    }


    ## OLD
    ## names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", colnames(modmat), sep = ' ')
    ## NEW 0.2.0
    names(indep.temp.xy) <- names(indep.temp.yx) <- paste("ModVal", colnames(modmat), sep = ' ')
    if (diff) names(diff.temp) <- names(indep.temp.yx)

  }

  response.name <- all.vars(formula(formula))[1]  # get name of response variable

  cdda.output <- list()
  cdda.output[[1]] <- indep.temp.yx
  cdda.output[[2]] <- indep.temp.xy
  cdda.output[[3]] <- list("target.model" = target.model, "alternate.model" = alternate.model)
  cdda.output[[4]] <- list("response.name" = response.name,
                           "pred.name"= pred,
                           "mod.name" = mod,
                           "mod.levels" = moderator.levels,
                           "mod.data" = data[,mod]
                          )
  cdda.output[[5]] <- list(
                          "function.call" = match.call(),
                          "function.name" = as.character(match.call()[[1]]),
                          "all.args" = as.list(match.call())[-1],  # Remove function name
                          "environment" = parent.frame(),
                          "formula" = if(exists("formula")) formula else NULL,
                          "data.name" = deparse(substitute(data))
                          )

  ## OLD
  ## names(cdda.output) <- c("cdda.target", "cdda.alternative", "models",
  ##                         "df.original", "call.info")
  ## NEW 0.2.0
  cdda.output[[6]] <- list(diff = diff.temp,
                           boot.args = c(boot.type, conf.level, B))

  names(cdda.output) <- c("cdda.target", "cdda.alternative", "models",
                          "df.original", "call.info", "cdda.diff")
  class(cdda.output) <- "cdda.indep"
  return(cdda.output)
}

#' @title Print Method for \code{cdda.indep} Objects.
#' @description \code{print} returns the output of standard linear model coefficients for causally competing target and alternative models.
#' @param x     An object of class \code{cdda.indep} when using \code{print} or \code{plot}.
#' @param ...   Additional arguments to be passed to the function.
#'
#' @export
#' @rdname cdda.indep
#' @method print cdda.indep
print.cdda.indep <- function(x, ...){
  cdda.output <- x
  cat("\n")
  cat(paste0("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Target Model", "\n")
  cat("\n")
  print(round(cdda.output[[3]]$target.model, 4))

  cat(paste0("-----------------------------------------------------"))
  cat("\n")

  cat("OLS Summary: Alternative Model", "\n")
  cat("\n")
  print(round(cdda.output[[3]]$alternate.model, 4))
  cat("\n")
}

