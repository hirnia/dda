% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dda_indep.r
\name{dda.indep}
\alias{dda.indep}
\title{Direction Dependence Analysis: Independence Distribution}
\usage{
dda.indep(
  formula,
  pred = NULL,
  data = list(),
  nlfun = NULL,
  hetero = FALSE,
  hsic.method = "gamma",
  diff = FALSE,
  B = 200,
  boot.type = "perc",
  conf.level = 0.95,
  parallelize = FALSE,
  cores = 1,
  ...
)
}
\arguments{
\item{formula:}{symbolic formula of the model to be tested or a \code{lm} object.}

\item{pred:}{a character indicating the variable name of the predictor which serves as the outcome in the alternative model.}

\item{data:}{an optional data frame containing the variables in the model (by default variables are taken from the environment which \code{dda.indep} is called from)}

\item{nlfun:}{Either a numeric value or a function of .Primitive type used for non-linear correlation tests. When \code{nlfun} is numeric the value is used in a power tranformation.}

\item{hetero:}{A logical value indicating whether separate homoscedasticity tests (i.e., standard and robust Breusch-Pagan tests) should be computed.}

\item{hsic.method:}{a character indicating the inference method for Hilbert-Schmidt Independence Criterion (HSIC). Must be one of the four values \code{c("gamma", "eigenvalue", "boot", "permutation")}. \code{hsic.method = "gamma"} is the default.}

\item{diff:}{A logical value indicating whether differences in HSIC, Distance Correlation (dCor), and MI values should be computed. Bootstrap confidence intervals are computed using \code{B} bootstrap samples.}

\item{B:}{Number of permutations for separate dCor tests and number of resamples if \code{hsic.method = c("boot", "permutation")} or \code{diff = TRUE}}

\item{boot.type:}{A vector of character strings representing the type of bootstrap confidence intervals required. Must be one of the two values \code{c("perc", "bca")}. \code{boot.type = "perc"} is the default.}

\item{conf.level:}{confidence level for bootstrap confidence intervals}

\item{parallelize:}{A logical value indicating whether boostrapping is performed on multiple cores. Only used if \code{diff = TRUE.}}

\item{cores:}{a numeric value indicating the number of cores. Only used if parallelize = TRUE}
}
\value{
An object of class \code{ddaindep} containing the results of the independence tests.
}
\description{
\code{dda.indep} tests the skewness and kurtosis of the variables of two competing models. It also tests the difference in skewness and kurtosis to assess independence properties of the two models. The function also provides bootstrap confidence intervals for the difference in skewness and kurtosis.
}
\examples{
dda.car.indep <- dda.indep(mpg ~ wt + qsec, pred = "wt",
                                     diff = TRUE, data = mtcars)
          dda.car.indep
          #OR
          car.test <- lm(mpg ~ wt + qsec, data = mtcars)
          dda.indep(car.test, pred = "wt", diff = TRUE, data = mtcars)

}
