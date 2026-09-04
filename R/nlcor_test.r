#' @title Non-linear correlation (nlcor) tests
#' @description The function \code{nlcor.test} computes non-linear
#'              correlation tests between two variables. The function
#'              returns the correlation coefficient, t-value, degrees
#'              of freedom, and p-value for three different tests:
#'              (1) Correlation between \code{x^fun} and \code{y},
#'              (2) Correlation between \code{x} and \code{y^fun},
#'              and (3) Correlation between \code{x^fun} and \code{y^fun},
#'			        with \code{fun} being a pre-selected non-linear function.
#'
#' @returns     A list of class \code{dda.nlcor} containing the results
#'              of the non-linear correlation tests.
#'
#' @param x     numeric vectors of data values. Must have the same length as y.
#' @param y     numeric vectors of data values. Must have the same length as x.
#' @param fun   a numeric value or a function of .Primitive type used to compute
#'              non-linear correlation tests. When \code{fun} is numeric,
#'              the value is used in a power transformation.
#'
#' @keywords internal
#' @noRd
#'

nlcor.test <- function(x, y, fun, fname = NULL)
  {
  varnames <- c(deparse(substitute(x)), deparse(substitute(y)))

  if (length(x) != length(y)) stop("Variables must have the same length")
  n <- length(x)

  x <- as.vector(scale(x))
  y <- as.vector(scale(y))

  if (is.numeric(fun)) {
    func <- if (!is.null(fname)) fname else as.character(fun)
    r1 <- cor(x^fun, y)
    r2 <- cor(x, y^fun)
    r3 <- cor(x^fun, y^fun)

    if (any(is.na(c(r1, r2, r3))) || any(is.nan(c(r1, r2, r3)))) {
      x <- x + abs(min(x)) + 0.1
      y <- y + abs(min(y)) + 0.1

      r1 <- cor(x^fun, y)
      r2 <- cor(x, y^fun)
      r3 <- cor(x^fun, y^fun)
    }
  } else {
    func <- if (!is.null(fname)) fname else paste(substitute(fun))
    test.run <- suppressWarnings(c(fun(x), fun(y)))

    if (any(is.na(test.run)) || any(is.nan(test.run))) {
      x <- x + abs(min(x)) + 0.1
      y <- y + abs(min(y)) + 0.1
    }

    r1 <- cor(fun(x), y)
    r2 <- cor(x, fun(y))
    r3 <- cor(fun(x), fun(y))
  }

  tval1 <- r1 * sqrt((n - 2) / (1 - r1^2))
  tval2 <- r2 * sqrt((n - 2) / (1 - r2^2))
  tval3 <- r3 * sqrt((n - 2) / (1 - r3^2))

  pval1 <- pt(abs(tval1), df = n - 2, lower.tail = FALSE) * 2
  pval2 <- pt(abs(tval2), df = n - 2, lower.tail = FALSE) * 2
  pval3 <- pt(abs(tval3), df = n - 2, lower.tail = FALSE) * 2

  output <- list(
    t1 = c(r1, tval1, n - 2, pval1),
    t2 = c(r2, tval2, n - 2, pval2),
    t3 = c(r3, tval3, n - 2, pval3),
    func = func,
    varnames = varnames
  )

  class(output) <- "dda.nlcor"
  output
}

#' @title Print Method for \code{nlcor.test} Objects
#' @description Displays non-linear correlation results between two variables.
#' @param x     An object of class \code{nlcor.test}.
#' @param ...   Additional arguments to be passed to the function.
#'
#' @returns \code{nlcor.test} test statistics and p-values.
#'
#' @rdname nlcor
#' @method print nlcor.test
#'
#' @keywords internal
#' @noRd
#'
print.dda.nlcor <- function(x, ...) {
  varnames <- x$varnames
  sigtests <- rbind(x$t1, x$t2, x$t3)
  sigtests <- round(sigtests, 4)

  cat("\n")
  if (is.na(suppressWarnings(as.numeric(x$func)))) {
    cat(paste("Non-linear correlation tests:", x$func, "transformation"))
    rownames(sigtests) <- c(
      paste("Cor[", x$func, "(", varnames[1], "), ", varnames[2], "]", sep = ""),
      paste("Cor[", varnames[1], ", ", x$func, "(", varnames[2], ")]", sep = ""),
      paste("Cor[", x$func, "(", varnames[1], "), ", x$func, "(", varnames[2], ")]", sep = "")
    )
  } else {
    cat(paste("Non-linear correlation tests: Power transformation using", x$func))
    rownames(sigtests) <- c(
      paste("Cor[", varnames[1], "^", x$func, ", ", varnames[2], "]", sep = ""),
      paste("Cor[", varnames[1], ", ", varnames[2], "^", x$func, "]", sep = ""),
      paste("Cor[", varnames[1], "^", x$func, ", ", varnames[2], "^", x$func, "]", sep = "")
    )
  }
  cat("\n")

  colnames(sigtests) <- c("estimate", "t-value", "df", "Pr(>|t|)")
  print.default(
    format(sigtests, digits = max(3L, getOption("digits") - 3L)),
    print.gap = 2L,
    quote = FALSE
  )
}


