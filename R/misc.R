# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

#' Convert a speed in knots to a speed in m/s
#'
#' See also [mps2knot()], which is the inverse of this function.
#'
#' @param knot Speed in knots.
#'
#' @return Speed in m/s.
#'
#' @author Dan Kelley
#'
#' @examples
#' library(whalestrike)
#' knots <- seq(0, 20)
#' plot(knots, knot2mps(knots), xlab = "Speed [knots]", ylab = "Speed [m/s]", type = "l")
#'
#' @family functions dealing with units
#'
#' @export
knot2mps <- function(knot) {
    knot * 1.852e3 / 3600
}

#' Convert a speed in m/s to a speed in knots
#'
#' This is done by dividing by the factor 1.852e3/3600,
#' See also [knot2mps()], which is the inverse of this function.
#'
#' @param mps Speed in metres per second.
#'
#' @return Speed in knots.
#'
#' @author Dan Kelley
#'
#' @examples
#' library(whalestrike)
#' mps <- seq(0, 10)
#' plot(mps, mps2knot(mps), xlab = "Speed [m/s]", ylab = "Speed [knots]", type = "l")
#'
#' @family functions dealing with units
#'
#' @export
mps2knot <- function(mps) {
    mps / (1.852e3 / 3600)
}

#' Pin numerical values between stated limits
#'
#' @param x Vector or matrix of numerical values
#'
#' @param lower Numerical values of minimum value allowed; set to `NULL`
#' to avoid trimming the lower limit.
#'
#' @param upper As for `lower`, but for the upper limit.
#'
#' @return Copy of `x`, with any value that exceeds `lim` having
#' been replaced by `lim`.
#'
#' @author Dan Kelley
#'
#' @export
pin <- function(x, lower = NULL, upper = NULL) {
    # Protect the ifelse() operation from getting riled by NAs
    na <- is.na(x)
    x[na] <- 0 # value is arbitrary because changed back to NA later
    if (!is.null(lower)) {
        x <- ifelse(x > lower, x, lower)
    }
    if (!is.null(upper)) {
        x <- ifelse(x < upper, x, upper)
    }
    x[na] <- NA
    x
}

#' Calculate derivative using first difference
#'
#' The derivative is estimated as the ratio of the first-difference of `var`
#' divided by the first-difference of `time`.  To make the results
#' have the same length as `time`, the final result is appended at
#' the end.
#'
#' @param var variable.
#'
#' @param t time in seconds.
#'
#' @return Derivative estimated by using [diff()] on both `var`
#' and `time`.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom utils tail
derivative <- function(var, t) {
    res <- diff(var) / diff(t)
    c(res, tail(res, 1))
}
