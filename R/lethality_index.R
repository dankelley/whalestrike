# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

#' Compute lethality index, based on compression stress
#'
#' The model used for this is the logistic model, fitting observed injury/lethality
#' statistics to the base-10 logarithm of the maximum compression stress during
#' a simulated impact event.
#'
#' @param stress numerical value or vector, giving whale compression stress in Pascals.
#'
#' @return threat of injury (in range 0 to 1)
#'
#' @examples
#' lethalityIndexFromStress(parameters()$logistic$tau50) # approx. 0.5
#'
#' @author Dan Kelley
#'
#' @family functions dealing with Whale Lethality index
#'
#' @export
lethalityIndexFromStress <- function(stress)
{
    logistic <- parameters()$logistic
    1 / (1 + exp(-(log10(stress) - logistic$logStressCenter) / logistic$logStressWidth))
}

#' Compute stress, based on lethality index
#'
#' The model used for this is the logistic model, fitting observed injury/lethality
#' statistics to the base-10 logarithm of the maximum compression stress during
#' a simulated impact event.
#'
#' @param injury numerical value or vector, giving threat of injury (in range 0 to 1).
#'
#' @return whale compression stress, in Pascals.
#'
#' @examples
#' stressFromLethalityIndex(0.5) # approx. 254000 Pa, i.e. parameters()$logistic$tau50
#'
#' @author Dan Kelley
#'
#' @family functions dealing with Whale Lethality index
#'
#' @export
stressFromLethalityIndex <- function(injury)
{
    logistic <- parameters()$logistic
    10^(logistic$logStressCenter - logistic$logStressWidth * log(1 / injury - 1)) # note natural log
}

#' Compute maximum Lethality Index for a strike simulation
#'
#' @param strike the value returned by a call to strike.
#'
#' @examples
#' library(whalestrike)
#' t <- seq(0, 0.7, length.out = 200)
#' state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0)
#' parms <- parameters()
#' s <- strike(t, state, parms)
#' max(lethalityIndexFromStress(s[["WCF"]][["stress"]]))
#' maximumLethalityIndex(s)
#'
#' @author Dan Kelley, wrapping code provided by Alexandra Mayette
#'
#' @family functions dealing with Whale Lethality index
#'
#' @export
maximumLethalityIndex <- function(strike)
{
    max(lethalityIndexFromStress(s[["WCF"]][["stress"]]))
}
