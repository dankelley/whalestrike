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

#' Find maximum Lethality Index during a strike
#'
#' This works by finding the maximum Lethality Index encountered during
#' a simulation created by calling [strike()], and so it is important to
#' use a detailed setting for the output times. In the example, the
#' results are reported every 0.7/200 seconds (i.e. 3.5 milliseconds),
#' which is likely sufficient (see the example, where a plot is
#' used for this assessment).
#'
#' @param strike the value returned by a call to strike.
#'
#' @examples
#' library(whalestrike)
#' t <- seq(0, 0.7, length.out = 200)
#' state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0)
#' parms <- parameters()
#' s <- strike(t, state, parms)
#' # Compute the desired value and (for context) show it on a plot
#' maximumLethalityIndex(s)
#' # For context, this is how this can be done "by hand"
#' max(lethalityIndexFromStress(s[["WCF"]][["stress"]]))
#' # Show the maximum on a plot (see also the plot title)
#' plot(s, which = "lethality index")
#' abline(h=maximumLethalityIndex(s), col=2)
#'
#' @author Dan Kelley, wrapping code provided by Alexandra Mayette
#'
#' @family functions dealing with Whale Lethality index
#'
#' @return The maximum value of the Lethality Index that is
#' involved in the simulation of the ship-whale collision event.
#' This is a unitless number; see Kelley et al. (2021).
#'
#' @export
#'
#' @references
#' Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the
#' Lethality of Ship Strikes on Whales Using Simple Biophysical Models." Marine
#' Mammal Science 37, no. 1 (January 2021): 251–67.
maximumLethalityIndex <- function(strike)
{
    max(lethalityIndexFromStress(strike[["WCF"]][["stress"]]))
}
