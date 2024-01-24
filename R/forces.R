# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

library(deSolve)

#' Whale blubber stress-strain relationship
#'
#' This is a data frame with elements `strain` and `stress`,
#' found by digitizing (accurate to perhaps 1 percent) the curve shown in Figure 2.13
#' of Raymond (2007). It is used to develop a stress-strain relationship used
#' by [parameters()], as shown in \dQuote{Examples}.
#'
#' @examples
#' data(raymond2007)
#' attach(raymond2007)
#' # Next yields \code{a=1.64e5} Pa and \code{b=2.47}.
#' m <- nls(stress ~ a * (exp(b * strain) - 1), start = list(a = 1e5, b = 1))
#' plot(strain, stress, xaxs = "i", yaxs = "i")
#' x <- seq(0, max(strain), length.out = 100)
#' lines(x, predict(m, list(strain = x)))
#'
#' @name raymond2007
#'
#' @references
#'
#' Raymond, J. J. "Development of a Numerical Model to Predict Impact Forces on a
#' North Atlantic Right Whale during Collision with a Vessel." University of New
#' Hampshire, 2007. \url{https://scholars.unh.edu/thesis/309/}.
#'
#' @docType data
NULL

#' Create a function for stress in laminated layers
#'
#' Denoting unforced layer thickness in the \eqn{i} layer as
#' \eqn{l_i} and strain there as \eqn{\epsilon_i=\Delta l_i/l_i},
#' we may write the stress-strain relationship as
#' \deqn{\sigma = a_i*(exp(b_i*\epsilon_i)-1)}
#' for each layer, where it is assumed that stress
#' \eqn{\sigma} is equal across layers.
#' Inverting this yields
#' \deqn{\epsilon_i= ln(1 + \sigma/a_i)/b_i}
#' where \eqn{ln} is the natural logarithm.  Therefore,
#' the change \eqn{\Delta L} in the total thickness \eqn{L=\sum l_i}
#' may be written
#' \deqn{0 = \Delta L - \sum((l_i/b_i) ln(1+\sigma/a_i))}.
#' Note that zero-thickness layers are removed from the calculation,
#' to avoid spurious forces.
#'
#' This expression is not easily inverted to get
#' \eqn{\sigma} in terms of \eqn{\Delta L}
#' but it may be solved
#' easily for particular numerical values, using [uniroot()].
#'
#' This is done for a sequence of `N` values of strain \eqn{\epsilon}
#' that range from 0 to 1. Then [approxfun()] is used to create
#' a piecewise-linear representation of the relationship between \eqn{\sigma} and \eqn{\Delta L},
#' which becomes the return value of the present function.
#' (The purpose of using a piecewise-linear representation to reduce
#' computation time.)
#'
#' @param l vector of layer thicknesses
#'
#' @param a vector of multipliers
#'
#' @param b vector of e-fold parameters
#'
#' @param N integer specifying how many segments to use in the spline
#'
#' @return A piecewise-linear function, created with [approxfun()],
#' that returns stress as a function of total strain of the
#' system of compressing layers. For the purposes of the whale-strike
#' analysis, the strain should be between 0 and 1, i.e. there is
#' no notion of compressing blubber, etc. to negative thickness.
#'
#' @examples
#' library(whalestrike)
#' # Set blubber parameters for each layer, to see if
#' # we recover the raymond2007 data.
#' param <- parameters(a = rep(1.64e5, 4), b = rep(2.47, 4))
#' x <- seq(0, 0.5, length.out = 100)
#' y <- param$stressFromStrain(x)
#' plot(x, y, type = "l", lwd = 4, col = "gray")
#' data("raymond2007")
#' points(raymond2007$strain, raymond2007$stress, col = 2)
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom stats approxfun uniroot
stressFromStrainFunction <- function(l, a, b, N = 1000) {
    use <- rep(TRUE, length(l))
    fcn <- function(sigma) {
        #> cat("fcn(): l=", paste(l,collapse=" "), ", a=", paste(a,collapse=" "), ", b=", paste(b,collapse=" "), "\n")
        DL - sum((l[use] / b[use]) * log(1 + sigma / a[use]))
    }
    L <- sum(l)
    sigma <- rep(NA, N)
    epsilon <- seq(0, 1, length.out = N)
    # We will limit the epsilon in any given layer to 1, i.e. we don't
    # permit layers to be compressed to negative thickness. This is expressed
    # by setting layer-by-layer conditions on maximum stress.
    sigmaMax <- max(a * (exp(b) - 1))
    #> cat("sigmaMax=", paste(sprintf("%.3g", sigmaMax), collapse=" "), "\n")
    sigmaLowerLimit <- 0
    sigmaUpperLimit <- 2 * sigmaMax
    for (i in seq_along(epsilon)) {
        # debug cat(sprintf("i=%3d, epsilon=%10.5f, ", i, epsilon[i]), ", use=", paste(use, collapse=" "), "\n")
        DL <- epsilon[i] * L
        #> cat("LINE 414. i=", i, ", about to call uniroot(fcn,...); use=",
        # paste(use, collapse=" "), "fcn(sigmaLowerLimit)=", fcn(sigmaLowerLimit),
        # ", fcn(big)=", fcn(sigmaUpperLimit), ", fcn(10*big)=", fcn(10*sigmaUpperLimit),
        # ", sigmaMax=", sigmaMax, "\n")
        trial <- try(uniroot(fcn, interval = c(sigmaLowerLimit, sigmaUpperLimit)), silent = TRUE)
        sigma[i] <- if (inherits(trial, "try-error")) 2 * sigmaMax else trial$root
        #> cat("  sigma[i]=", sigma[i], "\n")
        use <- sigma[i] < sigmaMax
        if (!any(use)) {
            sigma[i] <- sigma[i - 1] # probably good enough; this occurs only at sigma=1, I think
        }
    }
    approxfun(epsilon, sigma)
}


#' Whale compression force
#'
#' Calculate the total compression stress and force, along
#' with the thicknesses of skin, blubber, sublayer, and bone.
#' The stress is computed with the [stressFromStrainFunction()] function that
#' is created by [parameters()] and stored in `para`.
#' the force is computed by multiplying stess by area
#' computed as the product of `parms$Ly` and `parms$Lz`.
#' Any negative layer thicknesses are set to zero, as a way to
#' avoid problems with aphysical engineering compression strains that
#' exceed 1.
#'
#' @param xs Ship position (m).
#'
#' @param xw Whale position (m).
#'
#' @template parmsTemplate
#'
#' @return A list containing: `force` (N), the
#' compression-resisting force; `stress` (Pa), the ratio
#' of that force to the impact area; `strain`, the total
#' strain, and `compressed`, a four-column matrix (m)
#' with first column for skin compression, second for blubber
#' compression, third for sublayer compression, and fourth
#' for bone compression.
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
whaleCompressionForce <- function(xs, xw, parms) {
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    # Note that the denominator of the strain expression vanishes in the stress calculation,
    # so the next three lines could be simplified. However, retaining it might be clearer,
    # if a nonlinear stress-strain relationship becomes desirable in the future.
    strain <- dx / parms$lsum
    stress <- parms$stressFromStrain(strain)
    force <- stress * parms$Ly * parms$Lz
    stress <- ifelse(stress < 0, 0, stress) # just in case; we don't want log(negative number)
    compressed <- cbind(
        parms$l[1] * (1 - log(1 + stress / parms$a[1]) / parms$b[1]),
        parms$l[2] * (1 - log(1 + stress / parms$a[2]) / parms$b[2]),
        parms$l[3] * (1 - log(1 + stress / parms$a[3]) / parms$b[3]),
        parms$l[4] * (1 - log(1 + stress / parms$a[4]) / parms$b[4])
    )
    compressed <- pin(compressed, lower = 0)
    list(force = force, stress = stress, strain = strain, compressed = compressed)
}

#' Skin force
#'
#' The ship-whale separation is used to calculate the deformation of the skin. The
#' parameters of the calculation are `parms$Ly` (impact area width, m),
#' `parms$Lz` (impact area height, in m), `parms$Ealpha` (skin elastic modulus in Pa),
#' `parms$alpha` (skin thickness in m), and `parms$theta` (skin bevel angle
#' degrees, measured from a vector normal to undisturbed skin).
#'
#' @param xs Ship position (m).
#'
#' @param xw Whale position (m).
#'
#' @template parmsTemplate
#'
#' @return A list containing `force`, the normal force (N), along with
#' `sigmay` and `sigmaz`, which are stresses (Pa) in the y (beam)
#' and z (draft) directions.
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
whaleSkinForce <- function(xs, xw, parms) {
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    C <- cos(parms$theta * pi / 180) # NB: theta is in deg
    S <- sin(parms$theta * pi / 180) # NB: theta is in deg
    lambda <- dx * S / C # dek20180622_skin_strain eq 1; called l until 20180725
    Lambda <- dx / C # dek20180622_skin_strain eq 2; called s until 20180725
    # Strains in y and z
    epsilony <- 2 * (Lambda - lambda) / (parms$Ly + 2 * lambda) # dek20180622_skin_strain  eq 3
    epsilonz <- 2 * (Lambda - lambda) / (parms$Lz + 2 * lambda) # analogous to dek20180622 eq 3
    # Stresses in y and z
    sigmay <- parms$a[1] * (exp(parms$b[1] * epsilony) - 1)
    sigmaz <- parms$a[1] * (exp(parms$b[1] * epsilonz) - 1)
    # Net normal force in x; note the cosine, to resolve the force to the normal
    # direction, and the 2, to account for two sides of length
    # Ly and two of length Lz
    force <- 2 * parms$l[1] * (parms$Lz * sigmaz + parms$Ly * sigmay) * C # dek20180622_skin_strain eq 8
    list(force = force, sigmay = sigmay, sigmaz = sigmaz)
}


#' Ship water force
#'
#' Compute the retarding force of water on the ship, based on a drag law
#' \eqn{(1/2)*rho*Cs*A*vs^2}{(1/2)*rho*Cs*A*vs^2}
#' where `rho` is 1024 (kg/m^3), `Cs` is `parms$Cs` and
#' `A` is `parms$Ss`.
#
#' @param vs ship velocity (m/s).
#'
#' @template parmsTemplate
#'
#' @return Water drag force (N).
#'
#' @author Dan Kelley
#'
#' @export
shipWaterForce <- function(vs, parms) {
    -(1 / 2) * 1024 * parms$Cs * parms$Ss * vs * abs(vs)
}


#' Whale force
#'
#' Compute the retarding force of water on the whale, based on a drag law
#' \eqn{(1/2)*rho*Cw*A*vw^2}{(1/2)*rho*Cw*A*vw^2}
#' where `rho` is 1024 (kg/m^3), `Cw` is `parms$Cw` and
#' `A` is `parms$Sw`.
#'
#' @param vw Whale velocity (m/s).
#'
#' @template parmsTemplate
#'
#' @return Water drag force (N).
#'
#' @author Dan Kelley
#'
#' @export
whaleWaterForce <- function(vw, parms) {
    -(1 / 2) * 1024 * parms$Cw * parms$Sw * vw * abs(vw)
}

#' Ship water force
#'
#' Compute the retarding force of water on the ship, based on a drag law
#' \eqn{(1/2)*rho*Cs*A*vs^2}{(1/2)*rho*Cs*A*vs^2}
#' where `rho` is 1024 (kg/m^3), `Cs` is `parms$Cs` and
#' `A` is `parms$Ss`.
#
#' @param vs ship velocity (m/s).
#'
#' @template parmsTemplate
#'
#' @return Water drag force (N).
#'
#' @author Dan Kelley
#'
#' @export
shipWaterForce <- function(vs, parms) {
    -(1 / 2) * 1024 * parms$Cs * parms$Ss * vs * abs(vs)
}
