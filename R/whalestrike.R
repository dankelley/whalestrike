library(deSolve)
library(xtable)

### #' Version of parameter defaults
### #'
### #' A character value that indicates which defaults are to be
### #' used by [parameters()].
### #'
### #' There are two choices, as of 2019 May 28:
### #'
### #'\itemize{
### #' \item \code{"2018"} to get default parameter values from the work
### #' in summer of 2018
### #' \item \code{"2019a"} to get parameter values as of
### #' the start of summer, 2019.
### #'}
### #'
### #' @name versionOfDefaults
### #' @docType data
### NULL
### versionOfDefaults <- "2019a"


#' Convert a speed in knots to a speed in m/s
#'
#' This is done by multiplying by the factor 1.852e3/3600,
#' according to https://en.wikipedia.org/wiki/Knot_(unit).
#' See also [mps2knot()], which is the inverse of this function.
#'
#' @param knot Speed in knots.
#'
#' @return Speed in m/s.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @md
knot2mps <- function(knot)
{
    knot * 1.852e3 / 3600 # exact definition according to https://en.wikipedia.org/wiki/Knot_(unit)
}

#' Convert a speed in m/s to a speed in knots
#'
#' This is done by dividing by the factor 1.852e3/3600,
#' according to https://en.wikipedia.org/wiki/Knot_(unit).
#' See also [knot2mps()], which is the inverse of this function.
#'
#' @param mps Speed in metres per second.
#'
#' @return Speed in knots.
#'
#' @author Dan Kelley
#'
#' @export
#' @md
mps2knot <- function(mps)
{
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
#'
#' @md
pin <- function(x, lower=NULL, upper=NULL)
{
    ## Protect the ifelse() operation from getting riled by NAs
    na <- is.na(x)
    x[na] <- 0 # value is arbitrary because changed back to NA later
    if (!is.null(lower))
        x <- ifelse(x > lower, x, lower)
    if (!is.null(upper))
        x <- ifelse(x < upper, x, upper)
    x[na] <- NA
    x
}

#' Draw polygon between two xy curves
#'
#' This adds to an existing plot by filling the area between the
#' lower=lower(x) and upper=upper(x) curves.  In most cases, as
#' shown in \dQuote{Examples}, it is helpful
#' to use `xaxs="i"` in the preceding plot call, so that the
#' polygon reaches to the edge of the plot area.
#'
#' @param x Coordinate along horizontal axis
#'
#' @param lower Coordinates of the lower curve, of same length as `x`,
#' or a single value that gets repeated to the length of `x`.
#'
#' @param upper Coordinates of the upper curve, or a single value that gets
#' repeated to the length of `x`.
#'
#' @param ... passed to [polygon()]. In most cases, this
#' will contain `col`, the fill colour, and possibly `border`,
#' the border colour, although cross-hatching with `density`
#' and `angle` is also a good choice.
#'
#' @examples
#' ## 1. CO2 record
#' plot(co2, xaxs="i", yaxs="i")
#' fillplot(time(co2), min(co2), co2, col="pink")
#'
#' ## 2. stack (summed y) plot
#' x <- seq(0, 1, 0.01)
#' lower <- x
#' upper <- 0.5 * (1 + sin(2 * pi * x / 0.2))
#' plot(range(x), range(lower, lower+upper), type='n',
#'      xlab="x", ylab="y1, y1+y2",
#'      xaxs="i", yaxs="i")
#' fillplot(x, min(lower), lower, col="darkgray")
#' fillplot(x, lower, lower+upper, col="lightgray")
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom graphics polygon
#'
#' @md
fillplot <- function(x, lower, upper, ...)
{
    n <- length(x)
    if (length(lower) == 1)
        lower <- rep(lower, n)
    if (length(upper) == 1)
        upper <- rep(upper, n)
    if (n != length(lower)) stop("lengths of x and lower must match")
    if (n != length(upper)) stop("lengths of x and upper must match")
    xx <- c(x, rev(x), x[1])
    yy <- c(upper, rev(lower), upper[1])
    polygon(xx, yy, ...)
}


#' Whale blubber stress-strain relationship
#'
#' This is a data frame with elements `strain` and `stress`,
#' found by digitizing (accurate to perhaps 1 percent) the curve shown in Figure 2.13
#' of Raymond (2007). It is used to develop a stress-strain relationship used
#' by [parameters()], as shown in \dQuote{Examples}.
#'
#' @template ref_raymond
#'
#' @examples
#' data(raymond2007)
#' attach(raymond2007)
#' ## Next yields \code{a=1.64e5} Pa and \code{b=2.47}.
#' m <- nls(stress~a*(exp(b*strain)-1), start=list(a=1e5, b=1))
#' plot(strain, stress, xaxs="i", yaxs="i")
#' x <- seq(0, max(strain), length.out=100)
#' lines(x, predict(m, list(strain=x)))
#'
#' @name raymond2007
#'
#' @docType data
#'
#' @md
NULL

#' Whale side-view shape
#'
#' This is a data frame containing 45 points specifying the shape
#' of a right whale, viewed from the side. It was created
#' by digitizing the whale shape (ignoring fins) that is
#' provided in the necropsy reports of Daoust et al. (2018). The
#' data frame contains `x` and `y`, which are distances
#' nondimensionalized by the range in `x`; that is, `x`
#' ranges from 0 to 1. The point at the front of the whale is
#' designated as x=y=0.
#'
#' @template ref_daoust
#'
#' @examples
#' library(whalestrike)
#' data(whaleshape)
#' plot(whaleshape$x, whaleshape$y, asp=1, type="l")
#' polygon(whaleshape$x, whaleshape$y, col="lightgray")
#' lw <- 13.7
#' Rmax <- 0.5 * lw * diff(range(whaleshape$y))
#' mtext(sprintf("Max. radius %.2fm for %.1fm-long whale", Rmax, lw), side=3)
#'
#' @name whaleshape
#'
#' @docType data
#'
#' @md
NULL


#' whalestrike: A Package to Simulate Ship-Whale Collisions
#'
#' This package solves Newton's second law for a simple model of
#' a ship colliding with a whale. This is a stripped-down model
#' that does not attempt to simulate the biomechanical interactions
#' that can be simulated in finite-element treatments such
#' as that of Raymond (2007).  The goal is to establish a
#' convenient framework for rapid computation of impacts in a
#' wide variety of conditions. The model runs
#' quickly enough to keep up with mouse movements to select
#' conditions, in a companion R shiny app. That app lets
#' the user see the effects of changing contact
#' area, ship speed, etc., as a way to build intuition
#' for scenarios ranging from collision with a slow-moving
#' fishing boat to collision with a thin dagger-board of a
#' much swifter racing sailboat. Another advantage of the
#' simple formulation is that it makes it easy to
#' modify various dynamical and biomechanical parameters,
#' to add new forces, and to explore a range of
#' criteria for whale damage.
#'
#' The documentation for [strike()] provides
#' a practical example of using the main functions of this package,
#' while the package vignette provides a general overview.
#' A companion manuscript is intended to
#' provide more detail about the mathematical
#' framework of the package, along with a discussion of its
#' purpose and application to real-world problems of ship
#' strikes on whales.
#'
#' @section Further reading:
#' \itemize{
#'
#' \item
#'
#' Daoust, Pierre-Yves, Émilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' “Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017.” Canadian Wildlife Health Cooperative, Marine Animal Response
#' Socieity, and Fisheries and Oceans Canada, 2018.
#' \url{http://publications.gc.ca/site/eng/9.850838/publication.html}.
#'
#' \item
#' Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. “Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis).” Journal of Mammalogy 93,
#' no. 5 (2012): 1342–54. \url{https://doi.org/10.1644/11-MAMM-A-297.1}.
#'
#' \item
#' Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt, Adam
#' P. Summers, and Petra Ditsche. “Mechanical Properties of Harbor Seal Skin and
#' Blubber--a Test of Anisotropy.” Zoology 126 (2018): 137–44.
#' \url{https://doi.org/10.1016/j.zool.2017.11.002}.
#'
#' \item
#' Kelley, Dan E. “Composite Spring,” May 28, 2018. 20180528_composite_string. Dan Kelley’s working notes.
#'
#' \item
#' Kelley, Dan. “Whale Area,” June 23, 2018. 20180623_whale_area. Dan Kelley’s working notes.
#'
#' \item
#' Kelley, Dan. “Ship Propulsion,” July 1, 2018. 20180701_ship_propulsion. Dan Kelley’s working notes.
#'
#' \item
#' Kelley, Dan. “Whale Mass,” July 7, 2018. 20180707_whale_mass. Dan Kelley’s working notes.
#'
#' \item
#' MAN Diesel & Turbo. “Basic Principles of Propulsion.” MAN Diesel & Turbo, 2011.
#' \url{https://spain.mandieselturbo.com/docs/librariesprovider10/sistemas-propulsivos-marinos/basic-principles-of-ship-propulsion.pdf?sfvrsn=2}
#'
#' \item
#' Manen, J. D. van, and P. van Oossanen. “Resistance.” In Principles of Naval
#' Architecture (Second Revision), Volume II - Resistance, Propulsion and
#' Vibration, edited by Edward V Lewis, 2nd ed., 1–125. Jersey City, NJ: Society
#' of Naval Architects and Marine Engineers (U.S.), 1988.
#'
#' \item
#' Miller, Carolyn A., Desray Reeb, Peter B. Best, Amy R. Knowlton, Moira W.
#' Brown, and Michael J. Moore. “Blubber Thickness in Right Whales Eubalaena
#' Glacialis and Eubalaena Australis Related with Reproduction, Life History
#' Status and Prey Abundance.” Marine Ecology Progress Series 438 (2011): 267–83.
#'
#' \item
#' Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
#' “Morphometry, Gross Morphology and Available Histopathology in North Atlantic
#' Right Whale (Eubalaena Glacialis) Mortalities (1970 to 2002).” Journal of
#' Cetacean Research and Management 6, no. 3 (2005): 199–214.
#'
#' \item
#' Ng, Laurel J., Vladislav Volman, Melissa M. Gibbons, Pi Phohomsiri, Jianxia
#' Cui, Darrell J. Swenson, and James H. Stuhmiller. “A Mechanistic End-to-End
#' Concussion Model That Translates Head Kinematics to Neurologic Injury.”
#' Frontiers in Neurology 8, no. JUN (2017): 1–18.
#' \url{https://doi.org/10.3389/fneur.2017.00269}
#'
#' \item
#' Raymond, J. J. “Development of a Numerical Model to Predict Impact Forces on a
#' North Atlantic Right Whale during Collision with a Vessel.” University of New
#' Hampshire, 2007.
#' \url{https://scholars.unh.edu/thesis/309}.
#'
#'}
#'
#' @docType package
#'
#' @name whalestrike
#'
#' @md
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
#' a piecewise-linear represention of the relationship between \eqn{\sigma} and \eqn{\Delta L},
#' which becomes the return value of the present function.
#' (The purpose of using a the piecewise-linear representation os to shorten
#' computation times.)
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
#' ## Set blubber parameters for each layer, to see if
#' ## we recover the raymond2007 data.
#' param <- parameters(a=rep(1.64e5,4), b=rep(2.47,4))
#' x <- seq(0, 0.5, length.out=100)
#' y <- param$stressFromStrain(x)
#' plot(x, y, type='l', lwd=4, col="gray")
#' data("raymond2007")
#' points(raymond2007$strain, raymond2007$stress, col=2)
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom stats approxfun uniroot
#'
#' @md
stressFromStrainFunction <- function(l, a, b, N=1000)
{
    fcn <- function(sigma) {
        DL - sum((l/bb) * log(1 + sigma / aa))
    }
    L <- sum(l)
    sigma <- rep(NA, N)
    epsilon <- seq(0, 1, length.out=N)
    ## We will limit the epsilon in any given layer to 1, i.e. we don't
    ## permit layers to be compressed to negative thickness. This is expressed
    ## by setting layer-by-layer conditions on maximum stress.
    sigmaMax <- a * (exp(b) - 1)
    ##debug cat("sigmaMax=", paste(sprintf("%.3g", sigmaMax), collapse=" "), "\n")
    use <- rep(TRUE, length(a))
    for (i in seq_along(epsilon)) {
        aa <- a[use]
        bb <- b[use]
        ##debug cat(sprintf("i=%3d, epsilon=%10.5f, ", i, epsilon[i]), ", use=", paste(use, collapse=" "), "\n")
        DL <- epsilon[i] * L
        sigma[i] <- uniroot(fcn, interval=c(0, 1e9))$root
        ##debug cat("  sigma[i]=", sprintf("%.3g", sigma[i]), "\n")
        use <- sigma[i] < sigmaMax
        if (!any(use))
            sigma[i] <- sigma[i-1] # probably good enough; this occurs only at sigma=1, I think
    }
    approxfun(epsilon, sigma)
}


#' Control parameters for whale strike simulation
#'
#' Assembles control parameters into a list suitable for passing to [strike()]
#' and the functions that it calls. If `file` is provided, then all the other
#' arguments are read from that source.
#' Below are some sources cited in the discussion of the function arguments.
#' @template ref_campbell_malone
#' @template ref_daoust
#' @template ref_grear
#' @template ref_raymond
#'
#' @param ms Ship mass (kg).
#'
#' @param Ss Ship wetted area (m^2). This, together with `Cs`, is used by
#' used by [shipWaterForce()] to estimate ship drag force. If `Ss`
#' is not given, then an esimate is made by calling [shipAreaFromMass()] with
#' the provided value of `ms`.
#'
#' @param Ly Ship impact horizontal extent (m); defaults to 1.23m if not specified,
#' based on an analysis of the shape of the bow of typical coastal fishing boats
#' of the Cape Islander variety.
#'
#' @param Lz Ship impact vertical extent (m); defaults to 1.23m if not specified,
#' based on the same analysis as for Ly.
#'
#' @param lw Whale length (m). This is used by [whaleAreaFromLength()] to
#' calculate area, which is needed for the water drag calculation done by
#' [whaleWaterForce()].
#'
#' @param species A string indicating the whale species. For the permitted values,
#' see [whaleMassFromLength()].
#'
#' @param mw Whale mass (kg). If this value is not provided, then
#' it is calculated from whale length, using [whaleMassFromLength()]
#' with `type="wetted"`.
#'
#' @param Sw Whale surface area (m^2). If not provided, this is calculated
#' from whale length using [whaleAreaFromLength()].
#'
#' @param l Numerical vector of length 4, giving thickness (m) of skin, blubber,
#' sublayer, and bone. If not provided, this is set to
#' `c(0.025, 0.16, 1.12, 0.1)`.
#' The skin thicknes default of 0.025 m represents the 0.9-1.0 inch value
#' stated in Section 2.2.3 of Raymond (2007).
#' The blubber default of 0.16 m is a rounded average of the values inferred
#' by whale necropsy, reported in Appendix 2 of Daoust et al., 2018.
## > round(mean(c(17,14,18.13,18,21.25,16.75,13.33,7)/100),2)
## [1] 0.16
#' The sublayer default of 1.12 m may be reasonable at some spots on the whale body.
#' The bone default of 0.1 m may be reasonable at some spots on the whale body.
#' The sum of these default values, 1.40 m, is a whale radius that
#' is consistent with a half-circumferance of 4.4 m, reported in Table 2.2
#' of Raymond (2007).
#'
#' @param a,b Numerical vectors of length 4, giving values to use in the
#' stress-strain law `stress=a*(exp(b*strain)-1)`, where `a` is in Pa
#' and `b` is unitless. By construction, `a*b` is the local modulus at
#' low strain (i.e. at low `b*strain` values), and that `b` is the
#' efolding scale for nonlinear increase in stress with strain.
#' This exponential relationship has been mapped out
#' for whale blubber, using a curve fit to Figure 2.13 of Raymond (2007), and
#' these values are used for the second layer (blubber); see
#' the documentation for the [raymond2007] dataset, to see
#' for how that fit was done.
#' If not provided, `a` defaults to
#' `c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)`
#' and `b` defaults to
#' `c(0.1, 2.54, 2.54, 0.1)`.
#' The skin defaults are set up to give a linear shape (since `b` is small)
#' with the `a*b` product
#' being 17.8e6 Pa, which is the adult-seal value
#' given in Table 3 of Grear et al. (2017).
#' The blubber defaults are from a regression of the stress-strain
#' relationship shown in Figure 2.13 of Raymond (2007).
#' The sublayer defaults are set to match those of blubber, lacking
#' any other information.
#' The bone default for `b` is small, to set up a linear function,
#' and `a*b` is set to equal 8.54e8 Pa,
#' given in Table 2.3 of Raymond (2007) and Table 4.5 of
#' Campbell-Malone (2007).
#'
#' @param s Numerical vector of length 4, giving the ultimate strengths (Pa) of
#' skin, blubber, sublayer, and bone, respectively. If not provided, the
#' value is set to
#' `c(19.6e6, 0.437e6, 0.437e6, 22.9e6)`,
#' with reasoning as follows.
#' The skin default of 19.6 MPa
#' is a rounded value from Table 3 of Grear et al. (2018) for adult seal skin strength at
#' an orientation of 0 degrees.  The blubber value of
#' 0.437 MPa is inferred by
#' multiplying Raymond's (2007) Figure 2.13 elastic modulus of 0.636 MPa
#' by the ratio 0.97/1.41 determined for adult seal strength/modulus, as reported
#' in Table 3 of Grear et al. (2018).
#' The sublayer value is taken to be identical to the blubber value, lacking
#' more specific information.
#' The bone default o 22.9 MPa is from Table 2.3 of Raymond (2007) and
#' Table 4.5 of Campbell-Malone (2007).
#'
#' @param theta Whale skin deformation angle (deg); defaults to 55 degrees,
#' if not supplied, because that angle produces a good match to Raymond's (2007)
#' Figure 6.1 for the total force as a function of vessel speed, for large
#' vessels. (Note that the match works almost as well in the range 50 deg
#' to 70 deg.)
#'
#' @param Cs Drag coefficient for ship (dimensionless),
#' used by [shipWaterForce()] to estimate ship drag force. Defaults
#' to 1e-2, which is 4 times the frictional coefficient of 2.5e-3
#' inferred from Figure 4 of Manen and van Oossanen (1988), assuming
#' a Reynolds number of 5e7, computed from speed 5m/s, lengthscale 10m
#' and viscosity 1e-6 m^2/s. (The factor of 4 is under the assumption
#' that frictional drag is about a quarter of total drag.)
#' The drag force is computed with [shipWaterForce()].
#'
#' @param Cw Drag coefficient for whale (dimensionless),
#' used by [whaleWaterForce()] to estimate whale drag force.
#' Defaults to 2.5e-3, for Reynolds number 2e7, computed from speed
#' 2 m/s, lengthscale 5m (between radius and length) and
#' viscosity 1e-6 m^2/s.  The drag force is computed with
#' [whaleWaterForce()].
#'
#' @param file Optional name a comma-separated file that holds all of the
#' previous values, except `Cs` and `Cw`. If provided,
#' then other parameters (except `Cs` and `Cw`) are
#' ignored, because values are sought from the file. The purpose of
#' this is in shiny apps that want to save a simulation framework.
#' The file should be saved [write.csv()] with
#' `row.names=FALSE`.
#'
#' @return
#' A named list holding the parameters, with defaults and alternatives reconciled
#' according to the system described above, along with a function that computes
#' compression force, which created by [stressFromStrainFunction()].
#'
#' @examples
#' parms <- parameters()
#' epsilon <- seq(0, 1, length.out=100) # strain
#' sigma <- parms$stressFromStrain(epsilon) # stress
#' plot(epsilon, log10(sigma), xlab="Strain", ylab="log10(Stress [MPa])", type="l")
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom utils read.csv
#'
#' @md
parameters <- function(ms=45e3, Ss, Ly=1.15, Lz=1.15,
                       species="N. Atl. Right Whale",
                       lw=13.7, mw, Sw,
                       l, a, b, s,
                       theta=55,
                       Cs=0.01, Cw=0.0025, file)
{
    if (!missing(file)) {
        rval <- as.list(read.csv(file))
        rval$Ss <- shipAreaFromMass(rval$ms)
        rval$mw <- whaleMassFromLength(rval$lw, species=species)
        rval$Sw <- whaleAreaFromLength(rval$lw, species=species, type="wetted")
        rval$tmax <- NULL
        rval$vs <- NULL
        rval$Cs <- Cs
        rval$Cw <- Cw
        rval$l <- c(rval$l1, rval$l2, rval$l3, rval$l4)
        rval$l1 <- rval$l2 <- rval$l3 <- rval$l4 <- NULL
        rval$lsum <- sum(rval$l)
        ## the next are copied from below. The app doesn't let the user
        ## set these things, so we know their values.
        ## NOTE: keep in synch with 'BBBB' below!
        rval$a <- c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)
        rval$b <- c(0.1, 2.54, 2.54, 0.1)
        rval$s <- c(19.6e6, 0.437e6, 0.437e6, 22.9e6)
        o <- sort(names(rval))
        rval <- rval[o]
    } else {
        if (length(ms) != 1) stop("cannot handle more than one 'ms' at a time")
        if (ms <= 0)
            stop("ship mass (ms) must be a positive number, but it is ", ms)
        if (missing(Ss))
            Ss <- shipAreaFromMass(ms)
        if (length(Ss) != 1) stop("cannot handle more than one 'Ss' at a time")
        if (Ss <= 0)
            stop("ship wetted area (Ss) must be a positive number, but it is ", Ss)
        if (length(Ly) != 1) stop("cannot handle more than one 'Ly' at a time")
        if (Ly <= 0)
            stop("impact width (Ly) must be a positive number, but it is ", Ly)
        if (length(Lz) != 1) stop("cannot handle more than one 'Lz' at a time")
        if (Lz <= 0)
            stop("impact height (Lz) must be a positive number, but it is ", Lz)
        if (length(lw) != 1) stop("cannot handle more than one 'lw' at a time")
        if (lw <= 0)
            stop("Whale length (lw) must be a positive number")
        if (missing(mw))
            mw <- whaleMassFromLength(lw, species=species)
        if (length(mw) != 1) stop("cannot handle more than one 'mw' at a time")
        if (missing(Sw))
            Sw <- whaleAreaFromLength(lw, species=species, type="wetted")
        if (length(Sw) != 1) stop("cannot handle more than one 'Sw' at a time")
        if (missing(l))
            l <- c(0.025, 0.16, 1.12, 0.1)
        if (length(l) != 4) stop("'l' must be a vector of length 4")
        ## NOTE: keep in synch with 'AAAA' above!
        if (missing(a))
            a <- c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)
        if (length(a) != 4) stop("'a' must be a vector of length 4")
        ## NOTE: keep in synch with above!
        if (missing(b))
            b <- c(0.1, 2.54, 2.54, 0.1)
        if (length(b) != 4) stop("'b' must be a vector of length 4")
        ## NOTE: keep in synch with above!
        if (missing(s))
            s <- c(19.6e6, 0.437e6, 0.437e6, 22.9e6)
        if (length(s) != 4) stop("'s' must be a vector of length 4")
        ## Value checks
        if (any(l <= 0) || length(l) != 4)
            stop("'l' must be 4 positive numbers")
        if (any(a <= 0) || length(a) != 4)
            stop("'a' must be 4 positive numbers")
        if (any(b <= 0) || length(b) != 4)
            stop("'b' must be 4 positive numbers")
        if (length(theta) != 1) stop("cannot handle more than one 'theta' at a time")
        if (theta < 0 || theta > 89)
            stop("whale skin deformation angle (theta) must be between 0 and 89 deg, but it is ", theta)
        if (length(Cs) != 1) stop("cannot handle more than one 'Cs' at a time")
        if (Cs <= 0)
            stop("ship resistance parameter (Cs) must be positive, but it is ", Cs)
        if (length(Cw) != 1) stop("cannot handle more than one 'Cw' at a time")
        if (Cw <= 0)
            stop("ship resistance parameter (Cw) must be positive, but it is ", Cw)
        rval <- list(ms=ms, Ss=Ss,
                     Ly=Ly, Lz=Lz,
                     mw=mw, Sw=Sw, lw=lw,
                     l=l, lsum=sum(l), a=a, b=b, s=s,
                     theta=theta,
                     Cs=Cs, Cw=Cw)
    }
    ## For efficiency, create and store an overall stress-strain function
    rval$stressFromStrain <- stressFromStrainFunction(rval$l, rval$a, rval$b)
    class(rval) <- "parameters"
    rval
}


#' Whale mass inferred from length
#'
#' Calculate an estimate of the mass of different species of whale,
#' based on animal length, based on formulae as listed in
#' \dQuote{Details}.
#'
#' The permitted values for `model` are as follows.
#'
#' * `"moore2005"` (which only works if `species` is `"N. Atl. Right Whale"`) yields
#' \eqn{242.988 * exp(0.4 * length)}{242.988 * exp(0.4 * L)},
#' which (apart from a unit change on `L`) is the regression equation
#' shown above Figure 1d in Moore et al. (2005) for right whales. A
#' difficult in the Moore et al. (2005) use of a single nonzero digit
#' in the multiplier on `L` is illustrated in \dQuote{Examples}.
#'
#' * `"fortune2012"` with `species="N. Atl. Right Whale"` yields the formula
#' \eqn{exp(-10.095 + 2.825*log(100*L))}{exp(-10.095 + 2.825*log(100*L))}
#' for North Atlantic right whales, according to corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'
#' * `"fortune2012"` with `species="N. Pac. Right Whale"` yields the formula
#' \eqn{exp(-12.286 + 3.158*log(100*L))}{exp(-12.286 + 3.158*log(100*L))}
#' for North Pacific right whales, according to corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'
#' * `"lockyer1976"` uses formulae from table 1 of Lockyer (1976). The
#' permitted `species` and the formulae used are as follows.
#' * `"Pac. Right Whale": \eqn{1e3 * 0.013200 * L^3.06}{1e3 * 0.013200 * L^3.06}
#' * `"Blue Whale"`: \eqn{1e3 * 0.002899 * L^3.25}{1e3 * 0.002899 * L^3.25}
#' * `"Fin Whale"`: \eqn{1e3 * 0.007996 * L^2.90}{1e3 * 0.007996 * L^2.90}
#' * `"Sei Whale"`: \eqn{1e3 * 0.025763 * L^2.43}{1e3 * 0.025763 * L^2.43}
#' * `"Bryde Whale"`: \eqn{1e3 * 0.012965 * L^2.74}{1e3 * 0.012965 * L^2.74}
#' * `"Minke Whale"`: \eqn{1e3 * 0.049574 * L^2.31}{1e3 * 0.049574 * L^2.31}
#' * `"Humpback Whale"`: \eqn{1e3 * 0.016473 * L^2.95}{1e3 * 0.016473 * L^2.95}
#' * `"Sperm Whale"`: \eqn{1e3 * 0.006648 * L^3.18}{1e3 * 0.006648 * L^3.18}
#'
#' @param L Whale length in m.
#'
#' @param species character value specifying the species; see \dQuote{Details}.
#'
#' @param model character value specifying the model (see \dQuote{Details}).
#'
#' @return Mass in kg.
#'
#' @examples
#' library(whalestrike)
#' L <- seq(5, 15, length.out=100)
#' kpt <- 1000 # kg per tonne
#' plot(L, whaleMassFromLength(L, model="moore2005")/kpt, type='l',
#'      xlab="Right-whale Length [m]", ylab="Mass [tonne]")
#' # Demonstrate sensitivity involved in the single-digit parameter in Moore's formula
#' lines(L, 242.988 * exp(0.35 * L)/kpt, lty='dotted')
#' lines(L, 242.988 * exp(0.45 * L)/kpt, lty='dashed')
#' lines(L, whaleMassFromLength(L, model="fortune2012atlantic")/kpt, col=2)
#' lines(L, whaleMassFromLength(L, model="fortune2012pacific")/kpt, col=3)
#' legend("topleft", lwd=1, col=1:3,
#'        legend=c("moore2005", "fortune2012atlantic", "fortune2012pacific"))
#'
#' @references
#'
#' * Lockyer, C. “Body Weights of Some Species of Large Whales.” J. Cons. Int.
#' Explor. Mer. 36, no. 3 (1976): 259–73.
#'
#' * Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
#' “Morphometry, Gross Morphology and Available Histopathology in North Atlantic
#' Right Whale (Eubalaena Glacialis) Mortalities (1970 to 2002).” Journal of
#' Cetacean Research and Management 6, no. 3 (2005): 199–214.
#'
#' * Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. “Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis).” Journal of Mammalogy 93,
#' no. 5 (2012): 1342–54. https://doi.org/10.1644/11-MAMM-A-297.1.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @md
whaleMassFromLength <- function(L, species="N. Atl. Right Whale", model="fortune2012")
{
    if (model == "moore2005") {
        if (species == "N. Atl. Right Whale")
            242.988 * exp(0.4 * L)
        else
            stop("The 'moore2005' model only works if species is 'N. Atl. Right Whale'")
    } else if (model == "fortune2012") {
        if (species == "N. Atl. Right Whale")
            exp(-10.095 + 2.825*log(100*L))
        else if (species == "Pac. Right Whale")
            exp(-12.286 + 3.158*log(100*L))
        else
            stop("The 'moore2005' model only works if species is 'N. Atl. Right Whale' or 'N. Pac. Right Whale'")
    } else if (model == "lockyer1976") {
        if (species == "Pac. Right Whale")
            1e3 * 0.013200 * L^3.06
        else if (species == "Blue Whale")
            1e3 * 0.002899 * L^3.25
        else if (species == "Fin Whale")
            1e3 * 0.007996 * L^2.90
        else if (species == "Sei Whale")
            1e3 * 0.025763 * L^2.43
        else if (species == "Bryde Whale")
            1e3 * 0.012965 * L^2.74
        else if (species == "Minke Whale")
            1e3 * 0.049574 * L^2.31
        else if (species == "Humpback Whale")
            1e3 * 0.016473 * L^2.95
        else if (species == "Sperm Whale")
            1e3 * 0.006648 * L^3.18
        else
            stop('The "lockyer1976" model requires the species to be one of the following: `"Blue Whale"`, `"Fin Whale"`, `"Sei Whale"`, `"Bryde Whale"`, `"Minke Whale"`, or `"Humpback Whale"`')
    } else {
        stop('model "', model, '" is unknown. This must be one of: "moore2005", "fortune2012", or "lockyer1976"')
    }
}

#' Compute whale length from mass
#'
#' This works by inverting [whaleMassFromLength()] using [uniroot()].
#'
#' @param M Whale mass (kg).
#'
#' @param species A string indicating the whale species. For the permitted values,
#' see [whaleMassFromLength()].
#'
#' @param model Character string specifying the model, with permitted
#' values `"moore2005"`, `"fortune2012atlantic"` and
#" `"fortune2012pacific"`. See the documentation
#' for [whaleMassFromLength()] for the details of these
#' formulations.
#'
#' @return Whale length (m).
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom stats approxfun
#'
#' @md
whaleLengthFromMass <- function(M, species="N. Atl. Right Whale", model="fortune2012atlantic")
{
    rval <- rep(NA, length(M))
    for (i in seq_along(M))
        rval[i] <- uniroot(function(x)
                           M[i] - whaleMassFromLength(x, species=species, model=model),
                           c(0.1, 100))$root
    rval
}

#' Whale projected area, as function of length
#'
#' This depends on calculations based on the digitized shape of
#' a whale necropsy, which is provided as [whaleshape].
#' The results are
#' \preformatted{0.143 * L^2} for the projected area (see reference 1)
#' and
#' \preformatted{0.448 * (0.877 * L)^2} for the wetted area,
#' where the latter uses a correction related to whale mass (see reference 2).
#'
#' Note that multiple digitizations were done,
#' and that the coefficients used in the formulae
#' agreed to under 0.7 percent percent between these digitizations.
#'
#' @param L length in m.
#'
#' @param species A string indicating the whale species. For the permitted values,
#' see [whaleMassFromLength()].
#'
#' @param type character string indicating the type of area, with
#' `"projected"` for a side-projected area, and
#' `"wetted"` for the total wetted area. The wetted
#' area was computed by mathematically spinning a spline fit to the
#' side-view. In both cases, the original data source is the
#' necropsy side-view presented in Daoust et al. (2018).
#'
#' @references
#' 1. Dan Kelley's internal document `dek/20180623_whale_area.Rmd`, available
#' upon request.
#'
#' 2. Dan Kelley's internal document `dek/20180707_whale_mass.Rmd`, available
#' upon request.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @md
whaleAreaFromLength <- function(L, species="N. Atl. Right Whale", type="wetted")
{
    speciesAllowed <- c("N. Atl. Right Whale")
    if (!(species %in% speciesAllowed))
        stop("unknown species \"", species, "\"; use one of the following: \"",
             paste(speciesAllowed, collapse="\", \""), "\"")
    ## below from dek/20180623_whale_area.Rmd, updated 20180802 and inserted with
    ## cut/paste (changing bullet to asterisk, and using ^ for exponentiation).
    ##
    ## * Projected area, with fins: 0.1466 ∗ L^2 where L is body length in metres.
    ## * Projected area, without fins: 0.1391 ∗ L^2 where L is body length in metres.
    ## * Wetted area, with fins: 0.4631 ∗ L^2 where L is body length in metres.
    ## * Wetted area, without fins: 0.4336 ∗ L^2 where L is body length in metres.
    ##
    ## The relevant case (with or without fins) being dependent on the application,
    ## there may be merit in averaging the two estimates, yielding:
    ## * Projected area: 0.1429 ∗ L^2 where L is body length in metres.
    ## * Wetted area: 0.4484 ∗ L^2 where L is body length in metres.
    if (type == "projected")
        0.143 * L^2
    else if (type == "wetted")
        0.448 * (0.877 * L)^2
    else stop("'type' must be 'projected' or 'wetted', but it is '", type, "'")
}

#' Whale compression force
#'
#' Calculate the total compression stress and force, along
#' with the thicknesses of skin, blubber, sublayer, and bone.
#' The stess is computed with the [stressFromStrainFunction()] function that
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
#'
#' @md
whaleCompressionForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    ## Note that the denominator of the strain expression vanishes in the stress calculation,
    ## so the next three lines could be simplified. However, retaining it might be clearer,
    ## if a nonlinear stress-strain relationship becomes desirable in the future.
    strain <- dx / parms$lsum
    ##E <- (parms$alpha + parms$beta + parms$gamma) / (parms$alpha/parms$Ealpha + parms$beta/parms$Ebeta + parms$gamma/parms$Egamma)
    stress <- parms$stressFromStrain(strain)
    force <- stress * parms$Ly * parms$Lz
    stress <- ifelse(stress < 0, 0, stress) # just in case; we don't want log(negative number)
    compressed <- cbind(parms$l[1]*(1-log(1 + stress / parms$a[1]) / parms$b[1]),
                        parms$l[2]*(1-log(1 + stress / parms$a[2]) / parms$b[2]),
                        parms$l[3]*(1-log(1 + stress / parms$a[3]) / parms$b[3]),
                        parms$l[4]*(1-log(1 + stress / parms$a[4]) / parms$b[4]))
    ##. message("compressed=", paste(compressed, " "), "before zero trim")
    compressed <- pin(compressed, lower=0)
    ##. message("compressed=", paste(compressed, " "), "after zero trim")
    list(force=force, stress=stress, strain=strain, compressed=compressed)
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
#'
#' @md
whaleSkinForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    C <- cos(parms$theta * pi / 180) # NB: theta is in deg
    S <- sin(parms$theta * pi / 180) # NB: theta is in deg
    lambda <- dx * S / C               # dek20180622_skin_strain eq 1; called l until 20180725
    Lambda <- dx / C                   # dek20180622_skin_strain eq 2; called s until 20180725
    ## Strains in y and z
    epsilony <- 2 * (Lambda - lambda) / (parms$Ly + 2 * lambda) # dek20180622_skin_strain  eq 3
    epsilonz <- 2 * (Lambda - lambda) / (parms$Lz + 2 * lambda) # analogous to dek20180622 eq 3
    ## Stresses in y and z
    sigmay <- parms$a[1] * (exp(parms$b[1] * epsilony) - 1)
    sigmaz <- parms$a[1] * (exp(parms$b[1] * epsilonz) - 1)
    ## Net normal force in x; note the cosine, to resolve the force to the normal
    ## direction, and the 2, to account for two sides of length
    ## Ly and two of length Lz
    F <- 2*parms$l[1]*(parms$Lz*sigmaz + parms$Ly*sigmay)*C # dek20180622_skin_strain eq 8
    list(force=F, sigmay=sigmay, sigmaz=sigmaz)
}


#' Compute ship wetted area from mass
#'
#' The method is based on scaling up the results for a single Cape
#' Islander ship, of displacement 20.46 tonnes, length 11.73m,
#' beam 4.63m, and draft 1.58m, on the assumption that the wetted area
#' is length*(2*draft+beam). This reference area is scaled to
#' the specified mass, `ms`, by multiplying by the 2/3
#' power of mass ratio.
#' This is a crude calculation meant as a stop-gap measure, for
#' estimates values of the `Ss` argument to [parameters()].
#' It would be much preferable, for a particular simulation, to use the
#' wetted area for a particular ship.
#'
#' @param ms Ship mass (kg).
#'
#' @return Estimated area (m^2).
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @md
shipAreaFromMass <- function(ms)
{
    length <- 11.73                        # m
    beam <- 4.63                           # m
    draft <- 1.58                          # m
    displacement <- 20.46e3                # m^3
    factor <- (ms / displacement)^(1/3) # lengthscale factor
    length * (beam + 2 * draft) * factor^2
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
#'
#' @md
shipWaterForce <- function(vs, parms)
{
    - (1/2) * 1024 * parms$Cs * parms$Ss * vs * abs(vs)
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
#'
#' @md
whaleWaterForce <- function(vw, parms)
{
    - (1/2) * 1024 * parms$Cw * parms$Sw * vw * abs(vw)
}

#' Dynamical law
#'
#' @param t time (s).
#'
#' @param y model state, a vector containing ship position `xs` (m),
#' ship speed `vs` (m/s), whale position `xw` (m),
#' and whale speed `vw` (m/s).
#'
#' @template parmsTemplate
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @md
dynamics <- function(t, y, parms)
{
    xs <- y[1]                         # ship position
    vs <- y[2]                         # ship velocity
    xw <- y[3]                         # whale position
    vw <- y[4]                         # whale velocity
    Fcompression <- whaleCompressionForce(xs, xw, parms)$force
    Fextension <- whaleSkinForce(xs, xw, parms)$force
    Freactive <- Fcompression + Fextension
    Fship <- parms$engineForce + shipWaterForce(vs, parms) - Freactive
    ##. if ((t > 0.1 && t < 0.11) || (t > 0.5 && t < 0.51))
    ##.     cat("t=", t, " vs=", vs, " shipEngineForce=", parms$shipEngineForce, " shipWaterForce=", shipWaterForce(vs, parms), " Freactive=", Freactive, "\n")
    if (is.na(Fship[1]))
        stop("Fship[1] is NA, probably indicating a programming error.")
    Fwhale <- Freactive + whaleWaterForce(vw, parms)
    if (is.na(Fwhale[1])) stop("Fwhale[1] is NA, probably indicating a programming error.")
    list(c(dxsdt=vs, dvsdt=Fship/parms$ms, dxwdt=vw, dvwdt=Fwhale/parms$mw))
}

#' Calculate derivative using first difference
#'
#' @param var variable.
#'
#' @param t time in seconds.
#'
#' @return Derivative estimated by using [diff()] on both `x`
#' and `t`. In order to return a value of the same length as `x` and
#' `t`, the last value is repeated.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom utils tail
#'
#' @md
derivative <- function(var, t)
{
    res <- diff(var) / diff(t)
    c(res, tail(res, 1))
}

#' Simulate the collision of a ship and a whale
#'
#' Newtonian mechanics are used, taking the ship as undeformable,
#' and the whale as being cushioned by a skin layer and a blubber layer.
#' The forces are calculated by
#' [shipWaterForce()],
#' [whaleSkinForce()],
#' [whaleCompressionForce()], and
#' [whaleWaterForce()].
#'
#' @param t a suggested vector of times (s) at which the simulated state will be reported.
#' This is only a suggestion, however, because `strike` is set up to detect high
#' accelerations caused by bone compression, and may set a finer reporting interval,
#' if such accelerations are detected. The detection is based on thickness of
#' compressed blubber and sublayer; if either gets to zero thickness, then
#' a new time grid is constructed, with 10 points during the timescale for
#' bone compression, which is assumed to be
#' \eqn{2*sqrt(Ly*Lz*a[4]*b[4]/(l[4]*mw)}, with terms as discussed in
#' the documentation for [parameters()]. If this grid is finer
#' than the grid in the stated `t`, then the simulatoin is redone
#' using the new grid.
#'
#' @param state A list or named vector holding the initial state of the model:
#' ship position `xs` (m),
#' ship speed `vs` (m/s),
#' whale position `xw` (m),
#' and whale speed `vw` (m/s).
#'
#' @template parmsTemplate
#'
#' @template debugTemplate
#'
#' @return An object of class `"strike"`, consisting of a
#' list containing vectors for time (`t` (s)), ship position (`xs` (m)),
#' boat speed (`vs` (m/s)), whale position (`xw` (m)), whale speed (`vw` (m/s)),
#' boat acceleration (`dvsdt` (m/s^2)), and whale acceleration (`dvwdt` (m/s^2)),
#' a list containing the model parameters (`parms`), a list with the results of
#' the skin force calculation (`SWF`), a list with the results of the compression
#' force calculations (`WCF`), and a vector of whale water force (`WWF`).
#'
#' @examples
#' library(whalestrike)
#' # Example 1: graphs, as in the shiny app
#' t <- seq(0, 0.7, length.out=200)
#' state <- list(xs=-2, vs=knot2mps(10), xw=0, vw=0) # ship speed 10 knots
#' parms <- parameters()
#' sol <- strike(t, state, parms)
#' par(mfcol=c(1, 3), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(sol)
#'
#' # Example 2: time-series plots of blubber stress and stress/strength,
#' # for a 200 tonne ship moving at 10 knots
#' t <- seq(0, 0.7, length.out=1000)
#' state <- list(xs=-2, vs=knot2mps(10), xw=0, vw=0) # ship speed 10 knots
#' parms <- parameters(ms=200 * 1000) # 1 metric tonne is 1000 kg
#' sol <- strike(t, state, parms)
#' par(mfrow=c(2, 1), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(t, sol$WCF$stress / 1e6, type="l", xlab="Time [s]", ylab="Blubber stress [MPa]")
#' plot(t, sol$WCF$stress/sol$parms$s[2], type="l", xlab="Time [s]", ylab="Blubber stress / strength")
#'
#' # Example 3: max stress and stress/strength, for a 200 tonne ship moving at various speeds
#' # This is a slow calculation, so we do not run it by default
#' \dontrun{
#' knots <- seq(0, 20, 0.5)
#' maxStress <- NULL
#' maxStressOverStrength <- NULL
#' for (speed in knot2mps(knots)) {
#'     t <- seq(0, 10, length.out=1000)
#'     state <- list(xs=-2, vs=speed, xw=0, vw=0)
#'     parms <- parameters(ms=200 * 1000) # 1 metric tonne is 1000 kg
#'     sol <- strike(t, state, parms)
#'     maxStress <- c(maxStress, max(sol$WCF$stress))
#'     maxStressOverStrength <- c(maxStressOverStrength, max(sol$WCF$stress/sol$parms$s[2]))
#' }
#' par(mfrow=c(2, 1), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
#' nonzero <- maxStress > 0
#' plot(knots[nonzero], log10(maxStress[nonzero]), type="o", pch=20, xaxs="i", yaxs="i",
#'      xlab="Ship Speed [knots]", ylab="log10 peak blubber stress")
#' abline(h=log10(sol$parms$s[2]), lty=2)
#' plot(knots[nonzero], log10(maxStressOverStrength[nonzero]), type="o", pch=20, xaxs="i", yaxs="i",
#'      xlab="Ship Speed [knots]", ylab="log10 peak blubber stress / strength")
#' abline(h=0, lty=2)
#'}
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom deSolve lsoda
#'
#' @md
strike <- function(t, state, parms, debug=0)
{
    if (missing(t))
        stop("must supply t")
    ## Ensure that the state is well-configured, because the error messages
    ## otherwise will be too cryptic for many users to fathom.
    if (missing(state))
        stop("must supply state")
    if (4 != sum(c("xs", "vs", "xw", "vw") %in% names(state)))
        stop("state must hold \"xs\", \"vs\", \"xw\", and \"vw\"")
    if (is.list(state))
        state <- c(xs=state$xs, vs=state$vs, xw=state$xw, vw=state$vw)
    if (missing(parms))
        stop("must supply parms")
    if (!inherits(parms, "parameters"))
        stop("parms must be the output of parameters()")
    if (debug > 0) {
        print("state:")
        print(state)
        print("parms:")
        print(parms)
    }
    for (need in c("xs", "vs", "xw", "vw")) {
        if (!(need %in% names(state)))
            stop("state must contain item named '", need, "'; the names you supplied were: ", paste(names(state), collapse=" "))
    }
    parms["engineForce"] <- -shipWaterForce(state["vs"], parms) # assumed constant over time
    sol <- lsoda(state, t, dynamics, parms)
    ## Add extra things for plotting convenience.
    t <- sol[, 1]
    xs <- sol[, 2]
    vs <- sol[, 3]
    xw <- sol[, 4]
    vw <- sol[, 5]
    dvsdt <- derivative(vs, t)
    dvwdt <- derivative(vw, t)
    SWF <- shipWaterForce(vs=vs, parms=parms)
    WSF <- whaleSkinForce(xs=xs, xw=xw, parms=parms)
    WCF <- whaleCompressionForce(xs=xs, xw=xw, parms=parms)
    WWF <- whaleWaterForce(vw=vw, parms=parms)
    ##. message("max stress ", max(abs(WCF$stress))/1e6, " MPA")
    ##. message("blubber compression ratio=", min(WCF$compressed[,2])/WCF$compressed[1,2])
    ##. message("sublayer compression ratio=", min(WCF$compressed[,3])/WCF$compressed[1,3])
    refinedGrid <- min(WCF$compressed[,2])/WCF$compressed[1,2] < 0.01 || min(WCF$compressed[,3])/WCF$compressed[1,3] < 0.01
    if (refinedGrid) {
        NEED <- 10                     # desired number of points in peak
        dt <- (1/NEED) * 0.5 * sqrt(parms$l[4] * parms$mw / (parms$Ly*parms$Lz*parms$a[4]*parms$b[4]))
        tstart <- t[1]
        tend <- tail(t, 1)
        nold <- length(t)
        n <- floor(0.5 * (tend - tstart) / dt)
        if (n > length(t)) {
            ##message("strike() : increasing from ", nold, " to ", n, " time steps, to capture acceleration peak")
            warning("increasing from ", nold, " to ", n, " time steps, to capture acceleration peak\n")
            t <- seq(tstart, tend, length.out=n)
            sol <- lsoda(state, t, dynamics, parms)
            ## Add extra things for plotting convenience.
            t <- sol[, 1]
            xs <- sol[, 2]
            vs <- sol[, 3]
            xw <- sol[, 4]
            vw <- sol[, 5]
            dvsdt <- derivative(vs, t)
            dvwdt <- derivative(vw, t)
            SWF <- shipWaterForce(vs=vs, parms=parms)
            WSF <- whaleSkinForce(xs=xs, xw=xw, parms=parms)
            WCF <- whaleCompressionForce(xs=xs, xw=xw, parms=parms)
            WWF <- whaleWaterForce(vw=vw, parms=parms)
        }
    }
    res <- list(t=t, xs=xs, vs=vs, xw=xw, vw=vw,
                dvsdt=dvsdt, dvwdt=dvwdt,
                SWF=SWF, WSF=WSF, WCF=WCF, WWF=WWF,
                parms=parms,
                refinedGrid=refinedGrid)
    class(res) <- "strike"
    res
}

#' Plot a strike object
#'
#' Creates displays of various results of a simulation performed
#' with [strike()].
#'
#' @param x An object created by [strike()].
#'
#' @param which A character vector that indicates what to plot.
#' This choices for its entries are listed below, in no particular order.
#' \itemize{
#'
#' \item `"location"` for a time-series plot of boat location `xw` in
#' dashed black, whale centerline `xs` in solid gray,
#' blubber-interior interface in red, and skin in blue. The maximum
#' acceleration of ship and whale (in "g" units) are indicated in notes
#' placed near the horizontal axes. Those acceleration indications report
#' just a single value for each of ship and whale, but if the blubber
#' and sublayer have been squeezed to their limits, yielding a short and
#' intense force spike as the bone compresses, then the summaries will
#' also report on the spike duration and intensity. The spike is computed
#' based on using [runmed()] on the acceleration data, with a
#' `k` value that is set to correspond to 5 ms, or to k=11, whichever
#' is larger.
#'
#' \item `"section"` to plot skin thickness, blubber thickness and sublayer thickness
#' in one panel, creating a cross-section diagram.
#'
## \item \code{"injury"} a stacked plot showing time-series traces of health
## indicators for skin extension, blubber compression, and sublayer compression.
## These take the form of dotted horizontal lines, with labels above and at the
## left side of the panel.  During times when an injury criterion is halfway met,
## e.g. that blubber stress 1/2 the value of blubber strength, then the dotted
## line is overdrawn with a thick gray line. During times when the criterion
## is exceeded, the colour shifts to black.
#'
#' \item `"threat"` a stacked plot showing time-series traces of an
#' ad-hoc measure of the possible threat to skin, blubber and sublayer.
#' The threat level is computed as the ratio
#' of stress to ultimate strength, e.g. for blubber, it is
#' `x$WCF$stress/x$parms$s[2]`. The same vertical scale is used in
#' each of the subpanels that make up the stack. Any values exceeding
#' 10 are clipped to 10, and in such a case the overall label on the vertical
#' axis will note this clipping, although it should be easy to see, because
#' the way it most often occurs is if the soft layers "bottom out" onto
#' the bone, which yields a short period of very high stress, owing to the
#' very high compression modulus of bone. Each of the curves is filled in with a light
#' gray colour for stress/strength values up to 1, and with black
#' for higher values; this makes it easy to tell at a glance whether the
#' threat level is high.
#'
#' \item `"whale acceleration"` for a time-series plot of whale acceleration.
#'
#' \item `"blubber thickness"` for a time-series plot of blubber thickness.
#'
#' \item `"sublayer thickness"` for a time-series plot of the thickness
#' of the layer interior to the blubber.
#'
#' \item `"reactive forces"` for a time-series showing the reactive
#' forces associated with skin stretching (solid) and the compression of the
#' blubber and sublayer components (dashed).
#'
#' \item `"compression stress"` for a time-series plot of the compression stress on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#'
#' \item `"skin stress"` for a time-series of skin stress in the along-skin y and z directions.
#'
#' \item `"values"` for a listing of `param` values.
#'
#' \item `"all"` for all of the above.
#'
#' \item `"default"` for a three-element plot showing `"location"`,
#' `"section"`, and `"threat"`.
#'}
#'
#' @param drawEvents Logical, indicating whether to draw lines for some events,
#' such as the moment of closest approach.
#'
#' @param colwcenter Colour used to indicate the whale centre.
#'
#' @param colwinterface Colour used to indicate the interface
#' between whale blubber and sublayer.
#'
#' @param colwskin Colour used to indicate the whale skin.
#'
#' @param cols As `colw`, but the colour to be used for the ship bow location,
#' which is drawn with a dashed line.
#'
## @param colInjury Two-element colour specification used in `"injury"`
## panels. The first colour is used to indicate values that are halfway to the
## injury crition, and the second is used to indicte values that exceed the
## criterion.
#'
#' @param colThreat Two-element colour specification used in `"threat"`.
#' The first colour is used for threat levels between 0 and 1, the second
#' for levels exceeding 1. The defaults are to use light gray and black.
#'
#' @param lwd Line width used in plots for time intervals in which damage
#' criteria are not exceeded.
#'
#' @param D Factor by which to thicken lines during times during which damage
#' criteria are exceeded.
#'
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#'
#' @param ... Ignored.
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @examples
#' ## 1. default 3-panel plot
#' t <- seq(0, 0.7, length.out=200)
#' state <- c(xs=-2, vs=knot2mps(10), xw=0, vw=0) # 10 knot ship
#' parms <- parameters() # default values
#' sol <- strike(t, state, parms)
#' par(mar=c(3,3,1,1), mgp=c(2,0.7,0), mfrow=c(1, 3))
#' plot(sol)
#' ## 2. all 12 plot types
#' par(mar=c(3,3,1,1) ,mgp=c(2,0.7,0), mfrow=c(4,3))
#' plot(sol, "all")
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom graphics abline axis box lines legend mtext par plot text
#' @importFrom grDevices hcl
#' @importFrom stats runmed
#'
#' @md
plot.strike <- function(x, which="default", drawEvents=TRUE,
                        colwcenter="black", #Slate Gray",
                        colwinterface="black", #colwinterface="Firebrick",
                        colwskin="black", #colwskin="Dodger Blue 4",
                        cols="black",
                        colThreat=c("lightgray", "black"),
                        lwd=1, D=3, debug=0, ...)
{
    g <- 9.8 # gravity
    t <- x$t
    xs <- x$xs
    vs <- x$vs
    xw <- x$xw
    vw <- x$vw
    death <- xs >= xw
    if (any(death)) {
        firstDead <- which(death)[1]
        dead <- firstDead:length(t)
        xw[dead] <- xw[firstDead]
        x$WCF$compressed[dead,1] <- 0
        x$WCF$compressed[dead,2] <- 0
    }
    showEvents <- function(xs, xw) {
        if (drawEvents) {
            ##grid()
            death <- which(xs >= xw)[1]
            tdeath <- if (is.finite(death)) t[death] else NA
            if (is.finite(tdeath)) {
                abline(v=tdeath, lwd=lwd, col="blue")
                mtext("Fatality", at=tdeath, side=3, line=0, col="blue", font=2, cex=par("cex"))
            }
            ## tclosest <- t[which.min(abs(xs-xw))]
            ## abline(v=tclosest, col="darkgreen", lwd=lwd, lty=3)
        }
    }
    all <- "all" %in% which
    if (length(which) == 1 && which == "default") {
        which <- c("location", "section", "threat")
    }

    ## Ensure that the plot type is known.
    allowed <- c("all", "location", "section", "threat", "whale acceleration",
                 "blubber thickness", "sublayer thickness",
                 "whale water force", "reactive forces", "skin stress",
                 "compression stress", "values")
    for (w in which) {
        if (!(w %in% allowed) && !length(grep("NEW", w)))
            stop("which value \"", w, "\" is not allowed; try one of: \"",
                 paste(allowed, collapse="\" \""), "\"")
    }

    ## x(t) and xw(t)
    if (all || "location" %in% which) {
        ylim <- range(c(xs, xw), na.rm=TRUE)
        plot(t, xs, type="l", xlab="Time [s]", ylab="Location [m]", col=cols, ylim=ylim, lwd=lwd, lty="84", xaxs="i")
        lines(t, xw, lwd=lwd, col=colwcenter)
        compressed <- x$WCF$compressed
        y <- xw - compressed[, 4]
        lines(t, y, col=colwinterface, lwd=lwd)
        y <- y - compressed[, 3]
        lines(t, y, col=colwinterface, lwd=lwd)
        y <- y - compressed[, 2]
        lines(t, y, col=colwskin, lwd=lwd)
        y <- y - compressed[, 1]
        lines(t, y, col=colwskin, lwd=lwd)
        ## Accelerations (quite complicated; possibly too confusing to viewer)
        k <- round(0.005 / (t[2] - t[1]))
        k <- max(k, 11L)
        if (!(k %% 2))
            k <- k + 1
        as <- derivative(vs, t)
        asmax <- max(abs(as))
        asmaxs <- max(abs(runmed(as, k))) # smoothed
        if (asmax > 2 * asmaxs) {
            peakTime <- sum(abs(as) > 0.5*(asmax+asmaxs)) * (t[2] - t[1])
            label <- sprintf("ship accel: %.1fg (%.1fms spike to %.0fg)",
                             asmaxs/g, peakTime*1e3, asmax/g)
        } else {
            label <- sprintf("ship accel: %.1fg", asmax/g)
        }
        mtext(label, side=1, line=-1.25, cex=par("cex"), adj=0.5)
        aw <- derivative(vw, t)
        awmax <- max(abs(aw))
        awmaxs <- max(abs(runmed(aw, k)))
        if (awmax > 2 * awmaxs) {
            peakTime <- sum(abs(aw) > 0.5*(awmax+awmaxs)) * (t[2] - t[1])
            label <- sprintf("whale accel: %.1fg (%.1fms spike to %.0fg)",
                             awmaxs/g, peakTime*1e3, awmax/g)
        } else {
            label <- sprintf("whale accel: %.1fg", awmax/g)
        }
        mtext(label, side=3, line=-1.25, cex=par("cex"), adj=0.5)
        showEvents(xs, xw)
    }
    if (all || "section" %in% which) {
        ##REMOVE_CRITERIA <- TRUE
        skin <- x$WCF$compressed[,1]
        blubber <- x$WCF$compressed[,2]
        sublayer <- x$WCF$compressed[,3]
        bone <- x$WCF$compressed[,4]
        maxy <- max(c(skin+blubber+sublayer+bone))
        ylim <- c(-maxy*1.2, 0)
        plot(t, -(skin+blubber+sublayer+bone), xlab="Time [s]", ylab="Whale-centred location [m]",
             type="l", lwd=lwd, ylim=ylim, xaxs="i", yaxs="i", col=colwskin)# outside skin margin
        lines(t, -(blubber+sublayer+bone), lwd=lwd, col=colwskin)
        lines(t, -(sublayer+bone), lwd=lwd, col=colwinterface)# , lty="42")
        lines(t, -bone, lwd=lwd, col=colwinterface)# , lty="42")
        ##abline(h=0, col=colwcenter, lwd=D*lwd)
        showEvents(xs, xw)
        xusr <- par("usr")[1:2]
        x0 <- xusr[1] - 0.01*(xusr[2] - xusr[1]) # snuggle up to axis
        text(x0, -0.5*x$parms$l[4], "bone", pos=4)
        text(x0, -x$parms$l[4]-0.5*x$parms$l[3], "sublayer", pos=4)
        text(x0, -x$parms$l[4]-x$parms$l[3]-0.5*x$parms$l[2], "blubber", pos=4)
        if (x$parms$l[1] > 0.1 * sum(x$parms$l))
            text(x0, -x$parms$l[4]-x$parms$l[3]-x$parms$l[2]-0.5*x$parms$l[1], "skin", pos=4)
        text(x0, 0.5*(ylim[1] - x$parms$lsum), "", pos=4)
        ##hatchPolygon <- FALSE
        ## Blubber
        ## if (FALSE) {
        ##     risk <- x$WCF$stress >= 0.5 * x$parms$UTSbeta
        ##     px <- c(t, rev(t))
        ##     py <- c(sublayer+blubber,
        ##             ifelse(rev(risk), rev(sublayer), rev(sublayer+blubber)))
        ##     if (hatchPolygon) {
        ##         polygon(px, py, border=NA, density=18, angle=45, col="lightgray")
        ##         polygon(px, py, border=NA, density=18, angle=-45, col="lightgray")
        ##     } else {
        ##         polygon(px, py, col=colInjury[1], border=NA)
        ##     }
        ##     injury <- x$WCF$stress >= x$parms$s[2]
        ##     px <- c(t, rev(t))
        ##     py <- c(sublayer+blubber,
        ##             ifelse(rev(injury), rev(sublayer), rev(sublayer+blubber)))
        ##     if (hatchPolygon) {
        ##         polygon(px, py, border=NA, density=18, angle=45, col="lightgray")
        ##         polygon(px, py, border=NA, density=18, angle=-45, col="lightgray")
        ##     } else {
        ##         polygon(px, py, col=colInjury[2], border=NA)
        ##     }
        ## }
        ## Sublayer
        ## if (FALSE) {
        ##     risk <- x$WCF$stress >= 0.5 * x$parms$s[3]
        ##     px <- t
        ##     py <- ifelse(risk, sublayer, 0)
        ##     if (hatchPolygon) {
        ##         polygon(px, py, border=NA, density=15, angle=45, col="darkgray")
        ##         polygon(px, py, border=NA, density=15, angle=-45, col="darkgray")
        ##     } else {
        ##         polygon(px, py, col=colInjury[1], border=NA)
        ##     }
        ##     injury <- x$WCF$stress >= x$parms$UTSgamma
        ##     px <- t
        ##     py <- ifelse(injury, sublayer, 0)
        ##     if (hatchPolygon) {
        ##         polygon(px, py, border=NA, density=15, angle=45, col="darkgray")
        ##         polygon(px, py, border=NA, density=15, angle=-45, col="darkgray")
        ##     } else {
        ##         polygon(px, py, col=colInjury[2], border=NA)
        ##     }
        ## }
        ##OLD Redraw lines (to clean up polygons that might have been drawn)
        ##OLD lines(t, sublayer+blubber, lwd=lwd*1.25, col="white")
        ##OLD lines(t, sublayer, lwd=lwd*1.25, col="white")
        ##OLD lines(t, sublayer+blubber, lwd=lwd, col=colwskin)
        ##OLD lines(t, sublayer, lwd=lwd, col=colwinterface)
    }
    ## if (all || "injury" %in% which) {
    ##     skinzInjury <- x$WSF$sigmaz / x$parms$s[1]
    ##     skinyInjury <- x$WSF$sigmay / x$parms$s[1]
    ##     blubberInjury <- x$WCF$stress /  x$parms$s[2]
    ##     sublayerInjury <- x$WCF$stress /  x$parms$s[3]
    ##     plot(range(t), c(0.5, 4.5), type="n", xlab="Time [s]", ylab="", axes=FALSE, xaxs="i")
    ##     mtext("Risk of Injury", side=2, line=1, cex=par("cex"))
    ##     axis(1)
    ##     box()
    ##     ## axis(2)
    ##     dy <- -1
    ##     y0 <- 4
    ##     xlim <- par('usr')[1:2]
    ##     for (i in 1:5) {
    ##         abline(h=y0+(i-1)*dy, lwd=lwd/2, lty="dotted")
    ##     }
    ##     x0 <- xlim[1] #+ 0.04 * (xlim[2] - xlim[1])
    ##     dylab <- 0.2
    ##     text(x0, y0     +dylab, "skin z", pos=4)
    ##     text(x0, y0+  dy+dylab, "skin y", pos=4)
    ##     text(x0, y0+2*dy+dylab, "blubber", pos=4)
    ##     text(x0, y0+3*dy+dylab, "sublayer", pos=4)
    ##     cex <- 1
    ##     pch <- 20
    ##     del <- 3
    ##     n <- length(t)
    ##     ## Halfway to criterion (i.e. caution)
    ##     y <- ifelse(skinzInjury >= 0.5, y0, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[1])
    ##     y <- ifelse(skinyInjury >= 0.5, y0+dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[1])
    ##     y <- ifelse(blubberInjury >= 0.5, y0+2*dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[1])
    ##     y <- ifelse(sublayerInjury >= 0.5, y0+3*dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[1])
    ##     ## Exceeding criterion (i.e. likely injury)
    ##     y <- ifelse(skinzInjury >= 1, y0, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[2])
    ##     y <- ifelse(skinyInjury >= 1, y0+dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[2])
    ##     y <- ifelse(blubberInjury >= 1, y0+2*dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[2])
    ##     y <- ifelse(sublayerInjury >= 1, y0+3*dy, NA)
    ##     lines(t, y, lwd=5*lwd, col=colInjury[2])
    ##     showEvents(xs, xw)
    ## }

    ## This version became the default 20180801T1200; older
    ## and trial ones are kept, but not documented.
    if (all || "threat" %in% which) {
        ## tcol <- hcl(h=c(30, 120, 210, 300), c=20, l=90, fixup=FALSE)
        ## tcol <- hcl(h=c(30, 120, 210, 300))
        tcol <- rep(1, 4)
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress /  x$parms$s[2]
        sublayerThreat <- x$WCF$stress /  x$parms$s[3]
        boneThreat <- x$WCF$stress /  x$parms$s[4]
        totalThreat <- skinThreat + blubberThreat + sublayerThreat + boneThreat
        worst <- max(c(skinThreat, blubberThreat, sublayerThreat, boneThreat))
        trimThreat <- 10
        trimmed <- worst > trimThreat
        if (trimmed) {
            skinThreat <- pin(skinThreat, upper=trimThreat)
            blubberThreat <- pin(blubberThreat, upper=trimThreat)
            sublayerThreat <- pin(sublayerThreat, upper=trimThreat)
            boneThreat <- pin(boneThreat, upper=trimThreat)
            worst <- pin(worst, upper=trimThreat)
        }
        dy <- round(0.5 + worst)
        ylim <- c(0, 3*dy+worst)
        plot(range(t), ylim, type="n", xlab="Time [s]", ylab="", axes=FALSE, xaxs="i")
        axis(1)
        box()
        yTicks <- pretty(c(0, worst))
        mtext(paste("Threat (stress / strength)",
                    if (trimmed) paste(" trimmed to", trimThreat) else ""),
              side=2, line=2, cex=par("cex"))
        ## Skin
        mtext("Skin", side=4, at=0, cex=0.8*par("cex"))
        y0 <- 0 # if (log) -1 else 0
        fillplot(t, y0, skinThreat, col=colThreat[2]) # high threat
        fillplot(t, y0, ifelse(skinThreat<=1, skinThreat, 1), col=colThreat[1]) # low threat
        abline(h=0, lty=3)
        axis(2, at=y0+yTicks, labels=yTicks)
        ## Blubber
        mtext("Blubber", side=4, at=dy, cex=0.8*par("cex"))
        fillplot(t, y0+dy, dy+blubberThreat,  col=colThreat[2])
        fillplot(t, y0+dy, dy+ifelse(blubberThreat<=1, blubberThreat, 1), col=colThreat[1])
        abline(h=dy, lty=3)
        axis(2, at=y0+dy+yTicks, labels=rep("", length(yTicks)), tcl=0.5)
        ## Sublayer
        mtext("Sublayer", side=4, at=2*dy, cex=0.8*par("cex"))
        fillplot(t, y0+2*dy, 2*dy+sublayerThreat,  col=colThreat[2])
        fillplot(t, y0+2*dy, 2*dy+ifelse(sublayerThreat<=1, sublayerThreat, 1), col=colThreat[1])
        abline(h=2*dy, lty=3)
        axis(2, at=y0+2*dy+yTicks, labels=yTicks)
        ## Bone
        mtext("Bone", side=4, at=3*dy, cex=0.8*par("cex"))
        fillplot(t, y0+3*dy, 3*dy+boneThreat,  col=colThreat[2])
        fillplot(t, y0+3*dy, 3*dy+ifelse(boneThreat<=1, boneThreat, 1), col=colThreat[1])
        axis(2, at=y0+3*dy+yTicks, labels=rep("", length(yTicks)), tcl=0.5)
        showEvents(xs, xw)
    }

    if (all || "threat-style-1" %in% which) {
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress /  x$parms$s[2]
        sublayerThreat <- x$WCF$stress /  x$parms$s[3]
        boneThreat <- x$WCF$stress /  x$parms$s[4]
        ylim <- c(0, max(c(skinThreat, blubberThreat, sublayerThreat, boneThreat),na.rm=TRUE))
        plot(t, skinThreat, ylim=ylim, col=1, type="l")
        lines(t, blubberThreat, col=2)
        lines(t, sublayerThreat, col=3)
        lines(t, boneThreat, col=4)
        legend("topright", lwd=par("lwd"), col=1:4, legend=c("Skin", "Blubber", "Sublayer", "Bone"))
        showEvents(xs, xw)
    }
    if (all || "threat-style-2" %in% which) {
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress /  x$parms$s[2]
        sublayerThreat <- x$WCF$stress /  x$parms$s[3]
        boneThreat <- x$WCF$stress /  x$parms$s[4]
        totalThreat <- skinThreat + blubberThreat + sublayerThreat + boneThreat
        ylim <- c(0, max(totalThreat, na.rm=TRUE))
        plot(t, skinThreat, ylim=ylim, col=1, type="l")
        lines(t, skinThreat+blubberThreat, col=2)
        lines(t, skinThreat+blubberThreat+sublayerThreat, col=3)
        lines(t, skinThreat+blubberThreat+sublayerThreat+boneThreat, col=4)
        legend("topright", title="Cumulative threat",
               lwd=4*par("lwd"), col=rev(1:4),
               legend=rev(c("Skin", "Blubber", "Sublayer", "Bone")))
        showEvents(xs, xw)
    }
    if (all || "threat-style-3" %in% which) {
        ##tcol <- hcl(h=c(30, 120, 210, 300), c=20, l=90, fixup=FALSE)
        tcol <- hcl(h=c(30, 120, 210, 300))
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress /  x$parms$s[2]
        sublayerThreat <- x$WCF$stress /  x$parms$s[3]
        boneThreat <- x$WCF$stress /  x$parms$s[4]
        totalThreat <- skinThreat + blubberThreat + sublayerThreat + boneThreat
        ylim <- c(0, max(totalThreat, na.rm=TRUE))
        plot(range(t), ylim, type="n", xlab="Time [s]", ylab="Cumulative Threat")
        fillplot(t, 0, skinThreat, col=tcol[1])
        fillplot(t, skinThreat, skinThreat+blubberThreat, col=tcol[2])
        fillplot(t, skinThreat+blubberThreat, skinThreat+blubberThreat+sublayerThreat, col=tcol[3])
        fillplot(t, skinThreat+blubberThreat+sublayerThreat, skinThreat+blubberThreat+sublayerThreat+boneThreat, col=tcol[4])
        legend("topright", lwd=8*par("lwd"), col=rev(tcol),
               legend=rev(c("Skin", "Blubber", "Sublayer", "Bone")))
        showEvents(xs, xw)
    }
    if (all || "threat-style-4" %in% which) {
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress /  x$parms$s[2]
        sublayerThreat <- x$WCF$stress /  x$parms$s[3]
        boneThreat <- x$WCF$stress /  x$parms$s[4]
        dy <- -2.25
        y0 <- 8.0
        plot(range(t), c(1, 10),
             type="n", xlab="Time [s]", ylab="", axes=FALSE, xaxs="i")
        mtext("Threat of Injury", side=2, line=2, cex=par("cex"))
        axis(1)
        box()
        ## axis(2)
        xlim <- par('usr')[1:2]
        for (i in 1:4) {
            abline(h=y0+(i-1)*dy, lwd=lwd/2, lty="dotted")
        }
        x0 <- xlim[1] #+ 0.04 * (xlim[2] - xlim[1])
        ##dylab <- 0.2
        mtext("bone", side=2, at=y0-0.5*dy, line=0.25, cex=par("cex"))
        mtext("sublayer", side=2, at=y0+dy-0.5*dy, line=0.25, cex=par("cex"))
        mtext("blubber", side=2, at=y0+2*dy-0.5*dy, line=0.25, cex=par("cex"))
        mtext("skin", side=2, at=y0+3*dy-0.5*dy, line=0.25, cex=par("cex"))
        n <- length(t)

        for (l in 0:3) {
            abline(h=y0+l*dy + 0.5, lty='dotted', lwd=0.5*lwd)
            abline(h=y0+l*dy + 1.0, lty='dotted', lwd=0.5*lwd)
            abline(h=y0+l*dy + 1.5, lty='dotted', lwd=0.5*lwd)
        }

        ## Make polygon fill down to lowest value
        polygonWindow <- function(x, y, floor=min(y), col="pink")
        {
            n <- length(x)
            xx <- c(x[1], x, x[n])
            yy <- c(floor, y, floor)
            polygon(xx, yy, col=col)
        }
        thicker <- 2.0
        thinner <- 0.8
        ## Bone at top
        polygonWindow(t, y0+0*dy+boneThreat, col=colThreat[3], floor=y0)
        polygonWindow(t, y0+0*dy+ifelse(boneThreat<1, boneThreat, 1), col=colThreat[2], floor=y0)
        polygonWindow(t, y0+0*dy+ifelse(boneThreat<0.5, boneThreat, 0.5), col=colThreat[1], floor=y0)
        lines(t, boneThreat + y0 + 0 * dy, lwd=thicker*lwd, col="white")
        lines(t, boneThreat + y0 + 0 * dy, lwd=thinner*lwd)
        ## Sublayer below bone
        polygonWindow(t, y0+1*dy+sublayerThreat, col=colThreat[3], floor=y0+dy)
        polygonWindow(t, y0+1*dy+ifelse(sublayerThreat<1, sublayerThreat, 1), col=colThreat[2], floor=y0+dy)
        polygonWindow(t, y0+1*dy+ifelse(sublayerThreat<0.5, sublayerThreat, 0.5), col=colThreat[1], floor=y0+dy)
        lines(t, sublayerThreat + y0 + 1 * dy, lwd=thicker*lwd, col="white")
        lines(t, sublayerThreat + y0 + 1 * dy, lwd=thinner*lwd)
        ## Blubber second from bottom
        polygonWindow(t, y0+2*dy+blubberThreat, col=colThreat[3], floor=y0+2*dy)
        polygonWindow(t, y0+2*dy+ifelse(blubberThreat<1, blubberThreat, 1), col=colThreat[2], floor=y0+2*dy)
        polygonWindow(t, y0+2*dy+ifelse(blubberThreat<0.5, blubberThreat, 0.5), col=colThreat[1], floor=y0+2*dy)
        lines(t, blubberThreat + y0 + 2 * dy, lwd=thicker*lwd, col="white")
        lines(t, blubberThreat + y0 + 2 * dy, lwd=thinner*lwd)
        ## Skin at bottom
        polygonWindow(t, y0+3*dy+skinThreat, col=colThreat[3], floor=y0+3*dy)
        polygonWindow(t, y0+3*dy+ifelse(skinThreat<1, skinThreat, 1), col=colThreat[2], floor=y0+3*dy)
        polygonWindow(t, y0+3*dy+ifelse(skinThreat<0.5, skinThreat, 0.5), col=colThreat[1], floor=y0+3*dy)
        lines(t, skinThreat + y0 + 3 * dy, lwd=thicker*lwd, col="white")
        lines(t, skinThreat + y0 + 3 * dy, lwd=thinner*lwd)

        ## Axes
        axis(4,y0+0*dy+c(0,0.5,1,1.5),labels=c("L", "M", "H", "E"))
        axis(4,y0+1*dy+c(0,0.5,1,1.5),labels=c("L", "M", "H", "E"))
        axis(4,y0+2*dy+c(0,0.5,1,1.5),labels=c("L", "M", "H", "E"))
        axis(4,y0+3*dy+c(0,0.5,1,1.5),labels=c("L", "M", "H", "E"))

        showEvents(xs, xw)
    }
    if (all || "whale acceleration" %in% which) {
        a <- derivative(vw, t)
        plot(t, a, xlab="Time [s]", ylab="Whale accel. [m/s^2]", type="l", lwd=lwd, xaxs="i")
        showEvents(xs, xw)
    }

    if (all || "blubber thickness" %in% which) {
        WCF <- x$WCF
        y <- WCF$compressed[, 2]
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Blubber thickness [m]", type="l", lwd=lwd, ylim=ylim, xaxs="i")
        showEvents(xs, xw)
    }
    if (all || "sublayer thickness" %in% which) {
        WCF <- x$WCF
        y <- WCF$compressed[, 3]
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Sublayer thickness [m]", type="l", lwd=lwd, ylim=ylim, xaxs="i")
        showEvents(xs, xw)
    }
    if (all || "whale water force" %in% which) {
        plot(t, whaleWaterForce(vw, x$parms) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=lwd, xaxs="i")
        showEvents(xs, xw)
    }
    if (all || "reactive forces" %in% which) {
        SF <- x$WSF$force
        CF <- x$WCF$force
        WWF <- x$WWF
        ylim <- range(c(WWF, SF, CF), na.rm=TRUE)/1e6
        plot(t, SF/1e6, type="l", xlab="Time [s]", ylab="Forces [MN]", lwd=lwd, ylim=ylim, xaxs="i")
        ##lines(t, CF/1e6, lty="dotdash", lwd=lwd)
        lines(t, CF/1e6, lty="dotted", lwd=lwd)
        ## legend("topleft", col=1:2, lwd=lwd, legend=c("Skin", "Blubber"))
        mtext(expression(" "*F[E]), side=3, line=-1.2, adj=0, cex=par("cex"))
        mtext(" (solid)", side=3, line=-2.2, adj=0, cex=par("cex"))
        mtext(expression(F[C]*" "), side=3, line=-1.2, adj=1, cex=par("cex"))
        mtext(" (dotted) ", side=3, line=-2.2, adj=1, cex=par("cex"))
        ## lines(t, WWF / 1e6 , lwd=lwd)
        showEvents(xs, xw)
    }
    ## if (all || "skin force" %in% which) {
    ##     Fs <- whaleSkinForce(xs, xw, x$parms)
    ##     plot(t, Fs$force/1e6, type="l", xlab="Time [s]", ylab="Skin Force [MN]", lwd=lwd)
    ##     showEvents(xs, xw)
    ## }
    if (all || "skin stress" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        ylim <- range(c(Fs$sigmay, Fs$sigmaz)/1e6)
        plot(t, Fs$sigmay/1e6, type="l", xlab="Time [s]", ylab="Skin Stress [MPa]", lwd=lwd, ylim=ylim, xaxs="i")
        lines(t, Fs$sigmaz/1e6, lty="dotted", lwd=lwd)
        mtext(" horiz.", side=3, line=-1.2, adj=0, cex=par("cex"))
        mtext(" (solid)", side=3, line=-2.2, adj=0, cex=par("cex"))
        mtext("vert. ", side=3, line=-1.2, adj=1, cex=par("cex"))
        mtext("(dotted) ", side=3, line=-2.2, adj=1, cex=par("cex"))
        showEvents(xs, xw)
    }
    ##OLDif (all || "compression force" %in% which) {
    ##OLD    plot(t, x$WCF$force/1e6, type="l", xlab="Time [s]", ylab="Compress. Force [MN]", lwd=lwd, xaxs="i")
    ##OLD    showEvents(xs, xw)
    ##OLD}
    if (all || "compression stress" %in% which) {
        force <- x$WCF$force
        stress <- force / (x$parms$Lz*x$parms$Ly)
        plot(t, stress/1e6, type="l", xlab="Time [s]", ylab="Compress. Stress [MPa]", lwd=lwd, xaxs="i")
        showEvents(xs, xw)
    }
    if (all || "values" %in% which) {
        omar <- par("mar")
        par(mar=rep(0, 4))
        ## Only take the vectors, thus ignoring stressFromStrain, a function
        parms <- x$parms[unlist(lapply(x$parms, function(p) is.vector(p)))]
        parms <- lapply(parms, function(x) signif(x, 4))
        parms["engineForce"] <- NULL # inserted during calculation, not user-supplied
        parms["lsum"] <- NULL # inserted during calculation, not user-supplied

        parms <- lapply(parms, function(p) deparse(p))
        ## parms$ms <- x$parms$ms
        parms$vs_knots <- knot2mps(x$vs[1])
        names <- names(parms)
        values <- unname(unlist(parms))
        ##values[["mw"]] <- x$parms$mw
        ##values[["vs"]] <- x$vs

        ## values <- paste(deparse(parms), collapse=" ")
        ## values <- gsub("^list\\(", "", values)
        ## values <- gsub(")$", "", values)
        ## values <- strsplit(gsub(".*\\((.*)\\)$", "\\1", values), ",")[[1]]
        ## values <- gsub("^ *", "", values)
        ## values <- gsub("([0-9])L$", "\\1", values) # it's ugly to see e.g. 1 as 1L

        n <- 1 + length(values)
        plot(1:n, 1:n, type="n", xlab="", ylab="", axes=FALSE)
        o <- order(names(parms), decreasing=TRUE)
        for (i in seq_along(values))
            text(1, i+0.5, paste(names[o[i]], "=", values[o[i]]), pos=4, cex=1)
        par(mar=omar)
    }
}

#' Summarize the parameters of a simulation, and its results
#'
#' @param object An object created by [strike()].
#'
#' @param style Character value indicating the representation to be used.
#' If `style` is `"text"`, then the output is in a textual format. If it
#' is `"latex"`, then the output is in latex format.
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom xtable xtable
#'
#' @md
summarize <- function(object, style="text")
{
    parm <- object$parm
    parm$a1 <- parm$a[1]
    parm$a2 <- parm$a[2]
    parm$a3 <- parm$a[3]
    parm$a4 <- parm$a[4]
    parm$a <- NULL
    parm$b1 <- parm$b[1]
    parm$b2 <- parm$b[2]
    parm$b3 <- parm$b[3]
    parm$b4 <- parm$b[4]
    parm$b <- NULL
    parm$l1 <- parm$l[1]
    parm$l2 <- parm$l[2]
    parm$l3 <- parm$l[3]
    parm$l4 <- parm$l[4]
    parm$l <- NULL
    parm$s1 <- parm$s[1]
    parm$s2 <- parm$s[2]
    parm$s3 <- parm$s[3]
    parm$s4 <- parm$s[4]
    parm$s <- NULL

    names <- names(parm)
    ## Remove some calculated things.
    parm <- parm[!(names %in% c("engineForce", "stressFromStrain", "lsum"))]
    names <- names(parm)
    parm <- unname(parm)
    meaning <- rep("", length(parm))
    o <- order(names)
    parm <- parm[o]
    names <- as.character(names[o])
    ## FIXME: update to l, a, b, and s.
    meaning[names=="alpha"] <- "Whale skin thickness [m]."
    meaning[names=="beta"] <- "Blubber thickness [m]."
    meaning[names=="gamma"] <- "Sub-layer thickness [m]."
    meaning[names=="Ealpha"] <- "Whale skin extension modulus [MPa]."
    meaning[names=="Ebeta"] <- "Blubber compression modulus [MPa]."
    meaning[names=="Egamma"] <- "Sub-layer compression modulus [MPa]."
    meaning[names=="UTSalpha"] <- "Whale skin ultimate tension strength [MPa]."
    meaning[names=="UTSbeta"] <- "Blubber ultimate compression strength [MPa]."
    meaning[names=="UTSgamma"] <- "Sub-layer ultimate compression strength [MPa]."
    meaning[names=="Ly"] <- "Width of impact zone [m]."
    meaning[names=="Lz"] <- "Height of impact zone [m]."
    meaning[names=="lw"] <- "Whale length [m]."
    meaning[names=="mw"] <- "Whale mass [m]."
    meaning[names=="Sw"] <- "Whale surface area [sq. m]."
    meaning[names=="mw"] <- "Whale mass [m]."
    meaning[names=="ms"] <- "Ship mass [m]."
    meaning[names=="Ss"] <- "Ship wetted surface area [sq. m]."
    meaning[names=="theta"] <- "Angle of skin-compression bevel [deg]."
    meaning[names=="Cw"] <- "Whale drag coefficient [unitless]."
    meaning[names=="Cs"] <- "Ship drag coefficient [unitless]."
    meaning[names=="a1"] <- "Skin stress factor [Pa]."
    meaning[names=="a2"] <- "Blubber stress factor [Pa]."
    meaning[names=="a3"] <- "Sublayer stress factor [Pa]."
    meaning[names=="a4"] <- "Bone stress factor [Pa]."
    meaning[names=="b1"] <- "Skin stress nonlinearity [unitless]."
    meaning[names=="b2"] <- "Blubber stress nonlinearity [unitless]."
    meaning[names=="b3"] <- "Sublayer stress nonlinearity [unitless]."
    meaning[names=="b4"] <- "Bone stress nonlinearity [unitless]."
    meaning[names=="l1"] <- "Skin thickness [m]."
    meaning[names=="l2"] <- "Blubber thickness [m]."
    meaning[names=="l3"] <- "Sublayer thickness [m]."
    meaning[names=="l4"] <- "Bone thickness [m]."
    meaning[names=="s1"] <- "Skin strength [Pa]."
    meaning[names=="s2"] <- "Blubber strength [Pa]."
    meaning[names=="s3"] <- "Sublayer strength [Pa]."
    meaning[names=="s4"] <- "Bone strength [Pa]."
    if (style == "latex") {
        ## names[which(names == "Ealpha")] <- "E_alpha"
        ## names[which(names == "Ebeta")] <- "E_beta"
        ## names[which(names == "Egamma")] <- "E_gamma"
        ## names[which(names == "UTSalpha")] <- "UTS_alpha"
        ## names[which(names == "UTSbeta")] <- "UTS_beta"
        ## names[which(names == "UTSgamma")] <- "UTS_gamma"
        names[which(names == "Cs")] <- "$C_s$"
        names[which(names == "Cw")] <- "$C_w$"
        names[which(names == "lw")] <- "$l_w$"
        names[which(names == "Ly")] <- "$L_y$"
        names[which(names == "Lz")] <- "$L_z$"
        names[which(names == "Ss")] <- "$S_s$"
        names[which(names == "Sw")] <- "$S_w$"
        names[which(names == "ms")] <- "$M_s$"
        names[which(names == "mw")] <- "$M_w$"
        names[which(names == "a1")] <- "$a_1$"
        names[which(names == "a2")] <- "$a_2$"
        names[which(names == "a3")] <- "$a_3$"
        names[which(names == "a4")] <- "$a_4$"
        names[which(names == "b1")] <- "$b_1$"
        names[which(names == "b2")] <- "$b_2$"
        names[which(names == "b3")] <- "$b_3$"
        names[which(names == "b4")] <- "$b_4$"
        names[which(names == "l1")] <- "$l_1$"
        names[which(names == "l2")] <- "$l_2$"
        names[which(names == "l3")] <- "$l_3$"
        names[which(names == "l4")] <- "$l_4$"
        names[which(names == "s1")] <- "$s_1$"
        names[which(names == "s2")] <- "$s_2$"
        names[which(names == "s3")] <- "$s_3$"
        names[which(names == "s4")] <- "$s_4$"
        parm <- data.frame(Symbol=names,
                           Meaning=meaning,
                           Value=as.matrix(unname(unclass(parm))),
                           Reference="")
        colnames(parm) <- c("Item", "Meaning", "Value", "Source")
        print(xtable(parm, label="table:model_parameters",
                     caption="Model parameters.",
                     display=c("d", "s", "s", "G", "s")),
              math.style.exponents=TRUE,
              sanitize.text.function=function(x) x,
              include.rownames=FALSE)
    } else if (style == "text") {
        parm <- data.frame(Name=names, Meaning=meaning, Value=as.matrix(unname(unclass(parm))))
        colnames(parm) <- c("Item", "Meaning", "Value")
        print(parm, include.rownames=FALSE)
    }
}

