library(deSolve)
library(xtable)

#' Whale blubber stress-strain relationship
#'
#' This is a data frame with elements \code{strain} and \code{stess},
#' found by digitizing (accurate to perhaps 1 percent) the curve shown in Figure 2.13
#' of Raymond (2007). It is used to develop a stress-strain relationship used
#' by \code{\link{parameters}}, as shown in \dQuote{Examples}.
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
#' @docType data
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
#' The documentation for \code{\link{strike}} provides
#' an example of using the main functions of this package,
#' and so it is a good place to start. A companion manuscript
#' is intended to provide more detail about the analysis
#' and the context.
#'
#' @section Further reading:
#' \itemize{
#'
#' \item
#' Daoust, Pierre-Yves, Émilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' “Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017.” Canadian Wildlife Health Cooperative, Marine Animal Response
#' Socieity, and Fisheries and Oceans Canada, 2018.
#' http://publications.gc.ca/site/eng/9.850838/publication.html.
#'
#' \item
#' Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. “Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis).” Journal of Mammalogy 93,
#' no. 5 (2012): 1342–54. https://doi.org/10.1644/11-MAMM-A-297.1.
#'
#' \item
#' Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt, Adam
#' P. Summers, and Petra Ditsche. “Mechanical Properties of Harbor Seal Skin and
#' Blubber − a Test of Anisotropy.” Zoology 126 (2018): 137–44.
#' https://doi.org/10.1016/j.zool.2017.11.002.
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
#' MAN Diesel & Turbo. “Basic Principles of Propulsion.” MAN Diesel & Turbo, 2011. https://marine.mandieselturbo.com/docs/librariesprovider6/propeller-aftship/basic-principles-of-propulsion.pdf?sfvrsn=0.
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
#' Hampshire, 2007. https://scholars.unh.edu/thesis/309.
#'
#'}
#'
#' @docType package
#' @name whalestrike
NULL

#' Create a function for stress in laminated layers
#'
#' Assuming that unforced layer thicknesses are \eqn{l_i}, and that within each layer
#' the stress is given by
#' \deqn{a_i*[exp(b_i*\epsilon)-1]}
#' with strain \eqn{\epsilon}
#' being \eqn{\Delta l_i/l_i}, then the change in the
#' total thickness \eqn{L=\sum l_i} obeys
#' \deqn{0 = \Delta L - \sum (l_i /b_i) ln(1+\sigma / a_i)}
#' where \eqn{ln} is the natural logarithm. This formula rests on the assumption
#' that the stress, \eqn{\sigma}, is the same in each layer.
#' This expression is not easily inverted to get
#' \eqn{\sigma} in terms of \eqn{\Delta L}, but it may be solved
#' easily for particular numerical vaues, using \code{\link{uniroot}}.
#' This is done for a sequence of \code{N} values of strain \eqn{\epsilon}
#' that range from 0 to 1. Then a cubic spline is created to represent these
#' the relationship between \eqn{\sigma} and \eqn{\Delta L} values,
#' and this spline function is returned.
#' The purpose is to speed up processing of simulations carried out
#' with \code{\link{strike}}. In the authors tests, the speedup was by
#' an order of magnitude, and the loss of accuracy was
#' negligible, amounting to fractional errors in the 8th decimal place.
#'
#' @param l vector of layer thicknesses
#' @param a vector of multipliers
#' @param b vector of e-fold parameters
#' @param N integer specifying how many segments to use in the spline
#'
#' @return A spline function, created with \code{\link{splinefun}},
#' that returns stress as a function of total strain of the
#' system of compressing layers. For the purposes of the whale-strike
#' analysis, the strain should be between 0 and 1, i.e. there is
#' no notion of compressing blubber, etc. to negative thickess.
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
stressFromStrainFunction <- function(l, a, b, N=1e3)
{
    fcn <- function(sigma)
        DL - sum((l/b) * log(1 + sigma / a))
    L <- sum(l)
    sigma <- rep(NA, N)
    epsilon <- seq(0, 1, length.out=N)
    for (i in seq_along(epsilon)) {
        DL <- epsilon[i] * L
        sigma[i] <- uniroot(fcn, interval=c(0, 1e8))$root
        ##cat("i=", i, ", epsilon=", epsilon[i], ", sigma=", sigma[i], "\n")
    }
    splinefun(epsilon, sigma)
}


#' Control parameters for whale strike simulation
#'
#' Assembles control parameters into a list suitable for passing to \code{\link{strike}}
#' and the functions that it calls. If \code{file} is provided, then all the other
#' arguments are read from that source.
#' Below are some sources cited in the discussion of the function arguments.
#' @template ref_daoust
#' @template ref_grear
#' @template ref_raymond
#'
#' @param ms Ship mass [kg].
#' @param Ss Ship wetted area [m^2]. This, together with \code{Cs}, is used by
#' used by \code{\link{shipWaterForce}} to estimate ship drag force. If \code{Ss}
#' is not given, then an esimate is made by calling \code{\link{shipAreaFromMass}} with
#' the provided value of \code{ms}.
#' @param Ly Ship impact horizontal extent [m]; defaults to 0.5m if not specified.
#' @param Lz Ship impact vertical extent [m]; defaults to 1.5m if not specified.
#' @param lw Whale length [m]. If not supplied, \code{\link{whaleLengthFromMass}} is used
#' to calculate this, given \code{lm}, but if neither \code{mw} nor \code{ml} is provided,
#' an error is reported. The length is used by \code{\link{whaleAreaFromLength}} to
#' calculate area, which is needed for the water drag calculation done by
#' \code{\link{whaleWaterForce}}.
#' @param mw Whale mass [kg]. If not provided, this is calculated from whale
#' length, using \code{\link{whaleMassFromLength}} with \code{type="wetted"}.
#' @param Sw Whale surface area [m^2]. If not provided, this is calculated
#' from whale length using \code{\link{whaleAreaFromLength}}.
#' @param l Numerical vector of length 4, giving thickness [m] of skin, blubber,
#' sublayer, and bone. If not provided, this is set to
#' \code{c(0.025, 0.16, 0.5, 0.2)}.
#' The skin thicknes default of 0.025 m represents the 0.9-1.0 inch value
#' stated in Section 2.2.3 of Raymond (2007).
#' The blubber default of 0.16 m is a rounded average of the values inferred
#' by whale necropsy, reported in Appendix 2 of Daoust et al., 2018.
## mean(c(17,14,18.13,18,21.25,16.75,13.33,7)) # cm
## [1] 15.6825
#' The sublayer default of 0.5 m may be reasonable at some spots on the whale body.
#' The bone default of 0.1 m may be reasonable at some spots on the whale body.
#' @param a,b Numerical vectors of length 4, giving values to use in the
#' stress-strain law \code{stress=a*(exp(b*strain)-1)}. \code{a} is in Pa
#' and \code{b} is unitless.  Note that \code{a*b} is the local modulus at
#' low strain, and that \code{b} is the efolding scale for nonlinear increase
#' in stress with strain. This exponential relationship has been mapped out
#' for whale blubber, using a curve fit to Figure 2.13 of Raymond (2007), and
#' these values are used for the second layer (blubber); see
#' \link{raymond2007} for how the fit was done.
#' If not provided, \code{a} defaults to
#' \code{c(17.80e6/0.1, 1.64e5, 1.64e5, 854.2e6/0.1)}
#' and \code{b} defaults to
#' \code{c(0.1, 2.47, 2.47, 0.1)}.
#' The skin defaults are set up to give a linear shape (since \code{b} is small)
#' with the \code{a*b} product
#' being 17.8e6 Pa, which is the adult-seal value
#' given in Table 3 of Grear et al. (2017).
#' The blubber defaults are from a regression of the stress-strain
#' relationship shown in Figure 2.13 of Raymond (2007).
#' The sublayer defaults are set to match those of blubber, lacking
#' any other information.
#' The bone default for \code{b} is small, to set up a linear function,
#' and \code{a*b} is set to equal 8.54e8 Pa,
#' given in Table 2.3 of Raymond (2007).
#' @param s Numerical vector of length 4, giving the ultimate strengths [Pa] of
#' skin, blubber, sublayer, and bone, respectively. If not provided, the
#' value is set to
#' \code{c(19.56e6, 4.37e5, 4.37e5, 22.9e6)},
#' with reasoning as follows.
#' The skin default of 1.96e7 Pa
#' is a rounded value from Table 3 of Grear et al. (2018) for adult seal skin strength at
#' an orientation of 0 degrees.  The blubber value of
#' 4.37e5 Pa is inferred by
#' multiplying Raymond's (2007) Figure 2.13 elastic modulus of 6.36e5 Pa
#' by the ratio 0.97/1.41 determined for adult seal strength/modulus, as reported
#' in Table 3 of Grear et al. (2018).
#' The sublayer value is taken to be identical to the blubber value, lacking
#' more specific information.
#' The bone default o 2.29e7 Pa is from Table 2.3 of Raymond (2007).
#' @param theta Whale skin deformation angle [deg]; defaults to 45deg if not supplied.
#' @param Cs Drag coefficient for ship [dimensionless],
#' used by \code{\link{shipWaterForce}} to estimate ship drag force. Defaults
#' to 1e-2, which is 4 times the frictional coefficient of 2.5e-3
#' inferred from Figure 4 of Manen and van Oossanen (1988), assuming
#' a Reynolds number of 5e7, computed from speed 5m/s, lengthscale 10m
#' and viscosity 1e-6 m^2/s. (The factor of 4 is under the assumption
#' that frictional drag is about a quarter of total drag.)
#' The drag force is computed with \code{\link{shipWaterForce}}.
#' @param Cw Drag coefficient for whale [dimensionless],
#' used by \code{\link{whaleWaterForce}} to estimate whale drag force.
#' Defaults to 2.5e-3, for Reynolds number 2e7, computed from speed
#' 2 m/s, lengthscale 5m (between radius and length) and
#' viscosity 1e-6 m^2/s.  The drag force is computed with
#' \code{\link{whaleWaterForce}}.
#' @param file Optional name a comma-separated file that holds all of the
#' previous values, except \code{Cs} and \code{Cw}. If provided,
#' then other parameters (except \code{Cs} and \code{Cw}) are
#' ignored, because values are sought from the file. The purpose of
#' this is in shiny apps that want to save a simulation framework.
#' The file should be saved \code{\link{write.csv}} with
#' \code{row.names=FALSE}.
#'
#' @return
#' A named list holding the parameters, with defaults and alternatives reconciled
#' according to the system described above.
#'
#' @examples
#' parms <- parameters(ms=20e3, lw=13)
#' epsilon <- seq(0, 0.3, length.out=100)
#' strain <- parms$stressFromStrain(epsilon)
#' plot(epsilon, strain, type="l")
parameters <- function(ms=20e3, Ss, Ly=0.5, Lz=1.0,
                       lw=13, mw, Sw,
                       l, a, b, s,
                       theta=45,
                       Cs=0.01, Cw=0.0025, file)
{
    if (!missing(file)) {
        rval <- read.csv(file)
        rval$Ss <- shipAreaFromMass(rval$ms)
        rval$mw <- whaleMassFromLength(rval$lw)
        rval$Sw <- whaleAreaFromLength(rval$lw, type="wetted")
        if ("gammaType" %in% names(rval)) {
            if (rval$gammaType == 1) {
                rval$gamma <- rval$gamma1
                rval$Egamma <- rval$Egamma1
                rval$UTSgamma <- rval$UTSgamma1
            } else {
                rval$gamma <- rval$gamma2
                rval$Egamma <- rval$Egamma2
                rval$UTSgamma <- rval$UTSgamma2
            }
            rval$gamma1 <- rval$Egamma1 <- rval$UTSgamma1 <- NULL
            rval$gamma2 <- rval$Egamma2 <- rval$UTSgamma2 <- NULL
            rval$gammaType <- NULL
        }
        rval$tmax <- NULL
        rval$vs <- NULL
        rval$Cs <- Cs
        rval$Cw <- Cw
        o <- sort(names(rval))
        rval <- rval[o]
    } else {
        if (ms <= 0)
            stop("ship mass (ms) must be a positive number")
        if (missing(Ss))
            Ss <- shipAreaFromMass(ms)
        if (Ss <= 0)
            stop("ship wetted area (Ss) must be a positive number, not ", Ss)
        if (Ly <= 0)
            stop("impact width (Ly) must be a positive number, not ", Ly)
        if (Lz <= 0)
            stop("impact height (Lz) must be a positive number, not ", Lz)
        if (lw <= 0)
            stop("Whale length (lw) must be a positive number")
        if (missing(mw))
            mw <- whaleMassFromLength(lw)
        if (missing(Sw))
            Sw <- whaleAreaFromLength(lw, type="wetted")
        if (missing(l))
            l <- c(0.025, 0.16, 0.5, 0.1)
        if (missing(a))
            a <- c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)
        if (missing(b))
            b <- c(0.1, 2.54, 2.54, 0.1)
        if (missing(s))
            s <- c(1.96e7, 4.37e5, 4.37e5, 2.29e7)
        ## Value checks
        if (any(l <= 0) || length(l) != 4)
            stop("'l' must be 4 positive numbers")
        if (any(a <= 0) || length(a) != 4)
            stop("'a' must be 4 positive numbers")
        if (any(b <= 0) || length(b) != 4)
            stop("'b' must be 4 positive numbers")
        if (theta < 0 || theta > 89)
            stop("whale skin deformation angle (theta) must be between 0 and 89 deg, not ", theta)
        if (Cs <= 0)
            stop("ship resistance parameter (Cs) must be positive, not ", Cs)
        if (Cw <= 0)
            stop("ship resistance parameter (Cw) must be positive, not ", Cw)

        ##DELETE if (Ealpha< 0)
        ##DELETE     stop("whale skin elastic modulus (Ealpha) must be positive, not ", Ealpha)
        ##DELETE if (UTSalpha< 0)
        ##DELETE     stop("whale skin elastic modulus (UTSalpha) must be positive, not ", UTSalpha)
        ##DELETE if (beta < 0)
        ##DELETE     stop("whale blubber thickness (beta) must be positive, not ", beta)
        ##DELETE if (Ebeta < 0)
        ##DELETE     stop("whale blubber elastic modulus (Ebeta) must be positive, not ", Ebeta)
        ##DELETE if (UTSbeta < 0)
        ##DELETE     stop("whale blubber ultimate strength (UTSbeta) must be positive, not ", UTSbeta)
        ##DELETE if (gamma < 0)
        ##DELETE     stop("whale sublayer thickness (gamma) must be positive, not ", gamma)
        ##DELETE if (Egamma < 0)
        ##DELETE     stop("whale sublayer elastic modulus (Egamma) must be positive, not ", Egamma)
        ##DELETE if (UTSgamma < 0)
        ##DELETE     stop("whale sublayer ultimate strength (UTSgamma) must be positive, not ", UTSgamma)

        ## overall-stress from overall-strain function
        stressFromStrain <- stressFromStrainFunction(l, a, b)
        rval <- list(ms=ms, Ss=Ss,
                     Ly=Ly, Lz=Lz,
                     mw=mw, Sw=Sw, lw=lw,
                     l=l, lsum=sum(l), a=a, b=b, s=s,
                     ##alpha=alpha, Ealpha=Ealpha, UTSalpha=UTSalpha,
                     theta=theta,
                     ##Ebeta=Ebeta, beta=beta, UTSbeta=UTSbeta,
                     ##Egamma=Egamma, gamma=gamma, UTSgamma=UTSgamma,
                     Cs=Cs, Cw=Cw,
                     stressFromStrain=stressFromStrain)
    }
    class(rval) <- "parameters"
    rval
}


#' Whale mass inferred from length
#'
#' Calculate an estimate of the mass of a whale, based on its length.
#'
#' The permitted values for \code{model} are as follows.
#'\itemize{
#' \item \code{"moore2005"} yields
#' \eqn{242.988 * exp(0.4 * length)}{242.988 * exp(0.4 * L)},
#' which (apart from a unit change on \code{L}) is the regression equation
#' shown above Figure 1d in Moore et al. (2005) for right whales. A
#' difficult in the Moore et al. (2005) use of a single nonzero digit
#' in the multiplier on \code{L} is illustrated in \dQuote{Examples}.
#'
#' \item \code{"fortune2012atlantic"} yields the formula
#' \eqn{exp(-10.095 + 2.825*log(100*L))}{exp(-10.095 + 2.825*log(100*L))}
#' for North Atlantic right whales, according to corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'
#' \item \code{"fortune2012pacific"} yields the formula
#' \eqn{exp(-12.286 + 3.158*log(100*L))}{exp(-12.286 + 3.158*log(100*L))}
#' for North Pacific right whales, according to corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'}
#'
#' @param L Whale length in m.
#' @param model Character string specifying the model (see \dQuote{Details}).
#' @return Mass in kg.
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
#' See \link{whalestrike} for a list of references.
whaleMassFromLength <- function(L, model="fortune2012atlantic")
{
    if (model == "moore2005")
        242.988 * exp(0.4 * L)
    else if (model == "fortune2012atlantic")
        exp(-10.095 + 2.825*log(100*L))
    else if (model == "fortune2012pacific")
        exp(-12.286 + 3.158*log(100*L))
    else
        stop("unrecognized model '", model, "'")
}

#' Compute whale length from mass
#'
#' This works by inverting \code{\link{whaleMassFromLength}} using
#' \code{\link{uniroot}}.
#'
#' @param M Whale mass [kg].
#' @param model Character string specifying the model, with permitted
#' values \code{"moore2005"}, \code{"fortune2012atlantic"} and
#" \code{"fortune2012pacific"}. See the documentation
#' for \code{\link{whaleMassFromLength}} for the details of these
#' formulations.
#'
#' @return Whale length [m].
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleLengthFromMass <- function(M, model="fortune2012atlantic")
{
    rval <- rep(NA, length(M))
    for (i in seq_along(M))
        rval[i] <- uniroot(function(x) M[i] - whaleMassFromLength(x, model), c(0.1, 100))$root
    rval
}

#' Whale projected area, as function of length
#' @param L length in m
#' @param type character string indicating the type of area, with
#' \code{"projected"} for side-projected area using 0.143L^2,
#' and \code{"wetted"}
#' for submerged surface wetted, calculated by spinning
#' the necropsiy side-view presented in Daoust et al. (2018)
#' along the animal axis, yielding
#' 0.451 * (0.8715 * L)^2, where the direct factor on L
#' was developed by comparing similarly-predicted masses to
#' the results of \code{\link{whaleMassFromLength}}, as described
#' in [2].
#'
#' @references
#' 1. Dan Kelley's internal document \code{dek/20180623_whale_area.Rmd}, available
#' upon request.
#'
#' 2. Dan Kelley's internal document \code{dek/20180707_whale_mass.Rmd}, available
#' upon request.
whaleAreaFromLength <- function(L, type="wetted")
{
    ## below from dek/20180623_whale_area.Rmd
    ## Projected area, with fins: 0.1466L2 where L is body length in metres.
    ## Projected area, without fins: 0.1398L2 where L is body length in metres.
    ## > mean(c(.1466,.1398)) [1] 0.1432
    ##
    ## Wetted area, with fins: 0.4631L2 where L is body length in metres.
    ## Wetted area, without fins: 0.4389L2 where L is body length in metres.
    ## > mean(c(.4631,.4389)) # [1] 0.451
    if (type == "projected")
        0.143 * L^2
    else if (type == "wetted")
        0.451 * (0.8715 * L)^2
    else stop("'type' must be 'projected' or 'wetted', not '", type, "' as given")
}

#' Whale compression force
#'
#' Calculate the reaction force of blubber and the sublayer to
#' its interior, as the product of stress and area.
#' Linear compression dynamics is assumed for each layer, i.e. that
#' stress is strain times elastic modulus. Static dynamics is assumed,
#' so that the stresses match in the the two layers. The calculation
#' is done in a two-step process. First, the stress for an equivalent
#' layer is computed, based on a combined elastic modulus. If this
#' stress exceeds the blubber modulus, the blubber thickness is set to zero
#' and the calculation is redone with just the sublayer. If not,
#' the layer thickness are computed assuming equal stress in each layer.
#' Force is calculated as the product of stress and the area
#' given by the product of \code{parms$Ly} and \code{parms$Lz}.
#'
#' @param xs Ship position [m]
#'
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return A list containing: \code{force} [N], the
#' compression-resisting force; \code{stress} [Pa], the ratio
#' of that force to the impact area; \code{strain}, the total
#' strain, and \code{compressed}, a four-column matrix [m]
#' with first column for skin compression, second for blubber
#' compression, third for sublayer compression, and fourth
#' for bone compression.
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleCompressionForce <- function(xs, xw, parms)
{
    zeroTrim <- function(x) # turn negatives into zeros
        ifelse(0 < x, x, 0)
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    ## Note that the denominator of the strain expression vanishes in the stress calculation,
    ## so the next three lines could be simplified. However, retaining it might be clearer,
    ## if a nonlinear stress-strain relationship becomes desirable in the future.
    strain <- dx / parms$lsum
    ##E <- (parms$alpha + parms$beta + parms$gamma) / (parms$alpha/parms$Ealpha + parms$beta/parms$Ebeta + parms$gamma/parms$Egamma)
    stress <- parms$stressFromStrain(strain)
    force <- stress * parms$Ly * parms$Lz
    compressed <- cbind(parms$l[1]*(1-log(1 + stress / parms$a[1]) / parms$b[1]),
                        parms$l[2]*(1-log(1 + stress / parms$a[2]) / parms$b[2]),
                        parms$l[3]*(1-log(1 + stress / parms$a[3]) / parms$b[3]),
                        parms$l[4]*(1-log(1 + stress / parms$a[4]) / parms$b[4]))
    ##. message("compressed=", paste(compressed, " "), "before zero trim")
    compressed <- zeroTrim(compressed)
    ##. message("compressed=", paste(compressed, " "), "after zero trim")
    list(force=force, stress=stress, strain=strain, compressed=compressed)
}

#' Skin force
#'
#' The ship-whale separation is used to calculate the deformation of the skin. The
#' parameters of the calculation are \code{parms$Ly} (impact area width, m),
#' \code{parms$Lz} (impact area height, in m), \code{parms$Ealpha} (skin elastic modulus in Pa),
#' \code{parms$alpha} (skin thickness in m), and \code{parms$theta} (skin bevel angle
#' degrees, measured from a vector normal to undisturbed skin).
#'
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return A list containing \code{force}, the normal force [N], along with
#' \code{sigmay} and \code{sigmaz}, which are stresses [Pa] in the y (beam)
#' and z (draft) directions.
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleSkinForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xs > (xw - parms$lsum)
    dx <- ifelse(touching, xs - (xw - parms$lsum), 0) # penetration distance
    C <- cos(parms$theta * pi / 180) # NB: theta is in deg
    S <- sin(parms$theta * pi / 180) # NB: theta is in deg
    l <- dx * S / C                    # dek20180622_skin_strain eq 1
    s <- dx / C                        # dek20180622_skin_strain eq 2
    ## Strains in y and z
    epsilony <- 2 * (s - l) / (parms$Ly + 2 * l) # dek20180622_skin_strain  eq 3
    epsilonz <- 2 * (s - l) / (parms$Lz + 2 * l) # analogous to dek20180622 eq 3
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
#' the specified mass, \code{ms}, by multiplying by the 2/3
#' power of mass ratio.
#' This is a crude calculation meant as a stop-gap measure, for
#' estimates values of the \code{Ss} argument to \code{\link{parameters}}.
#' It would be much preferable, for a particular simulation, to use the
#' wetted area for a particular ship.
#'
#' @param ms Ship mass [kg].
#'
#' @return Estimated area in m^2.
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
#' where \code{rho} is 1024 (kg/m^3), \code{Cs} is \code{parms$Cs} and
#' \code{A} is \code{parms$Ss}.
#
#' @param vs ship velocity [m/s]
#'
#' @template parmsTemplate
#'
#' @return Water drag force [N]
shipWaterForce <- function(vs, parms)
{
    - (1/2) * 1024 * parms$Cs * parms$Ss * vs * abs(vs)
}


#' Whale force
#'
#' Compute the retarding force of water on the whale, based on a drag law
#' \eqn{(1/2)*rho*Cw*A*vw^2}{(1/2)*rho*Cw*A*vw^2}
#' where \code{rho} is 1024 (kg/m^3), \code{Cw} is \code{parms$Cw} and
#' \code{A} is \code{parms$Sw}.
#'
#' @param vw Whale velocity [m/s]
#'
#' @template parmsTemplate
#'
#' @return Water drag force [N]
whaleWaterForce <- function(vw, parms)
{
    - (1/2) * 1024 * parms$Cw * parms$Sw * vw * abs(vw)
}

#' Dynamical law
#'
#' @param t time [s].
#'
#' @param y model state, a vector containing ship position \code{xs} [m],
#' ship speed \code{vs} [m/s], whale position \code{xw} [m],
#' and whale speed \code{vw} [m/s].
#'
#' @template parmsTemplate
#'
#' @references
#' See \link{whalestrike} for a list of references.
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
    if (is.na(Fship[1])) stop("Fship[1] is NA, probably indicating a programming error.")
    Fwhale <- Freactive + whaleWaterForce(vw, parms)
    if (is.na(Fwhale[1])) stop("Fwhale[1] is NA, probably indicating a programming error.")
    list(c(dxsdt=vs, dvsdt=Fship/parms$ms, dxwdt=vw, dvwdt=Fwhale/parms$mw))
}

#' Calculate derivative using first difference
#' @param var variable.
#' @param t time in seconds.
#' @return Derivative estimated by using \code{\link{diff}} on both \code{x}
#' and \code{t}. In order to return a value of the same length as \code{x} and
#' \code{t}, the last value is repeated.
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
#' \code{\link{shipWaterForce}},
#' \code{\link{whaleSkinForce}},
#' \code{\link{whaleCompressionForce}}, and
#' \code{\link{whaleWaterForce}}.
#'
#' @param t time [s].
#'
#' @param state A named vector holding the initial state of the model:
#' ship position \code{xs} [m],
#' ship speed \code{vs} [m/s],
#' whale position \code{xw} [m]
#' and whale speed \code{vw} [m/s].
#'
#' @template parmsTemplate
#'
#' @template debugTemplate
#'
#' @return An object of class \code{"whalestrike"}, consisting of a
#' list containing vectors for time (\code{t} [s]), ship position (\code{xs} [m]),
#' boat speed (\code{vs} [m/s]), whale position (\code{xw} [m]), whale speed (\code{vw} [m/s]),
#' boat acceleration (\code{dvsdt} [m/s^2]), and whale acceleration (\code{dvwdt} [m/s^2]),
#' a list containing the model parameters (\code{parms}), a list with the results of
#' the skin force calculation (\code{SWF}), a list with the results of the compression
#' force calculations (\code{WCF}), and a vector of whale water force (\code{WWF}).
#'
#' @examples
#' library(whalestrike)
#' t <- seq(0, 0.7, length.out=200)
#' state <- c(xs=-2, vs=10*0.5144, xw=0, vw=0) # 10 knot ship
#' parms <- parameters(ms=20e3, lw=13)
#' sol <- strike(t, state, parms)
#' par(mfcol=c(1, 3), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(sol)
#'
#' @references
#' See \link{whalestrike} for a list of references.
strike <- function(t, state, parms, debug=0)
{
    if (missing(t))
        stop("must supply t")
    if (missing(state))
        stop("must supply state")
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
    res <- list(t=t,
                xs=xs,
                vs=vs,
                xw=xw,
                vw=vw,
                dvsdt=derivative(vs, t),
                dvwdt=derivative(vw, t),
                SWF=shipWaterForce(vs=vs, parms=parms),
                WSF=whaleSkinForce(xs=xs, xw=xw, parms=parms),
                WCF=whaleCompressionForce(xs=xs, xw=xw, parms=parms),
                WWF=whaleWaterForce(vw=vw, parms=parms),
                parms=parms)
    class(res) <- "strike"
    res
}

#' Plot a strike object
#'
#' @description
#' Creates displays of various results of a simulation performed
#' with \code{\link{strike}}.
#'
#' @param x An object inheriting from class \code{strike}
#' @param which A character vector that indicates what to plot.
#' This choices for its entries are listed below, in no particular order.
#' \itemize{
#'
#' \item \code{"location"} for a time-series plot of boat location \code{xw} in
#' dashed black, whale centerline \code{xs} in solid gray,
#' blubber-interior interface in red, and skin in blue. The maximum
#' acceleration of ship and whale (in "g" units) are indicated in notes
#' placed near the horizontal axes.
#'
#' \item \code{"section"} to plot skin thickness, blubber thickness and sublayer thickness
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
#' \item \code{"threat"} a stacked plot showing time-series traces of an
#' ad-hoc measure of the possible threat to skin, blubber and sublayer.
#' The threat level is computed as the ratio
#' of stress to ultimate strength, e.g. for blubber, it is
#' \code{x$WCF$stress/x$parms$s[2]}. Colour-coding indicates levels
#' from 0 to 1/2, from 1/2 to 1, and above 1. An axis categorizes
#' threat level 0 as "L", 0.5 as "M", and 1 as "H".
#'
#' \item \code{"whale acceleration"} for a time-series plot of whale acceleration.
#'
#' \item \code{"blubber thickness"} for a time-series plot of blubber thickness.
#'
#' \item \code{"sublayer thickness"} for a time-series plot of the thickness
#' of the layer interior to the blubber.
#'
#' \item \code{"reactive forces"} for a time-series showing the reactive
#' forces associated with skin stretching (solid) and the compression of the
#' blubber and sublayer components (dashed).
#'
#' \item \code{"compression stress"} for a time-series plot of the compression stress on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#'
#' \item \code{"skin stress"} for a time-series of skin stress in the along-skin y and z directions.
#'
#' \item \code{"values"} for a listing of \code{param} values.
#'
#' \item \code{"all"} for all of the above.
#'
#' \item \code{"default"} for a three-element plot showing \code{"location"},
#' \code{"section"}, and \code{"threat"}.
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
#' @param cols As \code{colw}, but the colour to be used for the ship bow location,
#' which is drawn with a dashed line.
#'
## @param colInjury Two-element colour specification used in \code{"injury"}
## panels. The first colour is used to indicate values that are halfway to the
## injury crition, and the second is used to indicte values that exceed the
## criterion.
#'
#' @param colThreat Three-element colour specification used in \code{"threat"}.
#' The first colour is used for threat levels between 0 and 0.5, the second
#' from 0.5 to 1, and the third above 1. If not provided, the colours
#' will be light-gray, gray, and black, each with an alpha setting that
#' lightens the colours and makes them semi-transparent.
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
#' See \link{whalestrike} for a list of references.
#'
#' @examples
#' ## 1. default 3-panel plot
#' t <- seq(0, 0.7, length.out=200)
#' state <- c(xs=-2, vs=10*0.5144, xw=0, vw=0) # 10 knot ship
#' parms <- parameters(ms=20e3, lw=13)
#' sol <- strike(t, state, parms)
#' par(mar=c(3,3,1,1), mgp=c(2,0.7,0), mfrow=c(1, 3))
#' plot(sol)
#' ## 2. all 12 plot types
#' par(mar=c(3,3,1,1) ,mgp=c(2,0.7,0), mfrow=c(4,3))
#' plot(sol, "all")
plot.strike <- function(x, which="default", drawEvents=TRUE,
                        colwcenter="black", #Slate Gray",
                        colwinterface="black", #colwinterface="Firebrick",
                        colwskin="black", #colwskin="Dodger Blue 4",
                        cols="black",
                        colThreat,
                        lwd=1, D=3, debug=0, ...)
{
    showLegend <- FALSE
    if (missing(colThreat))
        colThreat <- c(rgb(0.827, 0.827, 0.827, alpha=0.7), # lightgray
                       rgb(0.745, 0.745, 0.745, alpha=0.7), # gray
                       rgb(0, 0, 0, alpha=0.7)) # black
    g <- 9.8 # gravity
    t <- x$t
    xs <- x$xs
    vs <- x$vs
    xw <- x$xw
    vw <- x$vw
    dvsdt <- x$dvsdt
    dvwdt <- x$dvwdt
    death <- xs >= xw
    if (any(death)) {
        firstDead <- which(death)[1]
        dead <- firstDead:length(t)
        xw[dead] <- xw[firstDead]
        x$WCF$compressed[1][dead] <- 0
        x$WCF$compressed[2][dead] <- 0
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
                 "blubber thickness", "sublayer thickness", "whale water
                 force", "reactive forces", "skin stress", "compression force",
                 "compression stress", "values")
    for (w in which) {
        if (!(w %in% allowed))
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
        waccel <- derivative(vw, t)
        saccel <- derivative(vs, t)
        mtext(sprintf("ship: %.1fg", max(abs(saccel))/g),
                      side=1, line=-1.25, cex=par("cex"), adj=0.5)
        mtext(sprintf("whale: %.1fg", max(abs(waccel))/g),
                      side=3, line=-1, cex=par("cex"), adj=0.5)
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
        plot(t, -(skin+blubber+sublayer+bone), xlab="Time [s]", ylab="Cross Section [m]",
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
        text(x0, 0.5*(ylim[1] - x$parms$lsum), "water/ship", pos=4)
        hatchPolygon <- FALSE
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

    if (all || "threat" %in% which) {
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
        dylab <- 0.2
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
    if (all || "compression force" %in% which) {
        plot(t, x$WCF$force/1e6, type="l", xlab="Time [s]", ylab="Compress. Force [MN]", lwd=lwd, xaxs="i")
        showEvents(xs, xw)
    }
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
        names <- names(parms)
        values <- unname(unlist(parms))

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
#' @param object An object inheriting from class \code{strike}
#' @param style Character value indicating the representation to be used.
#' If \code{style} is \code{"text"}, then the output is in a textual format. If it
#' is \code{"latex"}, then the output is in latex format.
#'
#' @references
#' See \link{whalestrike} for a list of references.
summarize <- function(object, style="text")
{
    parm <- object$parm
    names <- names(parm)
    ## remove engineForce
    parm <- parm[!(names %in% c("engineForce"))]
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
    if (style == "latex") {
        ## names[which(names == "Ealpha")] <- "E_alpha"
        ## names[which(names == "Ebeta")] <- "E_beta"
        ## names[which(names == "Egamma")] <- "E_gamma"
        ## names[which(names == "UTSalpha")] <- "UTS_alpha"
        ## names[which(names == "UTSbeta")] <- "UTS_beta"
        ## names[which(names == "UTSgamma")] <- "UTS_gamma"
        names[which(names == "lw")] <- "l_w"
        names[which(names == "Ly")] <- "L_y"
        names[which(names == "Lz")] <- "L_z"
        names[which(names == "Ss")] <- "S_s"
        names[which(names == "Sw")] <- "S_w"
        names[which(names == "ms")] <- "M_s"
        names[which(names == "mw")] <- "M_w"
        parm <- data.frame(Symbol=names,
                           Meaning=meaning,
                           Value=as.matrix(unname(unclass(parm))),
                           Reference="")
        colnames(parm) <- c("Item", "Meaning", "Value", "Source")
        print(xtable(parm, label="table:model_parameters",
                     caption="Model parameters.",
                     digits=4,
                     display=c("d", "s", "s", "G", "s")), include.rownames=FALSE)
    } else if (style == "text") {
        parm <- data.frame(Name=names, Meaning=meaning, Value=as.matrix(unname(unclass(parm))))
        colnames(parm) <- c("Item", "Meaning", "Value")
        print(parm, include.rownames=FALSE)
    }
}

