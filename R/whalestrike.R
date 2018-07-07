library(deSolve)

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
#' @references
#' \itemize{
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
#'}
#'
#' @docType package
#' @name whalestrike
NULL

#' Control parameters for whale strike simulation
#'
#' Assembles control parameters into a list suitable for passing to \code{\link{strike}}
#' and the functions that it calls.
#'
#' @param ms Ship mass [kg].
#' @param Ss Ship wetted area [m^2]. This, together with \code{Cs}, is used by
#' used by \code{\link{shipWaterForce}} to estimate ship drag force.
#' @param impactWidth Ship impact horizontal extent [m]; defaults to 2m if not specified.
#' @param impactHeight Ship impact vertical extent [m]; defaults to 1.5m if not specified.
#' @param lw Whale length [m]. If not supplied, \code{\link{whaleLengthFromMass}} is used
#' to calculate this, given \code{lm}, but if neither \code{mw} nor \code{ml} is provided,
#' an error is reported. The length is used by \code{\link{whaleAreaFromLength}} to
#' calculate area, which is needed for the water drag calculation done by
#' \code{\link{whaleWaterForce}}.
#' @param mw Whale mass [kg]. If length is known, this could be estimated with
#' \code{\link{whaleMassFromLength}}.
#' @param Sw Whale surface area [m^2]. If length is known, this could be estimated with
#' \code{\link{whaleAreaFromLength}}.
#' @param gamma Whale skin thickness [m]. Defaults to 0.01 m, if not supplied.
#' @param Eskin Whale skin elastic modulus [Pa]. If not provided, defaults to 20e6 Pa,
#' a rounded value of the 19.56+-4.03MPa estimate for seal skin, provided in Table 3
#' of Grear et al. (2017).
#' @param theta Whale skin deformation angle [deg]; defaults to 45deg if not supplied.
#' @param beta Whale blubber thickness [m]; defaults to 0.3m if not supplied.
#' @param Ebeta Elastic modulus of blubber [Pa]; defaults to 0.6 MPa (the value
#' suggested Raymond (2007 fig 37), rounded to 1 digit), if not supplied.
#' @param UTSbeta Numerical value indicating the ultimate tensile strength
#' of the blubber; if not \code{NA}, then this is used to indicate
#' problems in plots made if \code{which} is \code{"compression stress"}.
#' The default is the product of the default whale blubber modulus (see
#' \code{Ebeta}, above) and the strength/modulus ratio for seal
#' blubber, given by (Grear et al. 2018 #' page 144).
#' @param alpha Thickness of interior region [m].
#' The default, 0.5m, may be in an appropriate range for soft tissue;
#' perhaps 0.05m would be more reasonable for bone.
#' @param Ealpha Elastic modulus of interior region [Pa].
#' The default, 0.4MPa, may be in an appropriate range for soft tissue;
#' perhaps 854MPa would be more reasonable for bone.
#' @param UTSalpha Numerical value indicating the ultimate tensile strength
#' of the sublayer.
#' The default, (0.8/1.2)*0.4MPa, may be in an appropriate range for soft tissue;
#' perhaps 22.9MPa would be more reasonable for bone.
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
#'
#' @return
#' A named list holding the parameters, with defaults and alternatives reconciled
#' according to the system described above.
#'
#' @references
#' See \link{whalestrike} for a list of references.
parameters <- function(ms, Ss, impactWidth=3, impactHeight=1.5,
                       lw, mw, Sw,
                       gamma=0.01, Eskin=20e6, theta=45,
                       beta=0.3, Ebeta=0.6e6, UTSbeta=(0.8/1.2)*0.6e6,
                       alpha=0.5, Ealpha=0.4e6, UTSalpha=(0.8/1.2)*0.4e6,
                       ##alpha=0.2, Ealpha=22.9e6, UTSalpha=22.9e6,
                       Cs=0.01, Cw=0.0025)
{
    if (missing(ms) || ms <= 0)
        stop("ship mass (ms) must be given, and a positive number")
    if (missing(Ss) || Ss <= 0)
        stop("ship wetted area (Ss) must be given, and a positive number")
    if (impactWidth <= 0)
        stop("impact width (impactWidth) must be a positive number, not ", impactWidth)
    if (impactHeight <= 0)
        stop("impact height (impactHeight) must be a positive number, not ", impactHeight)
    if (missing(lw))
        stop("Whale length (lm) must be given")
    if (missing(mw))
        stop("Whale mass (mw) must be given; try using whaleMassFromLength() to compute")
    if (missing(Sw))
        stop("Whale surface area (mw) must be given; try using whaleAreaFromLength() to compute")
    if (gamma < 0)
        stop("whale skin thickness (gamma) must be positive, not ", gamma)
    if (Eskin < 0)
        stop("whale skin elastic modulus (Eskin) must be positive, not ", Eskin)
    if (theta < 0 || theta > 89)
        stop("whale skin deformation angle (theta) must be between 0 and 89 deg, not ", theta)
    if (beta < 0)
        stop("whale blubber thickness (beta) must be positive, not ", beta)
    if (Ebeta < 0)
        stop("whale blubber elastic modulus (Ebeta) must be positive, not ", Ebeta)
    if (UTSbeta < 0)
        stop("whale blubber ultimate strength (UTSbeta) must be positive, not ", UTSbeta)
    if (alpha < 0)
        stop("whale sublayer thickness (alpha) must be positive, not ", alpha)
    if (Ealpha < 0)
        stop("whale sublayer elastic modulus (Ealpha) must be positive, not ", Ealpha)
    if (UTSalpha < 0)
        stop("whale sublayer ultimate strength (UTSalpha) must be positive, not ", UTSalpha)
    if (Cs < 0)
        stop("ship resistance parameter (Cs) must be positive, not ", Cs)
    if (Cw < 0)
        stop("ship resistance parameter (Cw) must be positive, not ", Cw)
    rval <- list(ms=ms, Ss=Ss, impactWidth=impactWidth, impactHeight=impactHeight, mw=mw, Sw=Sw, lw=lw,
                 gamma=gamma, Eskin=Eskin, theta=theta,
                 Ebeta=Ebeta, beta=beta, UTSbeta=UTSbeta,
                 Ealpha=Ealpha, alpha=alpha, UTSalpha=UTSalpha,
                 Cs=Cs, Cw=Cw)
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

##OLD #' Contact area
##OLD #'
##OLD #' Area of contact between vessel and whale skin
##OLD #'
##OLD #' @param xs Ship position [m]
##OLD #'
##OLD #' @param xw Whale position [m]
##OLD #'
##OLD #' @template parmsTemplate
##OLD #'
##OLD #' @return Contact area [m^2].
##OLD #' @section DEVELOPMENT NOTE: think about formulation (re theta).
##OLD contactArea <- function(xs, xw, parms)
##OLD {
##OLD     touching <- xs < xw & xw < (xs + parms$beta + parms$alpha)
##OLD     ifelse(touching, parms$impactWidth * parms$impactHeight, 0)
##OLD }

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
#' given by the product of \code{parms$impactWidth} and \code{parms$impactHeight}.
#'
#' @param xs Ship position [m]
#'
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return A list containing \code{force} [N], the
#' compression-resisting force, \code{stress} [Pa], the ratio
#' of that force to the impact area, \code{betaCompressed} [m],
#' the blubber thickness and \code{alphaCompressed} [m],
#' the sublayer thickness.
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleCompressionForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xs > (xw - parms$beta - parms$alpha)
    dx <- ifelse(touching, xs - (xw - parms$alpha - parms$beta), 0) # penetration distance
    strain <- dx / (parms$beta + parms$alpha)
    E <- (parms$alpha + parms$beta) / (parms$alpha / parms$Ealpha + parms$beta / parms$Ebeta)
    stress <- E * strain
    ## Assume equal stress in blubber and interior layers. Do not
    ## permit compression past zero thickness.
    betaTMP <- parms$beta * (1 - stress / parms$Ebeta)
    betaCompressed <- ifelse(0 < betaTMP, betaTMP, 0)
    blubberCrushed <- betaCompressed == 0
    if (any(blubberCrushed)) {         # no blubber left; redo calculation
        strain[blubberCrushed] <- ((dx - parms$beta) / parms$alpha)[blubberCrushed]
        stress[blubberCrushed] <- (parms$Ealpha * strain)[blubberCrushed]
    }
    alphaTMP <- parms$alpha * (1 - stress / parms$Ealpha)
    alphaCompressed <- ifelse(0 < alphaTMP, alphaTMP, 0)
    force <- stress * parms$impactWidth * parms$impactHeight
    list(force=force, stress=stress, strain=strain,
         alphaCompressed=alphaCompressed, betaCompressed=betaCompressed)
}

#' Skin force
#'
#' The ship-whale separation is used to calculate the deformation of the skin. The
#' parameters of the calculation are \code{parms$impactWidth} (impact area width, m),
#' \code{parms$impactHeight} (impact area height, in m), \code{parms$Eskin} (skin elastic modulus in Pa),
#' \code{parms$gamma} (skin thickness in m), and \code{parms$theta} (skin bevel angle
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
    touching <- xs < xw & xs > (xw - parms$beta - parms$alpha)
    dx <- ifelse(touching, xs - (xw - parms$alpha - parms$beta), 0) # penetration distance
    C <- cos(parms$theta * pi / 180) # NB: theta is in deg
    S <- sin(parms$theta * pi / 180) # NB: theta is in deg
    l <- dx * S / C                    # dek20180622_skin_strain eq 1
    s <- dx / C                        # dek20180622_skin_strain eq 2
    ## Strains in y and z
    epsilony <- 2 * (s - l) / (parms$impactWidth + 2 * l) # dek20180622_skin_strain  eq 3
    epsilonz <- 2 * (s - l) / (parms$impactHeight + 2 * l) # analogous to dek20180622 eq 3
    ## Stresses in y and z
    sigmay <- parms$Eskin * epsilony   # dek20180622_skin_strain eq 6
    sigmaz <- parms$Eskin * epsilonz   # dek20180622_skin_strain eq 7
    ## Net normal force in x; note the cosine, to resolve the force to the normal
    ## direction, and the 2, to account for two sides of length
    ## impactWidth and two of length impactHeight.
    F <- 2*parms$gamma*(parms$impactHeight*sigmaz+parms$impactWidth*sigmay)*C # dek20180622_skin_strain eq 8
    list(force=F, sigmay=sigmay, sigmaz=sigmaz)
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
#' along with a list containing the model parameters (\code{parms}).
#'
#' @examples
#' library(whalestrike)
#' t <- seq(0, 1, length.out=500)
#' state <- c(xs=-2, vs=5, xw=0, vw=0)
#' ls <- 10           # ship length [m]
#' draft <- 1.5       # ship draft [m]
#' beam <- 3          # ship beam [m]
#' lw <- 11           # whale length [m]
#' impactWidth <- 0.5 # impact region width [m]
#' impactHeight <- 1  # impact region height [m]
#' parms <- parameters(ms=20e3, Ss=ls*(2*draft+beam),
#'               impactWidth=impactWidth, impactHeight=impactHeight,
#'               lw=lw, mw=whaleMassFromLength(lw),
#'               Sw=whaleAreaFromLength(lw, "wetted"),
#'               gamma=0.02, Eskin=2e7, theta=45,
#'               beta=0.3, Ebeta=6e5,
#'               alpha=0.5, Ealpha=4e5)
#' sol <- strike(t, state, parms)
#' par(mfcol=c(2, 3), mar=c(2, 3, 1, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
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
    for (need in c("ms", "Ss", # ship mass and wetted area
                   "mw", "Sw", # whale mass and wetted area
                   "impactWidth", "impactHeight", # linear extents of impact region
                   "gamma", "Eskin", "theta", # skin properties
                   "beta", "Ebeta", # blubber properties
                   "alpha", "Ealpha")) { # inner-layer properties
        if (!(need %in% names(parms)))
            stop("parms must contain item named '", need, "'; the names you supplied were: ", paste(names(parms), collapse=" "))
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
#' This choices for its entries are:\itemize{
#'
#' \item \code{"location"} for a time-series plot of boat location \code{xw} in
#' dashed black, whale centerline \code{xs} in solid gray,
#' blubber-interior interface in red, and skin in blue.
#'
#' \item \code{"whale acceleration"} for a time-series plot of whale acceleration.
#' If \code{drawCriteria} is \code{TRUE}, the line is thickened during
#' times when the acceleration exceeds 433 m/s^2, which is an estimate of
#' the 10% probability of concussion for the National Football League
#' concussion data presented in Figure 9 of Ng et al. 2017.
#'
#' \item \code{"blubber thickness"} for a time-series plot of blubber thickness.
#' If \code{drawCriteria[1]} is \code{TRUE} then a red line drawn at 2/3 of the
#' initial thickness, based on the Grear et al. 2018 result (page 144, left column)
#' that mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
#' modulus is 1.2MPa, for a ratio of about 2/3.
#'
#' \item \code{"sublayer thickness"} for a time-series plot of the thickness
#' of the layer interior to the blubber.
#' If \code{drawCriteria[1]} is \code{TRUE} then a red line drawn at 2/3 of the
#' initial thickness, based on the Grear et al. 2018 result (page 144, left column)
#' that mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
#' modulus is 1.2MPa, for a ratio of about 2/3. This is really just a provisional
#' guess at a criterion for muscle or organs. (If the simulation resulted from a
#' call to \code{\link{strike}} with \code{Ealpha} corresponding to bone, it is
#' quite unlikely that this criterion will be met.)
#'
#' \item \code{"thickness"} to plot both blubber thickness and sublayer thickness
#' in one panel, creating a cross-section diagram.   If \code{drawCriteria[1]}
#' is \code{TRUE}, then the curve is drawn with a wider pen to
#' indicate times when the blubber is excessively thinned (see
#' \code{"blubber thickness"}).
#' If \code{drawCriteria} is of length 2 and the second
#' element is \code{TRUE}, then the same criterion is used to indicate excessive
#' thinning of the sublayer; this should only be used if the \code{Ealpha} value
#' provided to \code{\link{strike}} represented soft tissue, not bone, and
#' it should be noted that the notion of identical criteria is just a guess.
#'
#' \item \code{"reactive forces"} for a time-series showing the reactive
#' forces associated with skin stretching and the compression of the
#' blubber and sublayer layers undernaeath the skin.
#'
#' \item \code{"compression stress"} for a time-series plot of the compression stress on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#' If \code{drawCriteria[1]} is \code{TRUE} then these traces are thickened if the compression
#' stress exceeds 22.9MPa if the sublayer is bone (Table 2.3 of Raymond, 2007).
#'
#' \item \code{"skin stress"} for a time-series of skin stress in the along-skin y and z directions.
#' If \code{drawCriteria[1]} is \code{TRUE} then these traces are thickened
#' during any times when the stress exceeds 19.5MPa, which is taken as
#' the maximum stress that the skin can accommodate without damage, inferred from an
#' entry in Table 3 of Grear et al. (2018) for the strength of seal skin.
#'
#' \item \code{"values"} for a listing of \code{param} values.
#'
#' \item \code{"all"} for all of the above.
#'
#' \item \code{"default"} for all except acceleration (not needed, since
#' maximal value is in the \code{"location"} panel) and water force (not
#' needed, since it is usually 100s of times smaller than other forces).
#'}
#'
#' @param drawCriteria Logical value indicating whether to
#' indicate dangerous conditions by thickening lines in time series
#' plots. Those conditions are stated in the documentation for individual
#' panels, and for multivariate plots, \code{drawCriteria} can be of
#' length exceeding 1, to control individual curves in the plot.
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
#' @param lwd Line width used in plots for time intervals in which damage
#' criteria are not exceed (or if \code{drawCriteria} is \code{FALSE}).
#'
#' @param D Factor by which to thicken lines during times during which damage
#' criteria are exceeded. Ignored unless \code{drawCriteria} is \code{TRUE}.
#'
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#'
#' @param ... Ignored.
#'
#' @references
#' See \link{whalestrike} for a list of references.
plot.strike <- function(x, which="default", drawCriteria=rep(TRUE, 2), drawEvents=TRUE,
                        colwcenter="Slate Gray", colwinterface="Firebrick", colwskin="Dodger Blue 4",
                        cols="black",
                        lwd=1.4, D=3, debug=0, ...)
{
    showLegend <- FALSE
    if (!is.logical(drawCriteria))
        stop("drawCriteria must be a logical, not ", paste(drawCriteria, collapse=" "), ", as given")
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
        x$WCF$alphaCompressed[dead] <- 0
        x$WCF$betaCompressed[dead] <- 0
    }
    ## t[death] <- NA
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
    if (which == "default")
        which <- c("location", "thickness", "reactive forces", "skin stress",
                   "compression stress", "values")

    ## x(t) and xw(t)
    if (all || "location" %in% which) {
        ylim <- range(c(xs, xw), na.rm=TRUE)
        plot(t, xs, type="l", xlab="Time [s]", ylab="Location [m]", col=cols, ylim=ylim, lwd=lwd, lty=2)
        lines(t, xw, lwd=lwd, col=colwcenter)
        lines(t, xw - x$WCF$alphaCompressed, col=colwinterface, lwd=lwd)
        lines(t, xw - x$WCF$alphaCompressed - x$WCF$betaCompressed, col=colwskin, lwd=lwd)
        accel <- derivative(vw, t)
        mtext(sprintf("Max. accel. %.0f m/s^2 ", max(accel)), side=1, line=-1, cex=par("cex"), adj=1)
        showEvents(xs, xw)
    }
    if (all || "whale acceleration" %in% which) {
        a <- derivative(vw, t)
        plot(t, a, xlab="Time [s]", ylab="Whale accel. [m/s^2]", type="l", lwd=lwd)
        showEvents(xs, xw)
        if (drawCriteria) {
            tt <- t
            tt[a < 433/20] <- NA
            lines(tt, a, lwd=D*lwd)
        }
    }
    if (all || "thickness" %in% which) {
        blubber <- x$WCF$betaCompressed
        sublayer <- x$WCF$alphaCompressed
        maxy <- max(c(blubber+sublayer))
        ylim <- c(0, maxy*1.2) # put y=0 at bottom, so whale-centre is visible
        plot(t, sublayer+blubber, xlab="Time [s]", ylab="Thickness [m]", type="l", lwd=lwd, ylim=ylim, yaxs="i", col=colwskin)
        lines(t, sublayer, col=colwinterface, lwd=lwd)
        abline(h=0, col=colwcenter, lwd=D*lwd)
        showEvents(xs, xw)
        ##text(0, x$parms$alpha + 0.5*x$parms$beta, "Blubber", pos=4)
        ##text(0, 0.5*x$parms$alpha, "Sublayer", pos=4)
        mtext(" skin", side=3, line=-1, adj=0, col=colwskin, cex=par("cex"))
        mtext("interface ", side=3, line=-1, adj=1, col=colwinterface, cex=par("cex"))
        ## legend("bottomright", lwd=lwd, lty=c(1, 3), legend=c("Blubber", "Sublayer"))
        if (drawCriteria[1]) {
            ## Blubber
            tt <- t
            tt[blubber > (1 - 0.8/1.2) * x$parms$beta] <- NA
            lines(tt, sublayer+blubber, lwd=D*lwd, col=colwskin)
        }
        if (length(drawCriteria) > 1 && drawCriteria[2]) {
            ## Sublayer
            tt <- t
            tt[sublayer > (1 - 0.8/1.2) * x$parms$alpha] <- NA
            lines(tt, sublayer, lwd=D*lwd, col=colwinterface)
         }
    }
    if (all || "blubber thickness" %in% which) {
        y <- x$WCF$betaCompressed
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Blubber thickness [m]", type="l", lwd=lwd, ylim=ylim)
        showEvents(xs, xw)
        if (drawCriteria) {
            tt <- t
            tt[blubber > (1 - 0.8/1.2) * x$parms$beta] <- NA
            lines(tt, sublayer+blubber, lwd=D*lwd)
        }
    }
    if (all || "sublayer thickness" %in% which) {
        y <- x$WCF$alphaCompressed
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Sublayer thickness [m]", type="l", lwd=lwd, ylim=ylim)
        showEvents(xs, xw)
        if (drawCriteria) {
            tt <- t
            tt[sublayer > (1 - 0.8/1.2) * x$parms$alpha] <- NA
            lines(tt, sublayer+blubber, lwd=D*lwd)
        }
    }
    if (all || "whale water force" %in% which) {
        plot(t, whaleWaterForce(vw, x$parms) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=lwd)
        showEvents(xs, xw)
    }
    if (all || "reactive forces" %in% which) {
        SF <- x$WSF$force
        CF <- x$WCF$force
        ylim <- range(c(SF/1e6, CF/1e6), na.rm=TRUE)
        plot(t, SF/1e6, type="l", xlab="Time [s]", ylab="Forces [MN]", lwd=lwd, ylim=ylim)
        lines(t, CF/1e6, col=2, lwd=lwd)
        ## legend("topleft", col=1:2, lwd=lwd, legend=c("Skin", "Blubber"))
        mtext(" skin", side=3, line=-1, adj=0, cex=par("cex"))
        mtext("blubber ", side=3, line=-1, adj=1, col=2, cex=par("cex"))
        mtext("& sublayer ", side=3, line=-2, adj=1, col=2, cex=par("cex"))
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
        plot(t, Fs$sigmay/1e6, type="l", xlab="Time [s]", ylab="Skin Stress [MPa]", lwd=lwd, ylim=ylim)
        lines(t, Fs$sigmaz/1e6, col=2, lwd=lwd)
        mtext(" horiz.", side=3, line=-1, adj=0, cex=par("cex"))
        mtext("vert. ", side=3, line=-1, adj=1, col=2, cex=par("cex"))
        if (drawCriteria[1]) {
            tt <- t
            tt[Fs$sigmay < 19.5e6] <- NA
            lines(tt, Fs$sigmay/1e6, lwd=D*lwd)
            tt <- t
            tt[Fs$sigmaz < 19.5e6] <- NA
            lines(tt, Fs$sigmaz/1e6, col=2, lwd=D*lwd)
        }
        showEvents(xs, xw)
    }
    if (all || "compression force" %in% which) {
        plot(t, x$WCF$force/1e6, type="l", xlab="Time [s]", ylab="Compress. Force [MN]", lwd=lwd)
        showEvents(xs, xw)
    }
    if (all || "compression stress" %in% which) {
        force <- x$WCF$force
        stress <- force / (x$parms$impactHeight*x$parms$impactWidth)
        plot(t, stress/1e6, type="l", xlab="Time [s]", ylab="Compress. Stress [MPa]", lwd=lwd)
        if (drawCriteria[1] && !is.na(x$parms$UTSalpha) && !is.na(x$parms$UTSbeta)) {
            tt <- t
            tt[stress < min(x$parms$UTSalpha, x$parms$UTSbeta)] <- NA
            lines(tt, stress/1e6, lwd=D*lwd)
        }
        showEvents(xs, xw)
    }
    ## FIXME: add items for sub-blubber force and stress
    if (all || "values" %in% which) {
        omar <- par("mar")
        par(mar=rep(0, 4))
        p <- lapply(x$parms, function(x) signif(x, 4))
        names(p) <- names(x$parms)
        values <- paste(deparse(p), collapse=" ")
        values <- strsplit(gsub(".*\\((.*)\\)$", "\\1", values), ",")[[1]]
        values <- gsub("^ *", "", values)
        values <- gsub("([0-9])L$", "\\1", values) # it's ugly to see e.g. 1 as 1L
        n <- 1 + length(values)
        plot(1:n, 1:n, type="n", xlab="", ylab="", axes=FALSE)
        for (i in seq_along(values))
            text(1, i+0.5, values[i], pos=4, cex=1)
        par(mar=omar)
    }
}

