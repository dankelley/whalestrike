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
#' @param delta Whale skin thickness [m]. Defaults to 0.01 m, if not supplied.
#' @param Eskin Whale skin elastic modulus [Pa]. If not provided, defaults to 20e6 Pa,
#' a rounded value of the 19.56+-4.03MPa estimate for seal skin, provided in Table 3
#' of Grear et al. (2017).
#' @param theta Whale skin deformation angle [deg]; defaults to 45deg if not supplied.
#' @param beta Whale blubber thickness [m]; defaults to 0.3m if not supplied.
#' @param Ebeta Elastic modulus of blubber [Pa]; defaults to 6e5 Pa (the value
#' suggested Raymond (2007 fig 37), rounded to 1 digit), if not supplied.
#' @param alpha Thickness of interior region [m]. This defaults to 1m, an
#' estimate of the radius of soft tissue below the blubber. For bone,
#' a value in the centimeter range might be reasonable.
#' @param Ealpha Elastic modulus of interior region [Pa]; defaults to 4e5 Pa,
#' rounded from the soft-tissue value 425294 Pa stated by Raymond (2007) page 36.
#' For bone, the value 9e8 Pa might be reasonable (see Raymond 2007
#' Table 2.3, which lists the elastic modulus of cortical bone as 854.2 MPa).
#' @param Cs Drag coefficient for ship [dimensionless],
#' used by \code{\link{shipWaterForce}} to estimate ship drag force. Defaults
#' to 5e-3, which is two times the frictional coefficient of 2.5e-3
#' inferred from Figure 4 of Manen and van Oossanen (1988), assuming
#' a Reynolds number of 5e7, computed from speed 5m/s, lengthscale 10m
#' and viscosity 1e-6 m^2/s. (The factor of 2 is under the assumption
#' that frictional drag is about half of total drag.)
#' The drag force is computed with \code{\link{shipWaterForce}}.
#' @param Cw Drag coefficient for whale [dimensionless],
#' used by \code{\link{whaleWaterForce}} to estimate whale drag force.
#' Defaults to 3.0e-3, for Reynolds number 2e7, computed from speed
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
                       delta=0.01, Eskin=20e6, theta=45,
                       beta=0.3, Ebeta=6e5,
                       alpha=0.5, Ealpha=4e5,
                       Cs=5e-3, Cw=3e-3)
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
    if (delta < 0)
        stop("whale skin thickness (delta) must be positive, not ", delta)
    if (Eskin < 0)
        stop("whale skin elastic modulus (Eskin) must be positive, not ", Eskin)
    if (theta < 0 || theta > 89)
        stop("whale skin deformation angle (theta) must be between 0 and 89 deg, not ", theta)
    if (beta < 0)
        stop("whale blubber thickness (beta) must be positive, not ", beta)
    if (Ebeta < 0)
        stop("whale blubber elastic modulus (Ebeta) must be positive, not ", Ebeta)
    if (alpha < 0)
        stop("whale sub-blubber thickness (alpha) must be positive, not ", alpha)
    if (Ealpha < 0)
        stop("whale sub-blubber elastic modulus (Ealpha) must be positive, not ", Ealpha)
    if (Cs < 0)
        stop("ship resistance parameter (Cs) must be positive, not ", Cs)
    if (Cw < 0)
        stop("ship resistance parameter (Cw) must be positive, not ", Cw)
    rval <- list(ms=ms, Ss=Ss, impactWidth=impactWidth, impactHeight=impactHeight, mw=mw, Sw=Sw, lw=lw,
                 delta=delta, Eskin=Eskin, theta=theta,
                 Ebeta=Ebeta, beta=beta,
                 Ealpha=Ealpha, alpha=alpha,
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
#' and \code{"surface"}
#' for submerged surface wetted, calculated by spinning
#' the necropsiy side-view presented in Daoust et al. (2018)
#' along the animal axis, yielding 0.737L^2. The \code{"surface"}
#' version is suitable for use in \code{\link{whaleWaterForce}}.
#'
#' @references
#' R worksheet \code{dek/20180623_whale_area.Rmd}, available
#' upon request.
whaleAreaFromLength <- function(L, type="wetted")
{
    ## below from dek/20180623_whale_area.Rmd
    ## Projected area, with fins: 0.1466L2 where L is body length in metres.
    ## Projected area, without fins: 0.1398L2 where L is body length in metres.
    ## Wetted area, with fins: 0.0737L2 where L is body length in metres.
    ## Wetted area, without fins: 0.0698L2 where L is body length in metres.
    if (type == "projected")
        0.143 * L^2
    else if (type == "wetted")
        0.0737 * L^2
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
#' Calculate the reaction force of blubber, as the product of stress
#' and area. Stress is computed as the product of blubber elastic
#' modulus \code{parms$Ebeta} times strain, where the latter is computed
#' from ship-whale separation using \code{1-(xw-xs)/parms$beta} if
#' \code{xw} is between \code{xs} and \code{xs+parms$beta},
#' or zero otherwise. Area is computed as the product of
#' \code{parms$impactWidth} and \code{parms$impactHeight}.
#' FIXME: rewrite docs for the new blubber+inner scheme
#'
#' @param xs Ship position [m]
#'
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return Compression-resisting force of whale blubber and the layer underneath it [N].
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleCompressionForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$beta + parms$alpha)
    dx <- ifelse(touching, parms$alpha + parms$beta - (xw - xs), 0)
    strain <- dx / (parms$beta + parms$alpha)
    E <- (parms$alpha + parms$beta) / (parms$alpha / parms$Ealpha + parms$beta / parms$Ebeta)
    stress <- E * strain
    force <- stress * parms$impactWidth * parms$impactHeight
    ## Assume equal stress in blubber and interior layers
    alphaCompressed <- parms$alpha * (1 - stress / parms$Ealpha)
    betaCompressed <- parms$beta * (1 - stress / parms$Ebeta)
    list(force=force, stress=stress, strain=strain,
         alphaCompressed=alphaCompressed, betaCompressed=betaCompressed)
}

#' Skin force
#'
#' The ship-whale separation is used to calculate the deformation of the skin. The
#' parameters of the calculation are \code{parms$impactWidth} (impact area width, m),
#' \code{parms$impactHeight} (impact area height, in m), \code{parms$Eskin} (skin elastic modulus in Pa),
#' \code{parms$delta} (skin thickness in m), and \code{parms$theta} (skin bevel angle
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
    touching <- xs < xw & xw < (xs + parms$beta + parms$alpha)
    dx <- ifelse(touching, parms$alpha + parms$beta - (xw - xs), 0)
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
    F <- 2*parms$delta*(parms$impactHeight*sigmaz+parms$impactWidth*sigmay)*C # dek20180622_skin_strain eq 8
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
    Fblubber <- whaleCompressionForce(xs, xw, parms)$force
    Fskin <- whaleSkinForce(xs, xw, parms)$force
    Freactive <- Fblubber + Fskin
    Fship <- parms$shipPropulsiveForce + shipWaterForce(vs, parms) - Freactive
    ##. if ((t > 0.1 && t < 0.11) || (t > 0.5 && t < 0.51))
    ##.     cat("t=", t, " vs=", vs, " shipPropulsiveForce=", parms$shipPropulsiveForce, " shipWaterForce=", shipWaterForce(vs, parms), " Freactive=", Freactive, "\n")
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
#' Rw <- 1            # whale radius [m]
#' impactWidth <- 2   # impact region width [m]
#' impactHeight <- 1  # impact region height [m]
#' parms <- parameters(ms=20e3, Ss=ls*(2*draft+beam),
#'               impactWidth=impactWidth, impactHeight=impactHeight,
#'               lw=lw, mw=whaleMassFromLength(lw),
#'               Sw=whaleAreaFromLength(lw, "wetted"),
#'               delta=0.02, Eskin=2e7, theta=45,
#'               beta=0.3, Ebeta=6e5,
#'               alpha=0.5, Ealpha=4e5)
#' sol <- strike(t, state, parms)
#' par(mfcol=c(3, 3), mar=c(2, 3, 0.5, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
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
                   "delta", "Eskin", "theta", # skin properties
                   "beta", "Ebeta", # blubber properties
                   "alpha", "Ealpha")) { # inner-layer properties
        if (!(need %in% names(parms)))
            stop("parms must contain item named '", need, "'; the names you supplied were: ", paste(names(parms), collapse=" "))
    }
    parms["shipPropulsiveForce"] <- -shipWaterForce(state["vs"], parms) # assumed constant over time
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
#' \item \code{"location"} for a time-series plot of boat location \code{xw} in black,
#' whale location \code{x} in blue, and skin location in dashed blue.
#'
#' \item \code{"whale acceleration"} for a time-series plot of whale acceleration,
#' with horizontal lines drawn at 433 and 721 m/s^2,
#' which are estimates of the 10% and 50% probability of concussion
#' for the National Football League concussoin data presented in
#' Figure 9 of Ng et al. 2017.
#'
#' \item \code{"blubber thickness"} for a time-series plot of blubber thickness,
#' with a red line drawn at 2/3 of the
#' initial thickness, based on the Grear et al. 2018 result (page 144, left column)
#' that mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
#' modulus is 1.2MPa, for a ratio of about 2/3.
#'
#' \item \code{"sublayer thickness"} for a time-series plot of the thickness
#' of the layer interior to the blubber. A red line is drawn at 2/3 of the initial
#' thickness, as for \code{"blubber thickness"}.
#'
#' \item \code{"whale water force"} for timea -series of water force on the whale.
#'
#' \item \code{"compression force"} for a time-series of compression force on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#'
#' \item \code{"compression stress"} for a time-series plot of the compression stress on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#'
#' \item \code{"skin force"} for a time-series of normal force resulting from skin tension.
#'
#' \item \code{"skin stress"} for a time-series of skin stress in the along-skin y and z directions.
#' A red horizontal line is drawn at 19.56MPa, if that value is on-scale. This is a suggestion
#' of the maximum stress that the skin can accommodate without damage, inferred from an
#' entry in Table 3 of Grear et al. (2018) for the strength of seal skin.
#'
#' \item \code{"values"} for a listing of \code{param} values.
#'
#' \item \code{"all"} for all of the above.
#'
#' \item \code{"default"} for all except acceleration.
#'}
#' @param center Logical, indicating whether to center time-series plots.
#' on the time when the vessel and whale and in closest proximity.
#' @param drawCriteria Logical, indicating whether to draw coloured horizontal lines
#' on some panels, suggesting critical values (see \dQuote{Details}).
#' @param drawEvents Logical, indicating whether to draw lines for some events,
#' such as the moment of closest approach.
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#' @param ... Ignored.
#'
#' @references
#' See \link{whalestrike} for a list of references.
plot.strike <- function(x, which="default", center=FALSE, drawCriteria=TRUE, drawEvents=TRUE, debug=0, ...)
{
    showLegend <- FALSE
    g <- 9.8 # gravity
    t <- x$t
    xs <- x$xs
    vs <- x$vs
    xw <- x$xw
    vw <- x$vw
    dvsdt <- x$dvsdt
    dvwdt <- x$dvwdt
    if (center) {
        ## Trim so event is centered (maybe; not a big deal)
        end <- 3 * which.min(abs(xs-xw))
        look <- 1:end
        t <- t[look]
        xs <- xs[look]
        vs <- vs[look]
        xw <- xw[look]
        vw <- vw[look]
        dvsdt <- dvsdt[look]
        dvwddt <- dvwdt[look]
    }
    if (debug) {
        cat("t=",  paste(head(t,  3), collapse=" "), " ... ", paste(head(t,  3), collapse=" "), "\n")
        cat("xs=", paste(head(xs, 3), collapse=" "), " ... ", paste(head(xs, 3), collapse=" "), "\n")
        cat("vs=", paste(head(vs, 3), collapse=" "), " ... ", paste(head(vs, 3), collapse=" "), "\n")
        cat("xw=", paste(head(xw, 3), collapse=" "), " ... ", paste(head(xw, 3), collapse=" "), "\n")
        cat("vs=", paste(head(vw, 3), collapse=" "), " ... ", paste(head(vw, 3), collapse=" "), "\n")
    }
    showEvents <- function(xs, xw) {
        if (drawEvents) {
            ##grid()
            death <- which(xs >= xw)[1]
            tdeath <- if (is.finite(death)) t[death] else NA
            if (is.finite(tdeath)) {
                abline(v=tdeath, lwd=2, col="blue")
                mtext("Fatality", at=tdeath, side=3, line=0, col="blue", font=2)
            }
            ## tclosest <- t[which.min(abs(xs-xw))]
            ## abline(v=tclosest, col="darkgreen", lwd=2, lty=3)
        }
    }
    all <- "all" %in% which
    if (which == "default")
        which <- c("location", "blubber thickness", "sublayer thickness",
                   "whale water force", "skin force", "skin stress",
                   "compression force", "compression stress", "values")

    ## Compute some forces, just in case we need them. FIXME: move to strike()

    ## x(t) and xw(t)
    if (all || "location" %in% which) {
        ylim <- range(c(xs, xw), na.rm=TRUE)
        plot(t, xs, type="l", xlab="Time [s]", ylab="Location [m]", ylim=ylim, lwd=2)
        lines(t, xw, col="blue", lwd=2)
        lines(t, xw - x$WCF$alphaCompressed, col="blue", lty="dotted", lwd=0.8)
        lines(t, xw - x$WCF$alphaCompressed - x$WCF$betaCompressed, col="blue", lwd=1)
        accel <- derivative(vw, t)
        mtext(sprintf("max. %.0f m/s^2", max(accel)), side=3, line=-1, cex=par("cex"))
        showEvents(xs, xw)
        if (showLegend)
            legend("topleft", col=c(1, 2, 2), legend=c("Ship", "Whale", "Blubber"),
                   lwd=rep(2, 3), lty=c(1, 1, 3), bg="white", cex=0.8)
    }
    if (all || "whale acceleration" %in% which) {
        plot(t, derivative(vw, t), xlab="Time [s]", ylab="Whale accel. [m/s^2]", type="l", lwd=2)
        showEvents(xs, xw)
        if (drawCriteria) {
            NFL10 <- 433
            NFL50 <- 721
            abline(h=NFL10, col="orange", lwd=2)
            abline(h=NFL50, col="red", lwd=2)
        }
    }
    if (all || "blubber thickness" %in% which) {
        y <- x$WCF$betaCompressed
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Blubber thickness [m]", type="l", lwd=2, ylim=ylim)
        showEvents(xs, xw)
        if (drawCriteria)
            abline(h=x$parms$beta*(1-0.8/1.2), col="red")
    }
    if (all || "sublayer thickness" %in% which) {
        y <- x$WCF$alphaCompressed
        ylim <- c(min(0, min(y)), max(y)) # include 0 if not there by autoscale
        plot(t, y, xlab="Time [s]", ylab="Sublayer thickness [m]", type="l", lwd=2, ylim=ylim)
        showEvents(xs, xw)
        if (drawCriteria)
            abline(h=x$parms$alpha*(1-0.8/1.2), col="red")
    }
    if (all || "whale water force" %in% which) {
        plot(t, whaleWaterForce(vw, x$parms) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin force" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        plot(t, Fs$force/1e6, type="l", xlab="Time [s]", ylab="Skin Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin stress" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        ylim <- range(c(Fs$sigmay, Fs$sigmaz)/1e6)
        plot(t, Fs$sigmay/1e6, type="l", xlab="Time [s]", ylab="Skin Stress [MPa]", lwd=2, ylim=ylim)
        lines(t, Fs$sigmaz/1e6, lty=3)
        ## legend("topright", lty=c(1,3), legend=c("horiz.", "vert."))
        abline(h=19.56, col="red")
        showEvents(xs, xw)
    }
    if (all || "compression force" %in% which) {
        plot(t, x$WCF$force/1e6, type="l", xlab="Time [s]", ylab="Compress. Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "compression stress" %in% which) {
        stress <- x$WCF$force / (x$parms$impactHeight*x$parms$impactWidth)
        plot(t, stress/1e6, type="l", xlab="Time [s]", ylab="Compress. Stress [MPa]", lwd=2)
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
            text(1, i+0.5, values[i], pos=4, cex=0.75)
        par(mar=omar)
    }
}

