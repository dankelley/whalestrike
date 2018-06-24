library(deSolve)

#' whalestrike: A Package to Simulate Ship-Whale Collisions
#'
#' This package solves Newton's second law for a simple model of
#' a ship colliding with a whale. This is a stripped-down model
#' that does not attempt to simulate the biomechanical interactions
#' that are possible in the finite-element treatment of Raymond (2007).
#' Rather, the goal is to establish a convenient framework for rapid
#' computation of responses to a variety of conditions. Indeed, the model runs
#' so quickly that it easily keeps up with mouse movements in an R shiny
#' app, in which the user can test the effects of e.g. changing contact
#' area, ship speed, etc., to build intuition on scenarios ranging
#' from collision with a slow-moving fishing boat to collision
#' with a thin dagger-board of a rapidly-moving racing sailboat.
#' The simple formulation also makes it easy to modify the system
#' by altering various dynamical and biomechanical parameters,
#' and by adding new forces or criteria for whale damage.
#'
#' See the documentation for \code{\link{strike}} for a typical
#' example of using the main functions of this package.
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
#' @param as Ship area [m^2]. This, together with \code{CDs}, is used by
#' used by \code{\link{shipWaterForce}} to estimate ship drag force.
#' @param B  Ship impact horizontal extent [m]; defaults to 2m if not specified.
#' @param D  Ship impact vertical extent [m]; defaults to 1.5m if not specified.
#' @param mw Whale mass [kg]. If not supplied, \code{\link{whaleMassFromLength}} is used
#' to calculate this, given \code{lw}, but if neither \code{mw} nor \code{ml} is provided,
#' an error is reported.
#' @param lw Whale length [m]. If not supplied, \code{\link{whaleLengthFromMass}} is used
#' to calculate this, given \code{lm}, but if neither \code{mw} nor \code{ml} is provided,
#' an error is reported. The length is used by \code{\link{whaleAreaFromLength}} to
#' calculate area, which is needed for the water drag calculation done by
#' \code{\link{whaleWaterForce}}.
#' @param delta Whale skin thickness [m]. Defaults to 0.01 m, if not supplied.
#' @param Es Whale skin elastic modulus [Pa]. If not provided, defaults to 20e6 Pa,
#' a rounded value of the 19.56+-4.03MPa estimate for seal skin, provided in Table 3
#' of Grear et al. (2017).
#' @param theta Whale skin deformation angle [deg]; defaults to 45deg if not supplied.
#' @param beta Whale blubber thickness [m]; defaults to 0.3m if not supplied.
#' @param blubbermodel Blubber model (see \code{\link{whaleBlubberForce}}); defaults to \code{"nonlinear"}
#' if not supplied.
#' @param CDs Drag coefficient for ship [dimensionless], used by \code{\link{shipWaterForce}} to estimate ship drag force.
#' @param CDw Drag coefficient for whale [dimensionless], used by \code{\link{whaleWaterForce}} to estimate whale drag force.
#
#' @return
#' A named list holding the parameters, with defaults and alternatives reconciled
#' according to the system described above.
#'
#' @references
#' See \link{whalestrike} for a list of references.
parameters <- function(ms, as=20, B=3, D=1.5,
                       mw, lw, delta=0.01, Es, theta=45, beta=0.3, blubbermodel="nonlinear",
                       CDs=1, CDw=1)
{
    if (missing(ms))
        stop("ship mass must be specified")
    if (missing(mw)) {
        if (missing(lw))
            stop("Whale mass (mw) or length (lw) must be specified, or both")
        mw <- whaleMassFromLength(lw)
    } else {
        if (missing(lw))
            lw <- whaleLengthFromMass(mw)
    }
    if (missing(Es))
        stop("whale skin elastic modulus (Es) must be specified")
    list(ms=ms, as=as, B=B, D=D, mw=mw, lw=lw, delta=delta, Es=Es, theta=theta, beta=beta, blubbermodel=blubbermodel,
         CDs=CDs, CDw=CDw)
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
        ##exp(2.825-10.095*log(100*L))
    else if (model == "fortune2012pacific")
        exp(-12.286 + 3.158*log(100*L))
        ##exp(3.158-12.286*log(100*L))
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
#' @param view character string indicating the viewpoint; only
#' \code{"side"} is permitted at present, because that is what
#' is needed by \code{\link{whaleWaterForce}}.
whaleAreaFromLength <- function(L, view="side")
{
    if (view == "side")
        0.143 * L^2
    else stop("'view' must be 'side'")
}

#' Contact area
#'
#' Area of contact between vessel and whale skin
#'
#' @param xs Ship position [m]
#'
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return Contact area [m^2].
#' @section DEVELOPMENT NOTE: think about formulation (re theta).
contactArea <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$beta)
    ifelse(touching, parms$B * parms$D, 0)
}

#' Blubber force
#'
#' The stress model depends on \code{parms$blubbermodel}, a character
#' string with the following possibilities:
#' \itemize{
#' \item \code{"Raymond"}: Raymond's (2007) modulus is used, with
#' stress [Pa] taken to be 636273*strain, where strain is the relative
#' compression of the blubber.
#' \item \code{"linear"}: a similar
#' linear law is used, but with coefficient 580004.6, resulting from a linear
#' fit to the data in his figure 2.13, constrained to have zero stress for zero
#' strain.
#' \item \code{"exponential"}: an exponential model is used, with
#' stress=-1.712240e5*(1-exp(strain/0.4185693)). (This is a fit to the
#' data in Raymond's (2007) Figure 2.13.)
#'}
#'
#' @param xs Ship position [m]
#'
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return Compression-resisting force [N]
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleBlubberForce <- function(xs, xw, parms)
{
    ##. message("xs: ", paste(xs, collapse=" "))
    ##. message("xw: ", paste(xw, collapse=" "))
    if (is.na(xs[1])) stop("xs is NA")
    if (is.na(xw[1])) stop("xw is NA")
    if (is.na(parms$beta)) stop("parms$beta is NA")
    touching <- xs < xw & xw < (xs + parms$beta)
    strain <- 1 - (xw - xs) / parms$beta
    ##. message("strain: ", strain)
    ##. ##. if (!is.na(touching[1]) && touching[1]) {
    ##.     cat(sprintf("whaleBlubberForce(): xs %.3f, xw %.3f, B %.3f\n", xs, xw, parms$B))
    ##.     cat(sprintf("whaleBlubberForce(): xs %.3f, xw %.3f, touching %d, strain %.5f\n", xs, xw, touching[1], strain[1]))
    ##. }
    ## nodd <- sum(strain<0) && touching
    ## if (is.na(nodd)) {message('BAD'); browser()}
    ## if (nodd > 0)
    ##     message('x[1]=', x[1], ', xw[1]=', xw[1], ' # negative strains=', nodd)
    ifelse(touching,
           if (parms$blubbermodel=="Raymond") 636273*strain
           else if (parms$blubbermodel=="linear") 580004.6*strain
           else if (parms$blubbermodel=="exponential") -1.712240e+05*(1-exp(strain/4.185693e-01))
           else NA,
           0) * contactArea(xs, xw, parms)
}

#' Skin force
#'
#' The ship-whale separation is used to calculate the deformation of the skin. The
#' parameters of the calculation are \code{parms$B} (ship beam width, m),
#' \code{parms$D} (ship draft, in m), \code{parms$Es} (skin elastic modulus in Pa),
#' \code{parms$delta} (skin thickness in m), and \code{parms$theta} (skin bevel angle
#' degrees, measured from a vector normal to undisturbed skin).
#'
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#'
#' @template parmsTemplate
#'
#' @return A list containing \code{F}, the normal force [N], along with
#' \code{sigmay} and \code{sigmaz}, which are stresses [Pa] in the y (beam)
#' and z (draft) directions.
#'
#' @references
#' See \link{whalestrike} for a list of references.
whaleSkinForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$beta)
    ##> if (is.na(touching[1])) stop("skinForce(): touching is NA")
    dx <- ifelse(touching, parms$beta - (xw - xs), 0)
    l <- dx * tan(parms$theta * pi / 180) # NB: theta is in deg
    s <- dx / cos(parms$theta * pi / 180) # NB: theta is in deg
    ## Strains in y and z
    epsilony <- 2 * (s - l) / (parms$B + 2 * l)
    epsilonz <- 2 * (s - l) / (parms$D + 2 * l)
    ## Stresses in y and z
    sigmay <- parms$E * epsilony
    sigmaz <- parms$E * epsilonz
    ## Net normal force in x; note the cosine, to resolve the force to the normal
    ## direction, and the 2, to account for two sides of length B and two of length D..
    F <- 2 * parms$delta * cos(parms$theta * pi / 180) * (parms$D * sigmaz + parms$B * sigmay)
    ##> if (is.na(F[1])) stop("skinForce(): F is NA")
    list(F=F, sigmay=sigmay, sigmaz=sigmaz)
}

#' Ship water force
#'
#' Estimate the retarding force of water on the ship, based on a drag law
#' \eqn{(1/2)*rho*CD*A*v^2}{(1/2)*rho*CD*A*v^2}
#' where \code{rho} is 1024 (kg/m^3), \code{CD} is \code{parms$CDs} and
#' \code{A} is \code{parms$as}.
#'
#' This function is called by \code{\link{strike}}
#' with the initial ship velocity, to find an engine propulsive force that
#' is assumed to be constant throughout the impact.
#
#' @param vs ship velocity [m/s]
#'
#' @template parmsTemplate
#'
#' @return Water drag force [N]
shipWaterForce <- function(vs, parms)
{
    -0.5 * 1024 * parms$CDs * parms$as * vs * abs(vs)
}


#' Whale force
#'
#' Estimate the retarding force of water on the whale, based on a drag law
#' \eqn{(1/2)*rho*CD*A*v^2}{(1/2)*rho*CD*A*v^2}
#' where \code{rho} is 1024 (kg/m^3), \code{CD} is \code{parms$CDw} and
#' \code{A} is \code{parms$aw}.
#'
#' @param vw Whale velocity [m/s]
#'
#' @template parmsTemplate
#'
#' @return Water drag force [N]
whaleWaterForce <- function(vw, parms)
{
    -0.5 * 1024 * parms$CDw * parms$aw * vw * abs(vw)
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
    Fblubber <- whaleBlubberForce(xs, xw, parms)
    Fskin <- whaleSkinForce(xs, xw, parms)$F
    Freactive <- Fblubber + Fskin
    Fship <- parms$shipPropulsiveForce + shipWaterForce(vs, parms) - Freactive
    Fwhale <- Freactive + whaleWaterForce(vw, parms)
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
#' \code{\link{whaleBlubberForce}}, and
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
#' state <- c(xs=-1.5, vs=5, xw=0, vw=0)
#' parms <- parameters(ms=20e3, B=3, D=1.5,
#'               lw=10, delta=0.02, Es=2e7, theta=45,
#'               beta=0.3, blubbermodel="exponential")
#' sol <- strike(t, state, parms)
#' par(mfcol=c(3, 3), mar=c(2, 3, 1, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(sol, which="all")
#'
#' @references
#' See \link{whalestrike} for a list of references.
strike <- function(t, state, parms, debug=0)
{
    if (missing(t))
        stop("must supply t")
    if (missing(state))
        stop("must supply state")
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
    if (missing(parms))
        stop("must supply parms")
    for (need in c("ms", "B", "D", "mw", "delta", "Es", "theta", "beta", "blubbermodel")) {
        if (!(need %in% names(parms)))
            stop("parms must contain item named '", need, "'; the names you supplied were: ", paste(names(parms), collapse=" "))
    }
    ##> print("in strike(), state is:")
    ##> print(state)
    ##> print("in strike(), parms is:")
    ##> print(parms)
    parms["shipPropulsiveForce"] <- -shipWaterForce(state["vs"], parms) # assumed constant over time
    parms["aw"] <- whaleAreaFromLength(L=parms$lw, view="side")
    sol <- lsoda(state, t, dynamics, parms)
    ##. print("in strike(), head(sol) is:")
    ##. print(head(sol))
    ## Add extra things for plotting convenience. Perhaps should also
    ## calculate the forces here.
    res <- list(t=sol[, 1],
                xs=sol[, 2],
                vs=sol[, 3],
                xw=sol[, 4],
                vw=sol[, 5],
                dvsdt=derivative(sol[, 3], sol[, 1]),
                dvwdt=derivative(sol[, 5], sol[, 1]),
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
#' with horizontal lines drawn at 8g and 15g, which are rough estimates of 50 percent
#' lethality, determined by visual inspection of the lethality-acceleration diagram
#' Figure 9 of Ng et al. 2017, for NFL (National Football League)
#' concussions and combat concussions, respectively.
#'
#' \item \code{"blubber thickness"} for a time-series plot of blubber thickness,
#' with a red line drawn at 2/3 of the
#' initial thickness, based on the Grear et al. 2018 result (page 144, left column)
#' that mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
#' modulus is 1.2MPa, for a ratio of 2/3.
#'
#' \item \code{"whale water force"} for a time-series plot of the water force on the whale.
#'
#' \item \code{"blubber force"} for a time-series plot of the normal force resulting from blubber compression.
#'
#' \item \code{"blubber stress"} for a time-series plot of the normal stress on the blubber.
#'
#' \item \code{"skin force"} for a time-series plot of the normal force resulting from skin tension.
#'
#' \item \code{"skin stress"} for a time-series plot of the skin stress in the along-skin y and z directions.
#' A red horizontal line is drawn at 19.56MPa, if that value is on-scale. This is a suggestion
#' of the maximum stress that the skin can accommodate without damage, inferred from an
#' entry in Table 3 of Grear et al. (2018) for the strength of seal skin.
#'
#' \item \code{"values"} for a listing of \code{param} values.
#'
#' \item \code{"all"} for all of the above.
#'}
#' @param center Logical, indicating whether to center time-series plots.
#' on the time when the vessel and whale and in closest proximity.
#' @param drawCriteria Logical, indicating whether to draw coloured horizontal lines
#' on some panels, suggesting critical values (see \dQuote{Details}).
#' @param drawEvents Logical, indicating whether to draw lines for some events,
#' such as the moment of closest approach.
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#' @param ... Other arguments (ignored).
#'
#' @references
#' See \link{whalestrike} for a list of references.
plot.strike <- function(x, which="all", center=FALSE, drawCriteria=TRUE, drawEvents=TRUE, debug=0, ...)
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

    ## x(t) and xw(t)
    if (all || "location" %in% which) {
        ylim <- range(c(xs, xw), na.rm=TRUE)
        plot(t, xs, type="l", xlab="Time [s]", ylab="Location [m]", ylim=ylim, lwd=2)
        lines(t, xw, col="blue", lwd=2)
        lines(t, xw - x$parms$beta, col="blue", lty=3, lwd=2)
        showEvents(xs, xw)
        if (showLegend)
            legend("topleft", col=c(1, 2, 2), legend=c("Ship", "Whale", "Blubber"),
                   lwd=rep(2, 3), lty=c(1, 1, 3), bg="white", cex=0.8)
    }
    if (all || "whale acceleration" %in% which) {
        plot(t, derivative(vw, t), xlab="Time [s]", ylab="Whale acceleration [m/s^2]", type="l", lwd=2)
        showEvents(xs, xw)
        if (drawCriteria) {
            aCFL <- 8 * 9.8
            abline(h=aCFL, col="orange", lwd=2) # FIXME: determine values more precisely
            ##mtext("CFL", side=4, at=aCFL)
            aConflict <- 15 * 9.8
            abline(h=aConflict, col="red", lwd=2) # FIXME: determine values more precisely
            ##mtext("conflict", side=4, at=aConflict)
        }
    }
    if (all || "blubber thickness" %in% which) {
        touching <- xs < xw & xw < (xs + x$parms$beta)
        thickness <- ifelse(touching, xw-xs, x$parms$beta)
        ylim <- c(0, max(thickness))
        plot(t, thickness, xlab="Time [s]", ylab="Blubber thickness [m]", type="l", lwd=2, ylim=ylim)
        showEvents(xs, xw)
        if (drawCriteria)
            abline(h=x$parms$beta*0.8/1.2, col="red")
    }
    if (all || "whale water force" %in% which) {
        plot(t, whaleWaterForce(vw, x$parms) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin force" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        plot(t, Fs$F/1e6, type="l", xlab="Time [s]", ylab="Skin Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin stress" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        ylim <- range(c(Fs$sigmay, Fs$sigmaz)/1e6)
        plot(t, Fs$sigmay/1e6, type="l", xlab="Time [s]", ylab="Skin Stress [MPa]", lwd=2, ylim=ylim)
        lines(t, Fs$sigmaz/1e6, lty=3)
        legend("topright", lty=c(1,3), legend=c("horiz.", "vert."))
        abline(h=19.56e6, col="red")
        showEvents(xs, xw)
    }
    if (all || "blubber force" %in% which) {
        Fb <- whaleBlubberForce(xs, xw, x$parms)
        plot(t, Fb/1e6, type="l", xlab="Time [s]", ylab="Blubber Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "blubber stress" %in% which) {
        A <- contactArea(xs, xw, x$parms)
        stressb <- ifelse(A, whaleBlubberForce(xs, xw, x$parms) / A, 0)
        plot(t, stressb/1e6, type="l", xlab="Time [s]", ylab="Blubber Stress [MPa]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "values" %in% which) {
        values <- paste(deparse(x$parms),collapse=" ")
        values <- strsplit(gsub(".*\\((.*)\\)$", "\\1", values), ",")[[1]]
        values <- gsub("^ *", "", values)
        values <- gsub("([0-9])L$", "\\1", values) # it's ugly to see e.g. 1 as 1L
        n <- 1 + length(values)
        plot(1:n, 1:n, type="n", xlab="", ylab="", axes=FALSE)
        box()
        for (i in seq_along(values))
            text(1, i+0.5, values[i], pos=4)
    }
}


