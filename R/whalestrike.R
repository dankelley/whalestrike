# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

library(deSolve)

#' Reference strike() solution
#'
#' This was produced with the package as it existed on 2020-jul-8,
#' prior to the publication of Kelley et al. (2021).  It is used
#' in testing, to ensure that the package does not inadvertently
#' change in its predictions.
#'
#' @examples
#' data(sol20200708)
#' opar <- par(mfrow = c(1, 3))
#' plot(sol20200708)
#' par(opar)
#'
#' @references
#'
#' Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the Lethality
#' of Ship Strikes on Whales Using Simple Biophysical Models." Marine Mammal
#' Science, 37(1), 2021, mms.12745. \doi{10.1111/mms.12745}.
#'
#' @name sol20200708
#'
#' @docType data
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
#' @references
#'
#' Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' "Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine Animal Response
#' Society, and Fisheries and Oceans Canada, 2018.
#' \url{https://publications.gc.ca/site/eng/9.850838/publication.html}.
#'
#' @examples
#' library(whalestrike)
#' shape <- whaleShape()
#' plot(shape$x, shape$y, asp = 1, type = "l")
#' polygon(shape$x, shape$y, col = "lightgray")
#' lw <- 13.7
#' Rmax <- 0.5 * lw * diff(range(shape$y))
#' mtext(sprintf("Max. radius %.2fm for %.1fm-long whale", Rmax, lw), side = 3)
#'
#' @export
#'
#' @return `whaleShape` returns a data frame with entries named `x` and `y`, which trace
#' out the side-view shape of a Right Whale, using values digitized from a diagram
#' in Daoust et al. (2017).
#'
#' @author Dan Kelley
whaleShape <- function() {
    structure(list(x = c(
        0, 0.00790510432232446, 0.0316204955297698,
        0.0592882041769615, 0.104084003978059, 0.134387321162819, 0.160737147300755,
        0.191040464485515, 0.20685051664922, 0.225295707907663, 0.262186872829267,
        0.31225216841512, 0.372858020379922, 0.429511946107337, 0.462449619982116,
        0.496705645808983, 0.546770941394835, 0.607376793359637, 0.685111049440232,
        0.762845305520828, 0.840577214387265, 0.919631387229387, 0.96178735350743,
        1, 0.963109617483115, 0.920945827157878, 0.882741004712502, 0.844528358219932,
        0.772069857161847, 0.731224418788698, 0.683793479892864, 0.62450519747543,
        0.553358789131679, 0.496705645808983, 0.459814480887378, 0.403161337564682,
        0.334650068315668, 0.276679355445603, 0.21739107302817, 0.166007425490229,
        0.115942129904377, 0.0632409128190676, 0.0197627608058111, 0.00263498261379353,
        0.00131749130689676
    ), y = c(
        0, 0.0184404968301264, 0.030294710732848,
        0.0461008508729563, 0.0645413477030828, 0.075078774463155, 0.0790307007005417,
        0.0790307007005417, 0.0856162012232272, 0.0961536279832994, 0.101422341363336,
        0.10800862429074, 0.10932541143339, 0.111959768123408, 0.110642980980758,
        0.106691054743372, 0.101422341363336, 0.0948368408406501, 0.0816650573905598,
        0.0658589172504515, 0.0461008508729563, 0.0263435669001806, 0.0144885705927397,
        0.00395114383266745, -0.00658628292740479, -0.0105374267600722,
        -0.0184404968301264, -0.0276611364475493, -0.0487359899676937,
        -0.059273416727766, -0.0698108434878382, -0.0777139135578923,
        -0.0882513403179646, -0.0974711975306681, -0.101423123768055,
        -0.100105554220686, -0.0908856970079826, -0.0856169836279465,
        -0.0763963440105236, -0.0671764867978201, -0.061907773417784,
        -0.0553214904903792, -0.0461016332776756, -0.0276611364475493,
        -0.00526871338003609
    )), class = "data.frame", row.names = c(
        NA,
        -45L
    ))
}



#' whalestrike: A Package to Simulate Ship-Whale Collisions
#'
#' This package solves Newton's second law for a simple model of
#' a ship colliding with a whale. This is a stripped-down model
#' that does not attempt to simulate the biomechanical interactions
#' that can be simulated in finite-element treatments such
#' as that of Raymond (2007).  For an in-depth discussion
#' of the reason for writing the model, of the principles involved
#' in its framing, and its use in developing a criterion for
#' strike lethality, please see Kelley et al. (2021).
#'
#' The goal of the model is to establish a
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
#' Kelley et al (2021) provide more
#' detail about the mathematical
#' framework of the package, along with a discussion of its
#' purpose and application to real-world problems of ship
#' strikes on whales.
#'
#' @section Further reading:
#' \itemize{
#'
#' \item
#'
#' Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' "Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine Animal Response
#' Society, and Fisheries and Oceans Canada, 2018.
#' \url{https://publications.gc.ca/site/eng/9.850838/publication.html}.
#'
#' \item
#' Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. "Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis)." Journal of Mammalogy 93,
#' no. 5 (2012): 1342-54. \doi{10.1644/11-MAMM-A-297.1}.
#'
#' \item
#' Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt, Adam
#' P. Summers, and Petra Ditsche. "Mechanical Properties of Harbor Seal Skin and
#' Blubber--a Test of Anisotropy." Zoology 126 (2018): 137-44.
#' \doi{10.1016/j.zool.2017.11.002}.
#'
#' \item
#' Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the
#' Lethality of Ship Strikes on Whales Using Simple Biophysical Models."
#' Marine Mammal Science 37, no. 1 (January 2021): 251–67.
#' \doi{10.1111/mms.12745}.
#'
#' \item
#' Kelley, Dan E. "Composite Spring," May 28, 2018. 20180528_composite_string.
#' Dan Kelley's working notes.
#'
#' \item
#' Kelley, Dan. "Whale Area," June 23, 2018. 20180623_whale_area.
#' Dan Kelley's working notes.
#'
#' \item
#' Kelley, Dan. "Ship Propulsion," July 1, 2018. 20180701_ship_propulsion.
#' Dan Kelley's working notes.
#'
#' \item
#' Kelley, Dan. "Whale Mass," July 7, 2018. 20180707_whale_mass. Dan Kelley's working notes.
#'
#' \item
#' Kelley, Dan E."“Whalestrike: An R Package for Simulating Ship Strikes on Whales."
#' Journal of Open Source Software 9, no. 97 (2024): 6473.
#' https://doi.org/10.21105/joss.06473.
#'
#' \item
#' MAN Diesel & Turbo. "Basic Principles of Propulsion." MAN Diesel & Turbo, 2011.
# nolint start line_length_linter
#' \code{https://spain.mandieselturbo.com/docs/librariesprovider10/sistemas-propulsivos-marinos/basic-principles-of-ship-propulsion.pdf?sfvrsn=2}
# nolint end line_length_linter
#'
#' \item
#' Manen, J. D. van, and P. van Oossanen. "Resistance." In Principles of Naval
#' Architecture (Second Revision), Volume II - Resistance, Propulsion and
#' Vibration, edited by Edward V Lewis, Second Edition, 1-125. Jersey City, NJ: Society
#' of Naval Architects and Marine Engineers (U.S.), 1988.
#'
#' \item
#' Mayette, Alexandra. "Whale Layer Thickness." December 15, 2025. (Personal
#' communication of a 5-page document.)
#'
#' \item
#' Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method to Estimate
#' Vessel Mass for Use in Whale-Ship Strike Risk Models." PloS One 21, no. 1 (2026):
#' e0339760. https://doi.org/10.1371/journal.pone.0339760.
#'
#' \item
#' Miller, Carolyn A., Desray Reeb, Peter B. Best, Amy R. Knowlton, Moira W.
#' Brown, and Michael J. Moore. "Blubber Thickness in Right Whales Eubalaena
#' Glacialis and Eubalaena Australis Related with Reproduction, Life History
#' Status and Prey Abundance." Marine Ecology Progress Series 438 (2011): 267-83.
#'
#' \item
#' Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
#' "Morphometry, Gross Morphology and Available Histopathology in North Atlantic
#' Right Whale (Eubalaena Glacialis) Mortalities (1970 to 2002)." Journal of
#' Cetacean Research and Management 6, no. 3 (2005): 199-214.
#'
#' \item
#' Ng, Laurel J., Vladislav Volman, Melissa M. Gibbons, Pi Phohomsiri, Jianxia
#' Cui, Darrell J. Swenson, and James H. Stuhmiller. "A Mechanistic End-to-End
#' Concussion Model That Translates Head Kinematics to Neurologic Injury."
#' Frontiers in Neurology 8, no. JUN (2017): 1-18.
#' \doi{10.3389/fneur.2017.00269}
#'
#' \item
#' Raymond, J. J. "Development of a Numerical Model to Predict Impact Forces on a
#' North Atlantic Right Whale during Collision with a Vessel." University of New
#' Hampshire, 2007.
#' \url{https://scholars.unh.edu/thesis/309/}.
#'
#' \item
#' Soetaert, Karline, Thomas Petzoldt, and R. Woodrow Setzer.
#' "Solving Differential Equations in R: Package DeSolve."
#' Journal of Statistical Software; Vol 1, Issue 9, 2010.
#' \doi{10.18637/jss.v033.i09}.
#'
#' }
#'
#' @name whalestrike
#' @docType package
#' @keywords internal
"_PACKAGE"
NULL

#' Whale projected area, as function of length
#'
#' This depends on calculations based on the digitized shape of
#' a whale necropsy, which is provided by [whaleShape()].
#' The results are
#' \eqn{0.143 * L^2}{0.143 * L^2}
#' for the projected area (see reference 1)
#' and
#' \eqn{0.448 * (0.877 * L)^2}{0.448 * (0.877 * L)^2}
#' for the wetted area
#' (see reference 2, but note that we use a correction related to whale mass).
#'
#' Note that multiple digitizations were done,
#' and that the coefficients used in the formulae
#' agreed to under 0.7 percent percent between these digitizations.
#'
#' @param L whale length in metres.
#'
#' @param species a string indicating the whale species. In
#' the present version of the package, this parameter is ignored,
#' and it is assumed that the formula developed for North
#' Atlantic Right Whales will be applicable to other species.
#' This is not a large concern, because the area only affects
#' the water drag, which will not be large during the
#' short interval of a ship impact.
#'
#' @param type character string indicating the type of area, with
#' `"projected"` for a side-projected area, and
#' `"wetted"` for the total wetted area. The wetted
#' area was computed by mathematically spinning a spline fit to the
#' side-view. In both cases, the original data source is the
#' necropsy side-view presented in Daoust et al. (2018).
#'
#' @examples
#' L <- 3:20
#' A <- whaleAreaFromLength(L)
#' plot(L, A, xlab = "Length [m]", ylab = "Area [m^2]", type = "l")
#'
#' @references
#' 1. Dan Kelley's internal document `dek/20180623_whale_area.Rmd`, available
#' upon request.
#'
#' 2. Dan Kelley's internal document `dek/20180707_whale_mass.Rmd`, available
#' upon request.
#'
#' 3. Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' "Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine Animal Response
#' Society, and Fisheries and Oceans Canada, 2018.
#' \url{https://publications.gc.ca/site/eng/9.850838/publication.html}.
#' @author Dan Kelley
#'
#' @return `whaleAreaFromLength` returns the surface area of the whale,
#' in square metres.
#'
#' @export
whaleAreaFromLength <- function(L, species = "N. Atl. Right Whale", type = "wetted") {
    #<uncomment later, perhaps> speciesAllowed <- c("N. Atl. Right Whale")
    #<uncomment later, perhaps> if (!(species %in% speciesAllowed)) {
    #<uncomment later, perhaps>     stop(
    #<uncomment later, perhaps>         "unknown species \"", species, "\"; use one of the following: \"",
    #<uncomment later, perhaps>         paste(speciesAllowed, collapse = "\", \""), "\""
    #<uncomment later, perhaps>     )
    #<uncomment later, perhaps> }
    # below from dek/20180623_whale_area.Rmd, updated 20180802.
    #
    # * Projected area, with fins: 0.1466 * L^2 where L is body length in metres.
    # * Projected area, without fins: 0.1391 * L^2 where L is body length in metres.
    # * Wetted area, with fins: 0.4631 * L^2 where L is body length in metres.
    # * Wetted area, without fins: 0.4336 * L^2 where L is body length in metres.
    #
    # The relevant case (with or without fins) being dependent on the application,
    # there may be merit in averaging the two estimates, yielding:
    # * Projected area: 0.1429 * L^2 where L is body length in metres.
    # * Wetted area: 0.4484 * L^2 where L is body length in metres.
    if (type == "projected") {
        0.143 * L^2
    } else if (type == "wetted") {
        0.448 * (0.877 * L)^2
    } else {
        stop("'type' must be 'projected' or 'wetted', but it is '", type, "'")
    }
}

#' Dynamical law
#'
#' This function handles Newton's second law, which is the dynamical
#' law that relates the accelerations of whale and ship to the forces
#' upon each.  It is used by [strike()], as the latter integrates
#' the acceleration equations to step forward in time through
#' the simulation of a whale-strike event. Thus, [dynamics()]
#' is a core function of this package.  The code is very simple,
#' because the forces are determined by other functions, as
#' described in the \dQuote{Details} section.
#'
#' Given a present state (defined by the positions and
#' velocities of ship and whale) at the present time,
#' apply Newton's second law to find the time derivatives
#' of that state.  Forces are determined with
#' [whaleCompressionForce()],
#' [whaleSkinForce()],
#' [shipWaterForce()],
#' [whaleWaterForce()], while engine force
#' (assumed constant over the course of a collision) is
#' computed from initial [shipWaterForce()].  Whale and
#' ship masses are set by [parameters()], which also sets up
#' areas, drag coefficients, etc.
#'
#' @param t time (s).
#'
#' @param y model state, a vector containing ship position `xs` (m),
#' ship speed `vs` (m/s), whale position `xw` (m),
#' and whale speed `vw` (m/s).
#'
#' @template parmsTemplate

#' @return An List contain items named `dxsdt` (time derivative of
#' ship location), `dvsdt` (time derivative of ship speed),
#' `dxwdt` (time derivative of whale location) and `dvwdt`
#' (time derivative of whale speed). These are computed by
#' solving the dynamical system using Newton's second law,
#' based on the known masses of ship and whale, and the forces
#' involved in the collision.
#
#  list(c(dxsdt = vs, dvsdt = Fship / parms$ms,
#  dxwdt = vw, dvwdt = Fwhale / parms$mw))
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
dynamics <- function(t, y, parms) {
    xs <- y[1] # ship position
    vs <- y[2] # ship velocity
    xw <- y[3] # whale position
    vw <- y[4] # whale velocity
    Fcompression <- whaleCompressionForce(xs, xw, parms)$force
    Fextension <- whaleSkinForce(xs, xw, parms)$force
    Freactive <- Fcompression + Fextension
    Fship <- if (!is.null(parms$engineForce)) {
        parms$engineForce + shipWaterForce(vs, parms) - Freactive
    } else {
        -Freactive
    }
    if (is.na(Fship[1])) {
        stop("Fship[1] is NA, probably indicating a programming error.")
    }
    Fwhale <- Freactive + whaleWaterForce(vw, parms)
    # cat("t=", t, ", whaleWaterForce=", whaleWaterForce(vw, parms), "\n")
    if (is.na(Fwhale[1])) {
        stop("Fwhale[1] is NA, probably indicating a programming error.")
    }
    list(c(dxsdt = vs, dvsdt = Fship / parms$ms, dxwdt = vw, dvwdt = Fwhale / parms$mw))
}

#' Simulate the collision of a ship and a whale
#'
#' Newtonian mechanics are used, taking the ship as non-deformable,
#' and the whale as being cushioned by a skin layer and a blubber layer.
#' The forces are calculated by
#' [shipWaterForce()],
#' [whaleSkinForce()],
#' [whaleCompressionForce()], and
#' [whaleWaterForce()] and the integration is carried out with
#' [deSolve::lsoda()].
#'
#' @param t a suggested vector of times (s) at which the simulated state will be reported.
#' This is only a suggestion, however, because `strike` is set up to detect high
#' accelerations caused by bone compression, and may set a finer reporting interval,
#' if such accelerations are detected. The detection is based on thickness of
#' compressed blubber and sublayer; if either gets below 1 percent
#' of the initial(uncompressed) value, then
#' a trial time grid is computed, with 20 points during the timescale for
#' bone compression, calculated as
#' \eqn{0.5*sqrt(Ly*Lz*a[4]*b[4]/(l[4]*mw)}{0.5*sqrt(Ly*Lz*a[4]*b[4]/(l[4]*mw)},
#' with terms as discussed in
#' the documentation for [parameters()]. If this trial grid is finer
#' than the grid in the `t` parameter, then the simulation is redone
#' using the new grid. Note that this means that the output will
#' be finer, so code should not rely on the output time grid being
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
#' t <- seq(0, 0.7, length.out = 200)
#' state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0) # ship speed 10 knots
#' parms <- parameters()
#' sol <- strike(t, state, parms)
#' opar <- par(mfcol = c(1, 3), mar = c(3, 3, 0.5, 2), mgp = c(2, 0.7, 0), cex = 0.7)
#' plot(sol)
#' par(opar)
#'
#' # Example 2: time-series plots of blubber stress and stress/strength,
#' # for a 200 tonne ship moving at 10 knots
#' t <- seq(0, 0.7, length.out = 1000)
#' state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0) # ship 10 knots
#' parms <- parameters(ms = 200 * 1000) # 1 metric tonne is 1000 kg
#' sol <- strike(t, state, parms)
#' opar <- par(mfrow = c(2, 1), mar = c(3, 3, 0.5, 2), mgp = c(2, 0.7, 0), cex = 0.7)
#' plot(t, sol$WCF$stress / 1e6,
#'     type = "l",
#'     xlab = "Time [s]", ylab = "Blubber stress [MPa]"
#' )
#' plot(t, sol$WCF$stress / sol$parms$s[2],
#'     type = "l",
#'     xlab = "Time [s]", ylab = "Blubber stress / strength"
#' )
#' par(opar)
#'
#' # Example 3: max stress and stress/strength, for a 200 tonne ship
#' # moving at various speeds. This is a slow calculation, so we do
#' # not run it by default.
#' \dontrun{
#' knots <- seq(0, 20, 0.5)
#' maxStress <- NULL
#' maxStressOverStrength <- NULL
#' for (speed in knot2mps(knots)) {
#'     t <- seq(0, 10, length.out = 1000)
#'     state <- list(xs = -2, vs = speed, xw = 0, vw = 0)
#'     parms <- parameters(ms = 200 * 1000) # 1 metric tonne is 1000 kg
#'     sol <- strike(t, state, parms)
#'     maxStress <- c(maxStress, max(sol$WCF$stress))
#'     maxStressOverStrength <- c(maxStressOverStrength, max(sol$WCF$stress / sol$parms$s[2]))
#' }
#' opar <- par(mfrow = c(2, 1), mar = c(3, 3, 0.5, 2), mgp = c(2, 0.7, 0), cex = 0.7)
#' nonzero <- maxStress > 0
#' plot(knots[nonzero], log10(maxStress[nonzero]),
#'     type = "o", pch = 20, xaxs = "i", yaxs = "i",
#'     xlab = "Ship Speed [knots]", ylab = "log10 peak blubber stress"
#' )
#' abline(h = log10(sol$parms$s[2]), lty = 2)
#' plot(knots[nonzero], log10(maxStressOverStrength[nonzero]),
#'     type = "o", pch = 20, xaxs = "i", yaxs = "i",
#'     xlab = "Ship Speed [knots]", ylab = "log10 peak blubber stress / strength"
#' )
#' abline(h = 0, lty = 2)
#' par(opar)
#' }
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom deSolve lsoda
strike <- function(t, state, parms, debug = 0) {
    if (missing(t)) {
        stop("must supply t")
    }
    # Ensure that the state is well-configured, because the error messages
    # otherwise will be too cryptic for many users to fathom.
    if (missing(state)) {
        stop("must supply state")
    }
    if (4 != sum(c("xs", "vs", "xw", "vw") %in% names(state))) {
        stop("state must hold \"xs\", \"vs\", \"xw\", and \"vw\"")
    }
    if (is.list(state)) {
        state <- c(xs = state$xs, vs = state$vs, xw = state$xw, vw = state$vw)
    }
    if (missing(parms)) {
        stop("must supply parms")
    }
    if (!inherits(parms, "parameters")) {
        stop("parms must be the output of parameters()")
    }
    # Check parameters
    parmsRequired <- c(
        "a", "b", "Cs", "Cw", "l", "lsum", "lw", "Ly", "Lz",
        "ms", "mw", "s", "Ss", "stressFromStrain", "Sw",
        "theta"
    )
    # cat("parms$Sw=", parms$Sw, "\n")
    if (!all(parmsRequired %in% names(parms))) {
        stop('parms must hold: "', paste(parmsRequired, collapse = '", "'), '"')
    }
    # All required elements are present, but it's prudent to check some values that
    # a user might be setting.
    if (!is.finite(parms$ms) || parms$ms <= 0) {
        stop("parms$ms (ship mass, in kg) must be a positive number, not ", parms$ms)
    }
    if (!is.finite(parms$mw) || parms$mw <= 0) {
        stop("parms$mw (whale mass, in kg) must be a positive number, not ", parms$mw)
    }
    if (!is.finite(parms$Ly) || parms$Ly <= 0) {
        stop("parms$Ly (impact width, in m) must be a positive number, not ", parms$Ly)
    }
    if (!is.finite(parms$Lz) || parms$Lz <= 0) {
        stop("parms$Lz (impact height, in m) must be a positive number, not ", parms$Lz)
    }
    if (!is.finite(parms$Cs) || parms$Cs <= 0) {
        stop("parms$Cs (drag coefficient of ship, dimensionless) must be a positive number, not ", parms$Cs)
    }
    if (!is.finite(parms$Cw) || parms$Cw <= 0) {
        stop("parms$Cw (drag coefficient of whale, dimensionless) must be a positive number, not ", parms$Cw)
    }
    if (!is.function(parms$stressFromStrain)) {
        stop("parms$stressFromStrain must be a function")
    }
    if (debug > 0) {
        print("state:")
        print(state)
        print("parms:")
        print(parms)
    }
    for (need in c("xs", "vs", "xw", "vw")) {
        if (!(need %in% names(state))) {
            stop(
                "state must contain item named '", need, "'; the names you supplied were: ",
                paste(names(state), collapse = " ")
            )
        }
    }
    parms["engineForce"] <- -shipWaterForce(state["vs"], parms) # assumed constant over time
    sol <- lsoda(state, t, dynamics, parms)
    # Add extra things for plotting convenience.
    t <- sol[, 1]
    xs <- sol[, 2]
    vs <- sol[, 3]
    xw <- sol[, 4]
    vw <- sol[, 5]
    dvsdt <- derivative(vs, t)
    dvwdt <- derivative(vw, t)
    SWF <- shipWaterForce(vs = vs, parms = parms)
    WSF <- whaleSkinForce(xs = xs, xw = xw, parms = parms)
    WCF <- whaleCompressionForce(xs = xs, xw = xw, parms = parms)
    WWF <- whaleWaterForce(vw = vw, parms = parms)
    # nolint start line_length_linter
    refinedGrid <- (min(WCF$compressed[, 2], na.rm = TRUE) / WCF$compressed[1, 2] < 0.01) || (min(WCF$compressed[, 3], na.rm = TRUE) / WCF$compressed[1, 3] < 0.01)
    # nolint end line_length_linter
    if (refinedGrid) {
        NEED <- 20 # desired number of points in peak
        dt <- (1 / NEED) * 0.5 * sqrt(parms$l[4] * parms$mw / (parms$Ly * parms$Lz * parms$a[4] * parms$b[4]))
        tstart <- t[1]
        tend <- tail(t, 1)
        nold <- length(t)
        n <- floor(0.5 * (tend - tstart) / dt)
        if (n > length(t)) {
            warning("increasing from ", nold, " to ", n, " time steps, to capture acceleration peak\n")
            t <- seq(tstart, tend, length.out = n)
            sol <- lsoda(state, t, dynamics, parms)
            # Add extra things for plotting convenience.
            t <- sol[, 1]
            xs <- sol[, 2]
            vs <- sol[, 3]
            xw <- sol[, 4]
            vw <- sol[, 5]
            dvsdt <- derivative(vs, t)
            dvwdt <- derivative(vw, t)
            SWF <- shipWaterForce(vs = vs, parms = parms)
            WSF <- whaleSkinForce(xs = xs, xw = xw, parms = parms)
            WCF <- whaleCompressionForce(xs = xs, xw = xw, parms = parms)
            WWF <- whaleWaterForce(vw = vw, parms = parms)
        }
    }
    res <- list(
        t = t, xs = xs, vs = vs, xw = xw, vw = vw,
        dvsdt = dvsdt, dvwdt = dvwdt,
        SWF = SWF, WSF = WSF, WCF = WCF, WWF = WWF,
        parms = parms,
        refinedGrid = refinedGrid
    )
    class(res) <- "strike"
    res
}

#' Summarize the parameters of a simulation, and its results
#'
#' @param object an object of class `"strike"`, as created by [strike()].
#'
#' @param \dots ignored
#'
#' @examples
#' library(whalestrike)
#' # Example 1: graphs, as in the shiny app
#' t <- seq(0, 0.7, length.out = 200)
#' state <- list(xs = -2, vs = knot2mps(10), xw = 0, vw = 0) # ship speed 10 knots
#' parms <- parameters()
#' sol <- strike(t, state, parms)
#' summary(sol)
#'
#' @author Dan Kelley
#'
#' @return None. This function is called for its side-effect of
#' printing information about the ship-whale collision simulation.
#'
#' @export
summary.strike <- function(object, ...) {
    parm <- object$parm
    summary(parm)
    LI <- lethalityIndexFromStress(object$WCF$stress)
    peakLI <- which.max(LI)
    cat("\nSimulation results returned by strike()\n")
    cat(sprintf("  simulated time range: 0 to %g s\n", max(object$t)))
    cat(sprintf("  xs: %12g m        -- ship position at t=0 s\n", object$xs[1]))
    cat(sprintf("  vs: %12.3g m/s      -- ship speed at t=0 s\n", object$vs[1]))
    cat(sprintf("      %12.3g knot     -- above, in a nautical unit\n", mps2knot(object$vs[1])))
    cat(sprintf("  lethality index had maximum value %.4g, at time %.4g s\n", LI[peakLI], object$t[peakLI]))
    timeOfDanger <- diff(range(object$t[LI > 0.5]))
    cat(sprintf("  lethality index exceeded 0.5 for %.4g s\n", timeOfDanger))
}
