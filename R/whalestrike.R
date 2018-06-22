library(deSolve)

#' Whale mass inferred from length
#'
#' Calculate an estimate of the mass of a whale, based on its length.
#'
#' The permitted values for \code{model} are as follows.
#'\itemize{
#' \item \code{"moore2005"} yields
#' \eqn{242.988 * exp(0.4 * length)}{242.988 * exp(0.4 * L)},
#' which (apart from a unit change on \code{L}) is the regression equation
#' shown above Figure 1d in Moore et al. (2005) for right whales.
#'
#' \item \code{"fortune2012atlantic"} yields the formula
#' \eqn{exp(-10.095 + 2.825*log(100*L))}{exp(-10.095 + 2.825*log(100*L))}
#' for North Atlantic right whales, according to a guess made about
#' an erroroneous formula given in the caption
#' of Figure 4 in Fortune et al (2012).
#' (The assumed error is an exchange of slope and intercept.)
#'
#' \item \code{"fortune2012pacific"} yields the formula
#' for North Atlantic right whales, according to a guess made about
#' \eqn{exp(-12.286 + 3.158*log(100*L))}{exp(-12.286 + 3.158*log(100*L))}
#' an erroroneous formula given in the caption
#' of Figure 4 in Fortune et al (2012).
#' (The assumed error is an exchange of slope and intercept.)
#'}
#'
#' @param L Whale length in m.
#' @param model Character string specifying the model (see \dQuote{Details}).
#' @return Mass in kg.
#' @examples
#' library(whalestrike)
#' L <- seq(5, 15, length.out=100)
#' plot(L, massFromLength(L)/1000, type='l',
#'      xlab="Right-whale Length [m]", ylab="Mass [tonne]")
#' lines(L, massFromLength(L, model="fortune2012atlantic")/1000, col=2)
#' lines(L, massFromLength(L, model="fortune2012pacific")/1000, col=3)
#' legend("topleft", lwd=1, col=1:3,
#'        legend=c("moore2005", "fortune2012atlantic", "fortune2012pacific"))
#' grid()
#' @references
#' Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
#' “Morphometry, Gross Morphology and Available Histopathology in North Atlantic
#' Right Whale (Eubalaena Glacialis) Mortalities (1970 to 2002).” Journal of
#' Cetacean Research and Management 6, no. 3 (2005): 199–214.
#'
#' Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. “Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis).” Journal of Mammalogy 93,
#' no. 5 (2012): 1342–54. https://doi.org/10.1644/11-MAMM-A-297.1.
massFromLength <- function(L, model="moore2005")
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

#' Contact area
#'
#' Area of contact between vessel and whale skin
#'
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
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
#' @param xw Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}). The
#' only entries used are blubber thickness \code{beta} [m] and blubber model \code{blubbermodel};
#' the latter is a string which may be \code{"linear"}, \code{"exponential"} or \code{"Raymond"}).
#'
#' @return Compression-resisting force [N]
#' @references
#' J. J. Raymond. Development of a numerical model to predict impact forces on a
#' North Atlantic Right Whale during collision with a vessel.
#' PhD thesis, University of New Hampshire, 2007.
blubberForce <- function(xs, xw, parms)
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
    ##.     cat(sprintf("blubberForce(): xs %.3f, xw %.3f, B %.3f\n", xs, xw, parms$B))
    ##.     cat(sprintf("blubberForce(): xs %.3f, xw %.3f, touching %d, strain %.5f\n", xs, xw, touching[1], strain[1]))
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
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{strike}}). The
#' quantities used by \code{skinForce} are \code{B} (ship beam width, m),
#' \code{D} (ship draft, in m), \code{Es} (skin elastic modulus in Pa),
#' \code{delta} (skin thickness in m), and \code{theta} (skin bevel angle
#' degrees, measured from a vector normal to undisturbed skin).
#'
#' @return Stretch-resisting skin force [N]
skinForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$beta)
    ##> if (is.na(touching[1])) stop("skinForce(): touching is NA")
    dx <- parms$beta - (xw - xs)
    l <- dx * tan(parms$theta * pi / 180) # NB: theta is in deg
    s <- dx / cos(parms$theta * pi / 180) # NB: theta is in deg
    ## Strains in y and z
    epsilony <- 2 * (s - l) / (parms$B + 2 * l)
    epsilonz <- 2 * (s - l) / (parms$D + 2 * l)
    ## Stresses in y and z
    sigmay <- parms$E * epsilony
    sigmaz <- parms$E * epsilonz
    ## Net force in x
    Fx <- 2 * parms$delta * (parms$B * sigmay + parms$D * sigmaz)
    ifelse(touching, Fx, 0)
}

#' Water drag force
#'
#' Compute water drag using a quadratic law.
#'
#' @section Development note:
#' A fixed plan area of 20m^2 is assumed, although that is subject to change.
#' There are published mass-length relationships, and using a shape
#' parameter, we should be able to express the area in terms of whale mass.
#'
#' @param v Whale velocity [m/s]
#' @param CD Drag coefficient. The default corresponds to that of the Empire State building.
#'
#' @return Water drag force [N]
waterForce <- function(v, CD=1.4)
{
    rho <- 1024
    CD <- 1.4 ## Empire State Building 1.3-1.5 (https://en.wikipedia.org/wiki/Drag_coefficient)
    R <- 2                             # whale radius [m]
    L <- 10                            # whale length [m]
    area <- 2 * R * L                  # frontal area
    -0.5 * rho * CD * area * v^2
}

#' Dynamical law
#'
#' @param t time [s].
#' @param y model state, a vector containing ship position \code{xs} [m],
#' ship speed \code{vs} [m/s], whale position \code{xw} [m],
#' and whale speed \code{vw} [m/s].
#' @param parms model parameters.
dynamics <- function(t, y, parms)
{
    xs <- y[1]                         # ship position
    vs <- y[2]                         # ship velocity
    xw <- y[3]                         # whale position
    vw <- y[4]                         # whale velocity
    Fb <- blubberForce(xs, xw, parms)
    if (is.na(Fb)) stop("Fb is NA")
    Fs <- skinForce(xs, xw, parms)
    if (is.na(Fs)) stop("Fs is NA")
    Fw <- waterForce(vw)
    if (is.na(Fw)) stop("Fw is NA")
    if (is.na(parms$ms)) stop("parms$ms is NA")
    if (is.na(parms$mw)) stop("parms$mw is NA")
    F <- Fb + Fs + Fw
    list(c(dxsdt=vs, dvsdt=-F/parms$ms, dxwdt=vw, dvwdt=F/parms$mw))
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
#' \code{\link{skinForce}},
#' \code{\link{blubberForce}}, and
#' \code{\link{waterForce}}.
#'
#' @param t time [s].
#' @param state A named vector holding the initial state of the model:
#' ship position \code{xs} [m],
#' ship speed \code{vs} [m/s],
#' whale position \code{xw} [m]
#' and whale speed \code{vw} [m/s].
#' @param parms A named list holding model parameters:
#'\describe{
#'\item{ms}{Ship mass [kg]}
#'\item{B}{Ship width beam [m]}
#'\item{D}{Ship draft [m]}
#'\item{mw}{Whale mass [kg]. (Consider using \code{\link{massFromLength}},
#' if length data are easier to obtain than mass data.)}
#'\item{delta}{Whale skin thickness [m]. See e.g. Miller et al. (2011), for
#' right whales.}
#'\item{Es}{Whale skin elastic modulus [Pa]}
#'\item{theta}{Whale skin deformation angle [deg]}
#'\item{beta}{Whale blubber thickness [m]}
#'\item{blubbermodel}{Blubber model to use (passed to \code{\link{blubberForce}})}
#'}
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
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
#' parms <- list(ms=20e3, B=3, D=3,
#'               mw=20e3,
#'               delta=0.02, Es=2e7, theta=30,
#'               beta=0.3, blubbermodel="exponential")
#' sol <- strike(t, state, parms)
#' par(mfcol=c(3, 3), mar=c(2, 3, 1, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(sol, which="all")
#'
#' @references
#'\itemize{
#'\item Miller, Carolyn A., Desray Reeb, Peter B. Best, Amy R. Knowlton, Moira W. Brown, and Michael J. Moore. “Blubber Thickness in Right Whales Eubalaena Glacialis and Eubalaena Australis Related with Reproduction, Life History Status and Prey Abundance.” Marine Ecology Progress Series 438 (2011): 267–83.
#'}
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
    ##. print("in strike(), state is:")
    ##. print(state)
    ##. print("in strike(), parms is:")
    ##. print(parms)
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
#' @details
#' Some panels may have additional horizontal lines for thresholds:
#' \itemize{
#' \item for \code{which="whale acceleration"}, lines will be drawn at 8g and 15g,
#' which are \strong{provisional} estimates of 50 percent
#' lethality, determined by visual inspection of the lethality-acceleration diagram
#' Figure 9 of Ng et al. 2017, for NFL (National Football League)
#' concussions and combat concussions, respectively.
#' \item for \code{which="blubber thickness"}, a red line is drawn at 2/3 of the
#' initial thickness, based on the Grear et al. 2018 result (page 144, left column)
#' that mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
#' modulus is 1.2MPa, for a ratio of 2/3.
#'}
#'
#' @param x An object inheriting from class \code{strike}
#' @param which A character vector that indicates what to plot.
#' This choices for its entries are:
#' \code{"location"} for a time-series plot of boat location \code{xw} in black,
#' whale location \code{x} in blue, and skin location in dashed blue,
#' \code{"whale acceleration"} for a time-series plot of whale acceleration,
#' \code{"blubber thickness"} for a time-series plot of blubber thickness,
#' \code{"water force"} for a time-series plot of the water-drag force on the whale,
#' \code{"blubber force"} for a time-series plot of the normal force resulting from blubber compression,
#' \code{"blubber stress"} for a time-series plot of the normal stress on the blubber,
#' \code{"skin force"} for a time-series plot of the normal force resulting from skin tension,
#' \code{"skin tension"} for a time-series plot of the skin tension, in the along-skin direction,
#' \code{"values"} for a listing of \code{param} values.
#' or \code{"all"} for all of the above.
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
#'
#' Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt, Adam P. Summers, and Petra Ditsche. “Mechanical Properties of Harbor Seal Skin and Blubber − a Test of Anisotropy.” Zoology 126 (2018): 137–44. \url{https://doi.org/10.1016/j.zool.2017.11.002}
#'
#' Ng, Laurel J., Vladislav Volman, Melissa M. Gibbons, Pi Phohomsiri, Jianxia Cui, Darrell J. Swenson, and James H. Stuhmiller. “A Mechanistic End-to-End Concussion Model That Translates Head Kinematics to Neurologic Injury.” Frontiers in Neurology 8, no. JUN (2017): 1–18. \url{https://doi.org/10.3389/fneur.2017.00269}
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
            grid()
            death <- which(xs >= xw)[1]
            tdeath <- if (is.finite(death)) t[death] else NA
            if (is.finite(tdeath)) {
                abline(v=tdeath, lwd=2, col="blue")
                mtext("Fatality", at=tdeath, side=3, line=0, col="blue", font=2)
            }
            tclosest <- t[which.min(abs(xs-xw))]
            abline(v=tclosest, col="darkgreen", lwd=2, lty=3)
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
            text(par('usr')[2], aCFL+7, "CFL concussion", pos=2)
            aConflict <- 15 * 9.8
            abline(h=aConflict, col="red", lwd=2) # FIXME: determine values more precisely
            text(par('usr')[2], aConflict+7, "conflict concussion", pos=2)
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
    if (all || "water force" %in% which) {
        plot(t, waterForce(vw) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin force" %in% which) {
        Fs <- skinForce(xs, xw, x$parms)
        plot(t, Fs/1e6, type="l", xlab="Time [s]", ylab="Skin Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "blubber force" %in% which) {
        Fb <- blubberForce(xs, xw, x$parms)
        plot(t, Fb/1e6, type="l", xlab="Time [s]", ylab="Blubber Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin tension" %in% which) {
        circumferance <- 2 * x$parms$B + 2 * x$parms$D
        stresss <- skinForce(xs, xw, x$parms) / (x$parms$delta * circumferance)
        plot(t, stresss/1e6, type="l", xlab="Time [s]", ylab="Skin Tension [MPa]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "blubber stress" %in% which) {
        stresss <- skinForce(xs, xw, x$parms) / (x$parms$delta*x$parms$D)
        A <- contactArea(xs, xw, x$parms)
        stressb <- ifelse(A, blubberForce(xs, xw, x$parms) / A, 0)
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


