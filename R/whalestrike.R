library(deSolve)

#' Contact area
#'
#' Area of contact between vessel and whale skin
#'
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Contact area [m^2].
contactArea <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$B)
    ifelse(touching, parms$beam * parms$draft, 0)
}

#' Blubber force
#'
#' @param xs Ship position [m]
#' @param xw Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Compression-resisting force [N]
blubberForce <- function(xs, xw, parms)
{
    ##. if (is.na(xs[1])) stop("xs is NA")
    ##. if (is.na(xw[1])) stop("xw is NA")
    ##. if (is.na(parms$B)) stop("parms$B is NA")
    touching <- xs < xw & xw < (xs + parms$B)
    strain <- 1 - (xw - xs) / parms$B
    ##. if (!is.na(touching[1]) && touching[1]) {
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
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Stretch-resisting skin force [N]
skinForce <- function(xs, xw, parms)
{
    touching <- xs < xw & xw < (xs + parms$B)
    ##> if (is.na(touching[1])) stop("skinForce(): touching is NA")
    dx <- parms$B - (xw - xs)
    ##> if (is.na(dx[1])) stop("skinForce(): dx is NA")
    sarea <- parms$delta * parms$draft # skin cross-section area
    ##> if (is.na(sarea)) stop("skinForce(): sarea is NA")
    rval <- ifelse(touching,
                   sarea * parms$Eskin*(sqrt(dx^2+parms$gamma^2)-parms$gamma) / (parms$beam/2+parms$gamma),
                   0)
    ##. if (is.na(rval[1])) browser()
    ##. if (is.na(rval[1])) stop("skinForce(): rval[1] is NA")
    rval
}

#' Water drag force
#'
#' @param v Whale velocity [m/s]
#'
#' @return Water drag force [N]
waterForce <- function(v)
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
#' @param t time [s].
#' @param state A named vector describing the initial state of the model. This
#' must contain ship position \code{xs} [m],
#' ship speed \code{vs} [m/s], whale position \code{xw} [m],
#' and whale speed \code{vw} [m/s].
#' @param parms a named vector holding model parameters.
#' ship mass \code{"ms"} [kg],
#' ship beam \code{"beam"} [m],
#' ship draft \code{"draft"} [m],
#' whale mass \code{"mw"} [kg],
#' whale skin thickness \code{"delta"} [m],
#' whale skin elastic modulus \code{"Eskin"} [Pa],
#' whale skin deformation extension \code{"gamma"} [m],
#' whale blubber thickness \code{"B"} [m],
#' and \code{"blubbermodel"} (one of \code{"linear"}, \code{"exponential"}
#' or \code{"raymond"}).
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
#' parms <- list(ms=20e3, beam=3, draft=3,
#'               mw=20e3,
#'               delta=0.02, Eskin=2e7, gamma=0.2,
#'               B=0.3, blubbermodel="exponential")
#' sol <- whale_strike(t, state, parms)
#' par(mfcol=c(3, 3), mar=c(2, 3, 1, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
#' plot(sol, which="all")
whale_strike <- function(t, state, parms)
{
    if (missing(t)) stop("must supply t")
    if (missing(state)) stop("must supply state")
    if (missing(parms)) stop("must supply parms")
    if (!all(c("xs", "vs", "xw", "vw") %in% names(state)))
        stop("state must contain items named xs, vs, xw and vw")
    if (!all(c("ms", "beam", "draft", "mw", "delta", "Eskin", "gamma", "B", "blubbermodel") %in% names(parms)))
        stop("parms must contain items named ms, beam, draft, mw, delta, Eskin, gamma, B, and blubbermodel")
    ##. print("in simulate(), state is:")
    ##. print(state)
    ##. print("in simulate(), parms is:")
    ##. print(parms)
    sol <- lsoda(state, t, dynamics, parms)
    ##. print("in simulate(), head(sol) is:")
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
    class(res) <- "whale_strike"
    res
}

#' Plot a whale_strike object
#'
#' See \code{\link{whale_strike}} for examples.
#'
#' @param x An object inheriting from class \code{whale_strike}
#' @param which A character vector that indicates what to plot.
#' This choices for its entries are:
#' \code{"location"} for a panel showing time-series of boat location \code{xw} in black,
#' whale location \code{x} in red, and skin location in dashed red,
#' \code{"whale acceleration"} for a time-series plot of whale acceleration,
#' \code{"water force"} for a time-series plot of the water-drag force on the whale,
#' \code{"blubber thickness"} for a time-series plot of blubber thickness,
#' \code{"blubber force"} for a time-series plot of the normal force resulting from blubber compression,
#' \code{"blubber stress"} for a time-series plot of the normal stress on the blubber,
#' \code{"skin force"} for a time-series plot of the normal force resulting from skin tension,
#' \code{"skin tension"} for a time-series plot of the skin tension, in the along-skin direction,
#' or \code{"all"} for all of the above.
#' @param center Logical, indicating whether to center time-series plots
#' on the time when the vessel and whale and in closest proximity.
#' @param indicateEvents Logical, indicating whether to draw lines for some events,
#' such as the moment of closest approach.
#' @param ... Other arguments (ignored).
plot.whale_strike <- function(x, which="all", center=FALSE, indicateEvents=FALSE, ...)
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
    showEvents <- function(xs, xw) {
        if (indicateEvents) {
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
        lines(t, xw, col=2, lwd=2)
        lines(t, xw - x$parms$B, col=2, lty=3, lwd=2)
        showEvents(xs, xw)
        if (showLegend)
            legend("topleft", col=c(1, 2, 2), legend=c("Ship", "Whale", "Blubber"),
                   lwd=rep(2, 3), lty=c(1, 1, 3), bg="white", cex=0.8)
    }
    if (all || "whale acceleration" %in% which) {
        plot(t, derivative(vw, t) / g, xlab="Time [s]", ylab="Whale accel [g]", type="l", lwd=2)
        showEvents(xs, xw)
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
    if (all || "blubber thickness" %in% which) {
        touching <- xs < xw & xw < (xs + x$parms$B)
        thickness <- ifelse(touching, xw-xs, x$parms$B)
        ylim <- c(0, max(thickness))
        plot(t, thickness, xlab="Time [s]", ylab="Blubber thickness", type="l", lwd=2, ylim=ylim)
        showEvents(xs, xw)
    }
    if (all || "blubber force" %in% which) {
        Fb <- blubberForce(xs, xw, x$parms)
        plot(t, Fb/1e6, type="l", xlab="Time [s]", ylab="Blubber Force [MN]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "skin tension" %in% which) {
        stresss <- skinForce(xs, xw, x$parms) / (x$parms$delta*x$parms$draft)
        plot(t, stresss/1e6, type="l", xlab="Time [s]", ylab="Skin Tension [MPa]", lwd=2)
        showEvents(xs, xw)
    }
    if (all || "blubber stress" %in% which) {
        A <- contactArea(xs, xw, x$parms)
        stressb <- ifelse(A, blubberForce(xs, xw, x$parms) / A, 0)
        plot(t, stressb/1e6, type="l", xlab="Time [s]", ylab="Blubber Stress [MPa]", lwd=2)
        showEvents(xs, xw)
    }
}

