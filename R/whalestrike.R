library(deSolve)

#' Contact area
#'
#' Area of contact between vessel and whale skin
#'
#' @param x Boat position [m]
#' @param X Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Contact area [m^2].
contactArea <- function(x, X, parms)
{
    touching <- X < x & x < (X + parms$B)
    ifelse(touching, parms$beam * parms$draft, 0)
}

#' Blubber force
#'
#' @param x Boat position [m]
#' @param X Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Compression-resisting force [N]
blubberForce <- function(x, X, parms)
{
    touching <- X < x & x < (X + parms$B)
    strain <- 1 - (x - X) / parms$B
    ## nodd <- sum(strain<0) && touching
    ## if (is.na(nodd)) {message('BAD'); browser()}
    ## if (nodd > 0)
    ##     message('x[1]=', x[1], ', X[1]=', X[1], ' # negative strains=', nodd)
    ifelse(touching,
           if (parms$blubbermodel=="Raymond") 636273*strain
           else if (parms$blubbermodel=="linear") 580004.6*strain
           else if (parms$blubbermodel=="exponential") -1.712240e+05*(1-exp(strain/4.185693e-01))
           else NA,
           0) * contactArea(x, X, parms)
}

#' Skin force
#'
#' @param x Boat position [m]
#' @param X Whale position [m]
#' @param parms Parameters of the simulation (as provided to \code{\link{simulate}}).
#'
#' @return Stretch-resisting skin force [N]
skinForce <- function(x, X, parms)
{
    touching <- X < x & x < (X + parms$B)
    dx <- parms$B - (x - X)
    sarea <- parms$delta * parms$draft # skin cross-section area
    ifelse(touching,
           sarea * parms$ES*(sqrt(dx^2+parms$gamma^2)-parms$gamma) / (parms$beam/2+parms$gamma),
           0)
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
#' @param y model state, a vector containing boat position \code{X} [m],
#' boat speed \code{V} [m/s], whale position \code{x} [m],
#' and whale speed \code{v} [m/s].
#' @param parms model parameters.
dynamics <- function(t, y, parms)
{
    X <- y[1]                          # boat position
    V <- y[2]                          # boat velocity
    x <- y[3]                          # whale position
    v <- y[4]                          # whale velocity
    dxdt <- v
    dXdt <- V
    Fb <- blubberForce(x, X, parms)
    Fs <- skinForce(x, X, parms)
    Fw <- waterForce(v)
    ##Fw <- rep(0, length.out=length(v))
    F <- Fb + Fs + Fw
    list(c(dXdt=V, dVdt=-F/parms$shipMass, dxdt=v, dvdt=F/parms$whaleMass))
}

#' Calculate derivative using first difference
#' @param x variable.
#' @param t time in seconds.
#' @return Derivative estimated by using \code{\link{diff}} on both \code{x}
#' and \code{t}. In order to return a value of the same length as \code{x} and
#' \code{t}, the last value is repeated.
derivative <- function(x, t)
{
    res <- diff(x) / diff(t)
    c(res, tail(res, 1))
}

#' Perform dynamical simulation
#' @param t time [s].
#' @param state model state, a vector containing boat position \code{X} [m],
#' boat speed \code{V} [m/s], whale position \code{x} [m],
#' and whale speed \code{v} [m/s].
#' @param parms a named vector holding model parameters.
#'
#' @return An object of class \code{"whalestrike"}, consisting of a
#' list  containing vectors for time (\code{t} [s]), boat position (\code{X} [m]),
#' boat speed (\code{V} [m/s]), whale position (\code{x} [m]), whale speed (\code{v} [m/s]),
#' boat acceleration (\code{dVdt} [m/s^2]), and whale acceleration (\code{dvdt} [m/s^2]),
#' along with a list containing the model parameters (\code{parms}).
#'
#' @examples
#' t <- seq(0, 1, length.out=500)
#' state <- c(X=-1.5, V=5, x=0, v=0)
#' parms <- list(shipMass=20e3, beam=3, draft=3, V=4.6296,
#'               whaleMass=20e3, B=0.35, gamma=1, delta=0.02, ES=2e+07, blubbermodel="exponential")
#' sol <- whalestrike(t, state, parms)
#' plot(sol)
whalestrike <- function(t, state, parms)
{
    print("in simulate(), state is:")
    print(state)
    print("in simulate(), parms is:")
    print(parms)
    sol <- lsoda(state, t, dynamics, parms)
    print("in simulate(), head(sol) is:")
    print(head(sol))
    cat("head(time)=", paste(sol[1:6,"time"], collapse=" "), "\n")
    cat("head(X)=", paste(sol[1:6,"X"], collapse=" "), "\n")
    cat("head(V)=", paste(sol[1:6,"V"], collapse=" "), "\n")
    cat("head(x)=", paste(sol[1:6,"x"], collapse=" "), "\n")
    cat("head(v)=", paste(sol[1:6,"v"], collapse=" "), "\n")
    cat("state:\n")
    print(state)
    ## Add extra things for plotting convenience. Perhaps should also
    ## calculate the forces here.
    res <- list(t=sol[, "time"],
                X=sol[, "X"],
                V=sol[, "V"],
                x=sol[, "x"],
                v=sol[, "v"],
                dVdt=derivative(sol[, "V"], sol[, "time"]),
                dvdt=derivative(sol[, "v"], sol[, "time"]),
                parms=parms)
    cat("simulate() finished\n")
    class(res) <- "whalestrike"
    res
}

#' Plot a whalestrike object
#'
#' @param x An object inheriting from class \code{whalestrike}
#' @param ... Other arguments (ignored).
## .method plot whalestrike
## .S3method plot whalestrike
plot.whalestrike <- function(x, ...)
{
    showLegend <- FALSE
    g <- 9.8 # gravity
    sol <- x # avoid name clash
    t <- sol$t
    X <- sol$X
    V <- sol$V
    x <- sol$x
    v <- sol$v
    dVdt <- sol$dVdt
    dvdt <- sol$dvdt
    centre <- FALSE
    if (centre) {
        ## Trim so event is centred (maybe; not a big deal)
        end <- 3 * which.min(abs(x-X))
        look <- 1:end
        t <- t[look]
        X <- X[look]
        V <- V[look]
        x <- x[look]
        v <- v[look]
        A <- A[look]
        a <- a[look]
    }
    ## Plots
    showEvents <- function(x, X) {
        grid()
        death <- which(X >= x)[1]
        tdeath <- if (is.finite(death)) t[death] else NA
        if (is.finite(tdeath)) {
            abline(v=tdeath, lwd=2, col="blue")
            mtext("Fatality", at=tdeath, side=3, line=0, col="blue", font=2)
        }
        tclosest <- t[which.min(abs(x-X))]
        abline(v=tclosest, col="darkgreen", lwd=2, lty=3)
    }
    par(mfcol=c(3, 3), mar=c(2, 3, 1, 0.5), mgp=c(2, 0.7, 0), cex=1)
    ## Top-left: position
    ylim <- range(c(X, x), na.rm=TRUE)
    plot(t, X, type="l", xlab="Time [s]", ylab="Location [m]", ylim=ylim, lwd=2)
    lines(t, x, col=2, lwd=2)
    lines(t, x - sol$parms$B, col=2, lty=3, lwd=2)
    showEvents(x, X)
    if (showLegend)
        legend("topleft", col=c(1, 2, 2), legend=c("Ship", "Whale", "Blubber"),
               lwd=rep(2, 3), lty=c(1, 1, 3), bg="white", cex=0.8)

    if (FALSE) {
        plot(t, x-X, type="l", ylab="Separation [m]")
        showEvents(x, X)
        abline(h=0, col="blue", lwd=2, lty=3)
        strain <- 1 - (x - X) / sol$parms$B
        plot(t, ifelse(X<x&x<(X+sol$parms$B), strain, 0), type="l", lwd=2, ylab="Blubber Strain")
        showEvents(x, X)
    }

    if (FALSE) {
        plot(t, contactArea(x, X, sol$parms), type="l", lwd=2, ylab="Area [m^2]")
        showEvents(x, X)
    }


    plot(t, derivative(v, t) / g, xlab="Time [s]", ylab="Whale accel [g]", type="l", lwd=2)
    showEvents(x, X)

    plot.new() # skip panel

    plot(t, waterForce(v) / 1e6 , xlab="Time [s]", ylab="Water force [MN]", type="l", lwd=2)
    showEvents(x, X)

    Fs <- skinForce(x, X, sol$parms)
    plot(t, Fs/1e6, type="l", xlab="Time [s]", ylab="Skin Force [MN]", lwd=2)
    showEvents(x, X)

    Fb <- blubberForce(x, X, sol$parms)
    plot(t, Fb/1e6, type="l", xlab="Time [s]", ylab="Blubber Force [MN]", lwd=2)
    showEvents(x, X)

    A <- contactArea(x, X, sol$parms)
    stresss <- ifelse(A, skinForce(x, X, sol$parms) / (sol$parms$delta*sol$parms$draft), 0)
    plot(t, stresss/1e6, type="l", xlab="Time [s]", ylab="Skin Stress [MPa]", lwd=2)
    showEvents(x, X)

    stressb <- ifelse(A, blubberForce(x, X, sol$parms) / A, 0)
    plot(t, stressb/1e6, type="l", xlab="Time [s]", ylab="Blubber Stress [MPa]", lwd=2)
    showEvents(x, X)
}

