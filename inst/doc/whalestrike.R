## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 0.5, length.out=500)
state <- c(xs=-1.5, vs=5, xw=0, vw=0)

ms <- 20e3                             # ship mass
Ss <- 11.73*(2*1.58+4.63)              # ship wetted area (max)
impactWidth <- 4.6                     # impact width
impactHeight <- 1.6                    # impact height
lw <- 10                               # whale length
mw <- whaleMassFromLength(lw)          # whale mass
Sw <- whaleAreaFromLength(lw,"wetted") # whale area
parms <- parameters(ms=ms, Ss=Ss, impactWidth=impactWidth, impactHeight=impactHeight,
                    lw=lw, mw=mw, Sw=Sw)
sol <- strike(t, state, parms)
par(mfcol=c(3, 3), mar=c(2, 3, 0.5, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
plot(sol)

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=500)
state <- c(xs=-1.5, vs=5, xw=0, vw=0)
beta <- seq(0.1, 0.3, length.out=100)
maxAccel <- rep(NA, length(beta))
for (i in seq_along(beta)) {
    parms <- parameters(ms=ms, Ss=Ss, impactWidth=impactWidth, impactHeight=impactHeight,
                        lw=lw, mw=mw, Sw=Sw,
                        delta=0.02, Es=2e7, theta=45,
                        Eb=1e6, beta=beta[i])
    sol <- strike(t, state, parms)
    maxAccel[i] <- max(abs(diff(sol$vw))) / (t[2] - t[1])
}
plot(beta, maxAccel, type="l", xlab=expression("Blubber thickness [m]"), ylab="Max. Acceleration [m/s^2]")

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=500)
## Hint: making x and y of different lengths, to avoid row,col
## versus i,j confusion.
beta <- seq(0.1, 0.4, length.out=9)
speedK <- seq(2, 15, length.out=10) # in knots
speed <- 0.5144 * speedK
maxStrain <- matrix(NA, nrow=length(speed), ncol=length(beta))
for (i in seq_along(beta)) {
    for (j in seq_along(speed)) {
        state <- c(xs=-1.5, vs=speed[j], xw=0, vw=0)
        parms <- parameters(ms=ms, Ss=Ss, impactWidth=impactWidth, impactHeight=impactHeight,
                            lw=lw, mw=mw, Sw=Sw,
                            delta=0.02, Es=2e7, theta=45,
                            Eb=1e6, beta=beta[i])
        sol <- strike(t, state, parms)
        maxStrain[j, i] <- max(sol$WCF$strain)
    }
}
contour(speedK, beta, maxStrain, level=seq(0, 0.3, 0.05),
        xlab="Speed [knots]", ylab="Blubber thickness [m]")
contour(speedK, beta, maxStrain, level=seq(0.3, 1, 0.05), lwd=3, add=TRUE)
mtext("Contours of max strain (bold if unsafe)", side=3)

