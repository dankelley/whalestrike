## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 0.5, length.out=500)
state <- c(xs=-1.5, vs=5, xw=0, vw=0)
parms <- parameters(ms=20e3, Ss=15*pi*3, B=3, D=1.5,
              lw=10, Sw=10*2*pi*3)
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
    parms <- parameters(ms=20e3, Ss=15*pi*3, B=3, D=1.5,
                        lw=10, Sw=10*2*pi*3, delta=0.02, Es=2e7, theta=45,
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
beta <- seq(0.1, 0.4, length.out=19)
speedK <- seq(2, 15, length.out=20) # in knots
speed <- 0.5144 * speedK
maxAccel <- matrix(NA, nrow=length(speed), ncol=length(beta))
for (i in seq_along(beta)) {
    for (j in seq_along(speed)) {
        cat(i, ' ', j, '\n')
        state <- c(xs=-1.5, vs=speed[j], xw=0, vw=0)
        parms <- parameters(ms=20e3, Ss=15*pi*3, B=3, D=1.5,
                            lw=10, Sw=10*2*pi*3, delta=0.02, Es=2e7, theta=45,
                            Eb=1e6, beta=beta[i])
        sol <- strike(t, state, parms)
        maxAccel[j, i] <- max(abs(diff(sol$vw))) / (t[2] - t[1])
    }
}
contour(speedK, beta, maxAccel, xlab="Speed [knots]", ylab="Blubber thickness [m]")
mtext("Contours of acceleration [m/s^2]", side=3)

