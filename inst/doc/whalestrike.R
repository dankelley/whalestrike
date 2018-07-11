## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=500)
state <- c(xs=-2.5, vs=10*0.5144, xw=0, vw=0) # 10 knot ship
parms <- parameters(ms=20e3,
                    impactWidth=0.5, impactHeight=1,
                    lw=10,
                    alpha=0.025, Ealpha=2e7, theta=45,
                    beta=0.25, Ebeta=6e5,
                    gamma=0.5, Egamma=4e5)
sol <- strike(t, state, parms)
par(mfcol=c(2, 2), mar=c(2, 3, 0.5, 0.5), mgp=c(2, 0.7, 0), cex=0.7)
plot(sol)

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=500)
state <- c(xs=-1.5, vs=5, xw=0, vw=0)
beta <- seq(0.1, 0.3, length.out=100)
maxAccel <- rep(NA, length(beta))
ms <- 20e3
for (i in seq_along(beta)) {
    parms <- parameters(ms=ms, impactWidth=0.5, impactHeight=1,
                        lw=10,
                        alpha=0.025, Ealpha=2e7, theta=45,
                        beta=beta[i], Ebeta=6e5)
    sol <- strike(t, state, parms)
    maxAccel[i] <- max(abs(diff(sol$vw))) / (t[2] - t[1])
}
plot(beta, maxAccel, type="l", xlab=expression("Blubber thickness [m]"), ylab="Max. Acceleration [m/s^2]")

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=500)
## Hint: making x and y of different lengths, to avoid row,col
## versus i,j confusion.
beta <- seq(0.1, 0.3, length.out=9)
speedK <- seq(2, 15, length.out=10) # in knots
speed <- 0.5144 * speedK
maxStrain <- matrix(NA, nrow=length(speed), ncol=length(beta))
for (i in seq_along(beta)) {
    for (j in seq_along(speed)) {
        state <- c(xs=-1.5, vs=speed[j], xw=0, vw=0)
        parms <- parameters(ms=ms, impactWidth=0.5, impactHeight=1,
                            lw=10,
                            alpha=0.025, Ealpha=2e7, theta=45,
                            beta=beta[i], Ebeta=6e5)
        sol <- strike(t, state, parms)
        maxStrain[j, i] <- max(sol$WCF$strain)
    }
}
danger <- 0.8 / 1.2
contour(speedK, beta, maxStrain, level=seq(0, danger, 0.1),
        xlab="Speed [knots]", ylab="Blubber thickness [m]")
contour(speedK, beta, maxStrain, level=seq(1, danger, -0.1), lwd=3, add=TRUE)
contour(speedK, beta, maxStrain, level=2/3, lwd=3, lty=2, add=TRUE)
## 2/3 = 0.8/1.2 is Grear et al. 2038 result (page 144, left column) that
## mean (seal) blubber strength is 0.8MPa, whereas mean (seal) blubber
## modulus is 1.2MPa.
mtext("Contours of max strain (dangerous if > 2/3)", side=3)

