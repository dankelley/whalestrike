## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----results="hide", fig.width=6, fig.height=3---------------------------
library(whalestrike)
t <- seq(0, 1, length.out=200)
state <- c(xs=-2.5, vs=10*0.5144, xw=0, vw=0) # 10 knot ship
parms <- parameters(ms=20e3, Ly=0.5, Lz=1, lw=13)
sol <- strike(t, state, parms)
par(mfcol=c(1, 3), mar=c(2, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
plot(sol)

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=200)
state <- c(xs=-1.5, vs=10*0.5144, xw=0, vw=0) # 10 knots
area <- seq(0.25, 2.5, length.out=50)
stress <- rep(NA, length.out=length(area)) # compressive stress [MPa]
for (i in seq_along(area)) {
    L <- sqrt(area[i])
    parms <- parameters(Ly=L, Lz=L)
    sol <- strike(t, state, parms)
    stress[i] <- max(sol$WCF$stress) / 1e6
}
danger <- parms$s[2] / 1e6
plot(area, stress, type="l", xlab="Area [m^2]", ylab="Stress [MPa]",
     ylim=c(0, max(stress)))
lines(area[stress>=danger], stress[stress>=danger], lwd=3)
abline(h=danger, lty="dashed")
mtext(sprintf("Compression stress [MPa]\n(injurious if > %.2f MPa)", danger),
      side=3, line=1)

## ----results="hide"------------------------------------------------------
library(whalestrike)
t <- seq(0, 1, length.out=200)
## Hint: the following creates x and y of different lengths,
## so that mismatches between row/col and i/j values will
## yield errors.
blubberThickness <- seq(0.1, 0.25, length.out=9)
speedK <- seq(4, 12, length.out=10) # in knots
speed <- 0.5144 * speedK
## stress = peak stress during each simulation, in MPa
stress <- matrix(NA, nrow=length(speed), ncol=length(blubberThickness))
for (i in seq_along(blubberThickness)) {
    for (j in seq_along(speed)) {
        state <- c(xs=-1.5, vs=speed[j], xw=0, vw=0)
        parms <- parameters(l=c(0.025, blubberThickness[i], 0.5, 0.200))
        sol <- strike(t, state, parms)
        stress[j, i] <- max(sol$WCF$stress) / 1e6
    }
}
danger <- parms$s[2] / 1e6
contour(speedK, beta, stress, levels=seq(0, danger, 0.1),
        xlab="Speed [knots]", ylab="Blubber thickness [m]")
contour(speedK, beta, stress, level=danger, lty=2, add=TRUE)
contour(speedK, beta, stress, level=seq(3, danger, -0.1), lwd=2, add=TRUE)
mtext(sprintf("Compression stress [MPa]\n(injurious if > %.2f MPa)", danger),
      side=3, line=1)

