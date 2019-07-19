library(whalestrike)
if (!interactive()) pdf("large_and_small_overall.pdf")
knots <- seq(0, 20, 0.5)
maxStressLarge <- NULL
maxStressSmall <- NULL
maxStressOverStrengthLarge <- NULL
maxStressOverStrengthSmall <- NULL
large <- 400                           # tonnes
small <- 20                            # tonnes
pchLarge <- 1
pchSmall <- 2
parmsLarge <- parameters(ms=large * 1000)
parmsSmall <- parameters(ms=small * 1000)
t <- seq(0, 10, length.out=5000)
for (speed in knot2mps(knots)) {
    state <- list(xs=-2, vs=speed, xw=0, vw=0)
    solLarge <- strike(t, state, parmsLarge)
    maxStressLarge <- c(maxStressLarge, max(solLarge$WCF$stress))
    maxStressOverStrengthLarge <- c(maxStressOverStrengthLarge, max(solLarge$WCF$stress/solLarge$parms$s[2]))
    solSmall <- strike(t, state, parmsSmall)
    maxStressSmall <- c(maxStressSmall, max(solSmall$WCF$stress))
    maxStressOverStrengthSmall <- c(maxStressOverStrengthSmall, max(solSmall$WCF$stress/solSmall$parms$s[2]))
}
par(mfrow=c(2, 1), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
nonzero <- maxStressLarge > 0
plot(knots[nonzero], log10(maxStressLarge[nonzero]), pch=pchLarge, xaxs="i",
     xlab="Ship Speed [knots]", ylab="log10 peak blubber stress")
nonzero <- maxStressSmall > 0
points(knots[nonzero], log10(maxStressSmall[nonzero]), pch=pchSmall)
legend("topleft", pch=c(pchLarge, pchSmall), legend=paste(c(large, small), "ship"))
abline(h=log10(solSmall$parms$s[2]), lty=2)
nonzero <- maxStressLarge > 0
plot(knots[nonzero], log10(maxStressOverStrengthLarge[nonzero]), pch=pchLarge, xaxs="i",
     xlab="Ship Speed [knots]", ylab="log10 peak blubber stress / strength")
nonzero <- maxStressSmall > 0
points(knots[nonzero], log10(maxStressOverStrengthSmall[nonzero]), pch=pchSmall)
abline(h=0, lty=2)
legend("topleft", pch=c(pchLarge, pchSmall), legend=paste(c(large, small), "ship"))
if (!interactive()) dev.off()

