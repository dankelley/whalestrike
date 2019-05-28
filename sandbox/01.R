## look at a 'cut whale' scenario ... why does bone not thin more?

library(whalestrike)
P<-function(x,y,...){plot(x,y,ylim=c(0,max(y)),...);abline(h=0,lty=3,col=2)}

t <- seq(0, 0.2, length.out=200)
state <- c(xs=-2, vs=25*0.5144, xw=0, vw=0) # 10 knot ship
parms <- parameters(ms=20e3, lw=13)
sol <- strike(t, state, parms)
par(mfcol=c(2, 2), mar=c(3, 3, 0.5, 2), mgp=c(2, 0.7, 0), cex=0.7)
plot(sol$t, sol$WCF$force/1e6, ylab="Force [MN]", type="l")
plot(sol$t, sol$WCF$compressed[,2], type="l");abline(h=0, lty=3, col=2)
plot(sol$t, sol$WCF$compressed[,3], type="l");abline(h=0, lty=3, col=2)
plot(sol$t, sol$WCF$compressed[,4], type="l");abline(h=0, lty=3, col=2)

