## strike_accel_auto.R

library(whalestrike)

if (!interactive()) pdf("strike_accel_auto.pdf", pointsize=9)
par(mar=c(3,3,3,1), mgp=c(2,0.7,0), mfcol=c(3,4))

a <- '"l1","l2","l3","l4","lw","Ly","Lz","ms","theta","tmax","vs"
0.025,0.16,1.08,0.1,13.2,0.5,1,311000,50,1,8.231104'
d <- read.csv(text=a)
parm <- parameters(l=c(d$l1, d$l2, d$l3, d$l4), lw=d$lw, Ly=d$Ly, Lz=d$Lz, ms=d$ms, theta=d$theta)


knot2ms <- 0.5144 # m/s per knot
for (knot in c(10, 15, 20, 25)) {
     state <- list(xs=-1.4, vs=knot*knot2ms, xw=0, vw=0)
     t <- seq(0, d$tmax, length.out=500)
     userTime <- system.time(sol <- strike(t, state, parm))[[1]]
     max <- which.max(sol$dvwdt)
     plot(t, sol$dvwdt/9.8, xlim=t[max+c(-10,10)],
          type='o', pch=20,
          xlab="Time [s]", ylab="Whale accel. [g]")
     accelmax <- max(abs(sol$dvwdt))
     accelmed <- median(abs(sol$dvwdt))
     mtext(sprintf("vs=%.1fknot: dt=%.4fs, userTime=%.1fs\nmaxAccel=%.1f g, medAccel=%.1f g",
                   state$vs/knot2ms, t[2]-t[1], userTime, accelmax, accelmed),
           side=3, line=0, cex=0.7, col=2)
     if (accelmax > 100 * accelmed) {
         ## sqrt(A*E/(l*m)) is the bounce timescale
         NEED <- 5 # hoped-for number of points in peak
         dt <- (1/NEED) * 0.5 * sqrt(parm$l[4] * parm$mw / (parm$Ly*parm$Lz*parm$a[4]*parm$b[4]))
         t <- seq(t[1], tail(t, 1), dt)
         warning("redoing with ", length(t), " report times\n")
         userTime <- system.time(sol <- strike(t, state, parm))[[1]]
         max <- which.max(sol$dvwdt)
         plot(t, sol$dvwdt/9.8, xlim=t[max+c(-10,10)],
              type='o', pch=20,
              xlab="Time [s]", ylab="Whale accel. [g]")
         accelmax <- max(abs(sol$dvwdt))
         accelmed <- median(abs(sol$dvwdt))
         mtext(sprintf("vs=%.1fknot: dt=%.4fs, userTime=%.1fs\nmaxAccel=%.1f g, medAccel=%.1f g",
                       state$vs/knot2ms, t[2]-t[1], userTime, accelmax, accelmed),
               side=3, line=0, cex=0.7, col=2)

         NEED <- 10 # hoped-for number of points in peak
         dt <- (1/NEED) * 0.5 * sqrt(parm$l[4] * parm$mw / (parm$Ly*parm$Lz*parm$a[4]*parm$b[4]))
         t <- seq(t[1], tail(t, 1), dt)
         warning("redoing with ", length(t), " report times\n")
         userTime <- system.time(sol <- strike(t, state, parm))[[1]]
         max <- which.max(sol$dvwdt)
         plot(t, sol$dvwdt/9.8, xlim=t[max+c(-10,10)],
              type='o', pch=20,
              xlab="Time [s]", ylab="Whale accel. [g]")
         accelmax <- max(abs(sol$dvwdt))
         accelmed <- median(abs(sol$dvwdt))
         mtext(sprintf("vs=%.1fknot: dt=%.4fs, userTime=%.1fs\nmaxAccel=%.1f g, medAccel=%.1f g",
                       state$vs/knot2ms, t[2]-t[1], userTime, accelmax, accelmed),
               side=3, line=0, cex=0.7, col=2)





     } else {
         plot(0:1, 0:1, axes=FALSE, xlab="", ylab="", type="n")
         text(0.5, 0.5, "Above is OK since\nmax_accel/med_accel < 100")
         box()
         plot(0:1, 0:1, axes=FALSE, xlab="", ylab="", type="n")
         text(0.5, 0.5, "Above is OK since\nmax_accel/med_accel < 100")
         box()
     }
}
if (!interactive()) dev.off()

