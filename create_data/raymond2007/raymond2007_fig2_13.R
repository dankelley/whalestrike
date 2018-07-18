d <- read.table("raymond2007_fig2_13.dat", skip=7, header=TRUE)
raymond2007 <- list(strain=d$strain, stress=d$stress)
save(raymond2007, file="raymond2007.rda")
tools::resaveRdaFiles("raymond2007.rda", compress="auto")

