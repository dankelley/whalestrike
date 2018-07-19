d <- read.csv("raymond2007_fig2_13_imager.dat", skip=7, header=TRUE)
raymond2007 <- list(strain=d$strain, stress=d$stress)
save(raymond2007, file="raymond2007.rda")
tools::resaveRdaFiles("raymond2007.rda", compress="auto")
plot(raymond2007$strain, raymond2007$stress, type="o")

