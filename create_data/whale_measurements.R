# Data from Alexandra Mayette, in a MSoft document sent 2026-01-28 to Dan Kelley
# See the header in whale_measurements.csv and also the documentation
# provided by typing '?whale_measurements' for more details.
whale_measurements <- read.csv("whale_measurements.csv", comment="#")
whale_measurements$girth <- NULL # not used in package
# Convert to metres, which is the unit used by 'parameters()'.
whale_measurements$skin <- 0.01 * whale_measurements$skin
whale_measurements$blubber <- 0.01 * whale_measurements$blubber
whale_measurements$sublayer <- 0.01 * whale_measurements$sublayer
whale_measurements$bone <- 0.01 * whale_measurements$bone
save(whale_measurements, file="whale_measurements.rda")
