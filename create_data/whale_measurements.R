# Data from Alexandra Mayette, in a MSoft document sent 2026-01-28 to Dan Kelley
# See the header in whale_measurements.csv and also the documentation
# provided by typing '?whale_measurements' for more details.
whaleMeasurements <- read.csv("whale_measurements.csv", comment="#")
whaleMeasurements$girth <- NULL # not used in package
# Convert to metres, which is the unit used by 'parameters()'.
whaleMeasurements$skin <- 0.01 * whaleMeasurements$skin
whaleMeasurements$blubber <- 0.01 * whaleMeasurements$blubber
whaleMeasurements$sublayer <- 0.01 * whaleMeasurements$sublayer
whaleMeasurements$bone <- 0.01 * whaleMeasurements$bone
#save(whale_measurements, file="whale_measurements.rda")
dput(whaleMeasurements)
