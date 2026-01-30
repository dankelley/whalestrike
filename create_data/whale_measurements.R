# Data from Alexandra Mayette, in a MSoft document sent 2026-01-28 to Dan Kelley
# See the header in whale_measurements.csv and also the documentation
# provided by typing '?whale_measurements' for more details.
whale_measurements <- read.csv("whale_measurements.csv", comment="#")
whale_measurements$girth <- NULL # not used in package
save(whale_measurements, file="whale_measurements.rda")
