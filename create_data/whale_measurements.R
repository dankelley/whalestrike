# Data from Alexandra Mayette, in a MSoft document sent 2026-01-28 to Dan Kelley
# Names: I changed some names to match the package (e.g. gray/grey).
# Values: I divided bone thickness by 2 (the variable used in the package)
# Units: length is in metres; all others are in centimetres.
whale_measurements <- read.csv("whale_measurements.csv", comment="#")
whale_measurements$bone <- 0.5 * whale_measurements$bone # FIX_mE: check with A_m
save(whale_measurements, file="whale_measurements.rda")
