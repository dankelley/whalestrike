#' Ship displacement in kg based on vessel type and length
#'
#' This is done using formulae in Table 3 of Mayette and
#' Brillant (2026).
#'
#' The formulae used are as follows.
#'
#' | `type`                | Formula         |
#' | :---                  | ---:            |
#' | "Bulk Carrier"        |   5.64 * L^3.06 |
#' | "Container Ship"      |  86.40 * L^2.46 |
#' | "Cruise"              |  97.51 * L^2.28 |
#' | "Ferry"               |  25.15 * L^2.62 |
#' | "Fishing"             |   0.71 * L^3.79 |
#' | "Government/Research" |   2.95 * L^3.22 |
#' | "Other"               |   2.64 * L^3.35 |
#' | "Passenger"           |   4.32 * L^3.08 |
#' | "Pleasure Craft"      |  34.47 * L^2.68 |
#' | "Sailing"             |   1.23 * L^3.53 |
#' | "Tanker"              |   7.25 * L^3.03 |
#' | "Tug"                 | 104.48 * L^2.51 |
#'
#' @param type either (1) a string identifying the ship type, in
#' which case the average overall length of the named vessel is
#' returned, or (2) NULL, in which case a vector of permitted
#' values of `type` is returned.
#'
#' @param L vessel length in metres.
#'
#' @return `shipMassFromLength` returns ship displacement mass
#' (in kg), according to Mayette and Brillant (2026) Table 3.
#'
#' @examples
#' library(whalestrike)
#' shipMassFromLength("Tug", 50) / 1e3 # 1920.648
#'
#' @references
#'
#' Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method
#' to Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models."
#' PloS One 21, no. 1 (2026): e0339760.
#' https://doi.org/10.1371/journal.pone.0339760.
#'
#' @export
#'
#' @author Dan Kelley, with help from Alexandra Mayette
#' @family functions relating to ship characteristics
shipMassFromLength <- function(type = NULL, L) {
    knownTypes <- c(
        paste("Bulk", "Carrier"),
        paste("Container", "Ship"),
        "Cruise",
        "Ferry",
        "Fishing",
        "Government/Research",
        "Other",
        "Passenger",
        paste("Pleasure", "Craft"),
        "Sailing",
        "Tanker",
        "Tug"
    )
    if (is.null(type)) {
        return(knownTypes)
    }
    if (!(type %in% knownTypes)) {
        stop(
            "type=\"", type, "\" not handled; try one of \"",
            paste(knownTypes, collapse = "\", \""), "\""
        )
    }
    switch(type,
        "Bulk Carrier" = 5.64 * L^3.06, # ∆ = 5.64 ∗ LOA3.06
        "Container Ship" = 86.40 * L^2.46, # ∆ = 86.40 ∗ LOA2.46
        "Cruise" = 97.51 * L^2.28, # ∆ = 97.51 ∗ LOA2.28
        "Ferry" = 25.15 * L^2.62, # ∆ = 25.15 ∗ LOA2.62
        "Fishing" = 0.71 * L^3.79, # ∆ = 0.71∗LOA3.79
        "Government/Research" = 2.95 * L^3.22, # ∆ = 2.95∗LOA3.22
        "Other" = 2.64 * L^3.35, # ∆ = 2.64 ∗ LOA3.35
        "Passenger" = 4.32 * L^3.08, # ∆ = 4.32 ∗ LOA3.08
        "Pleasure Craft" = 34.47 * L^2.68, # ∆ = 34.47∗LOA2.68
        "Sailing" = 1.23 * L^3.53, # ∆ = 1.23∗LOA3.53
        "Tanker" = 7.25 * L^3.03, # ∆ = 7.25 ∗ LOA3.03
        "Tug" = 104.48 * L^2.51 # ∆ = 104.48 ∗ LOA2.51
    )
}


# Next is from copy/paste of Mayette and Brillant (2026) Table 1. I added
# the column headings. I removed the DWT columns because I think the unit
# is wrongly stated in the publication, and I only want this table for the
# 'avg' (meaning average LOA in m).
# type                  n1 n2   n   low  high   avg   SD   DWT DWTsd
# Bulk Carrier         100  0 100  85.9 300.0 221.0 50.1 84444 61877
# Container Ship       100  0 100  84.0 336.6 240.7 61.1 50843 29726
# Cruise                92  0  92  90.6 344.3 219.1 69.7  6059  3961
# Ferry                 47  0  47  25.5 203.3 100.1 45.0  1586  2022
# Fishing               39 61 100   4.6 104.5  25.9 22.2   794   806
# Government/ Research  48  4  52   7.5 182.5  63.7 31.1  1617  3519
# Other                 43  6  49   6.2 178.8  83.3 45.3  5284  5766
# Passenger              4 26  30   6.2  95.1  23.2 20.9   609   733
# Pleasure Craft        42 20  62   4.9  81.2  42.2 23.0   286   389
# Sailing                7 32  39   4.6  76.0  20.3 18.2   512   440
# Tanker               100  0 100 110.0 277.0 198.4 44.2 63904 43735
# Tug                  102  0 102  25.2  95.0  44.1 18.6   952  1204

#' Nominal ship length in m
#'
#' This is based on the "Average LOA in m" column in Table 1
#' of Mayette and Brillant (2026).
#'
#' @param type either (1) a string identifying the ship type, in
#' which case the average overall length of the named vessel is
#' returned, or (2) NULL, in which case a data frame containing
#' type and length is returned.
#'
#' @return
#' `shipLength` returns ship length in m, as defined in Mayette
#' and Brillant (2026).
#'
#' @references
#'
#' Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method
#' to Estimate Vessel Mass for Use in Whale-Ship Strike Risk Models."
#' PloS One 21, no. 1 (2026): e0339760.
#' https://doi.org/10.1371/journal.pone.0339760.
#'
#' @examples
#' library(whalestrike)
#' # An individual length
#' shipLength("Fishing")
#' # A table of lengths
#' shipLength()
#'
#' @export
#'
#' @author Dan Kelley, with help from Alexandra Mayette
#' @family functions relating to ship characteristics
shipLength <- function(type = NULL) {
    knownTypes <- c(
        paste("Bulk", "Carrier"),
        paste("Container", "Ship"),
        "Cruise",
        "Ferry",
        "Fishing",
        "Government/Research",
        "Other",
        "Passenger",
        paste("Pleasure", "Craft"),
        "Sailing",
        "Tanker",
        "Tug"
    )
    if (is.null(type)) {
        return(data.frame(type = knownTypes, length = unlist(lapply(knownTypes, shipLength))))
    }
    if (!(type %in% knownTypes)) {
        stop(
            "type=\"", type, "\" not handled; try one of \"",
            paste(knownTypes, collapse = "\", \""), "\""
        )
    }
    switch(type,
        "Bulk Carrier" = 221.0,
        "Container Ship" = 240.7,
        "Cruise" = 219.1,
        "Ferry" = 100.1,
        "Fishing" = 25.9,
        "Government/Research" = 63.7,
        "Other" = 83.3,
        "Passenger" = 23.2,
        "Pleasure Craft" = 42.2,
        "Sailing" = 20.3,
        "Tanker" = 198.4,
        "Tug" = 44.1
    )
}

#' Compute ship wetted area from mass
#'
#' Estimate the wetted area of a Cape Islander boat,
#' given the vessel mass.
#'
#' The method is based on scaling up the results for a single Cape
#' Islander ship, of displacement 20.46 tonnes, length 11.73m,
#' beam 4.63m, and draft 1.58m, on the assumption that the wetted area
#' is proportional to
#' \eqn{length*(2*draft+beam)}{length*(2*draft+beam)}.
#' This reference area is scaled to
#' the specified mass, `ms`, by multiplying by the 2/3
#' power of the mass ratio.
#'
#' Note that this is a crude calculation meant as a stop-gap measure, for
#' estimates values of the `Ss` argument to [parameters()].
#' It should not be used in preference to inferences
#' made from architectural drawings of a given ship under study.
#'
#' @param ms Ship mass (kg).
#'
#' @return Estimated area (m^2).
#'
#' @author Dan Kelley
#'
#' @export
#' @family functions relating to ship characteristics
shipAreaFromMass <- function(ms) {
    length <- 11.73 # m
    beam <- 4.63 # m
    draft <- 1.58 # m
    displacement <- 20.46e3 # m^3
    factor <- (ms / displacement)^(1 / 3) # lengthscale factor
    length * (beam + 2 * draft) * factor^2
}

