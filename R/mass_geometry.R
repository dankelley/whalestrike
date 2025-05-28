# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

#' Whale mass inferred from length
#'
#' Calculate an estimate of the mass of different species of whale,
#' based on animal length, based on formulae as listed in
#' \dQuote{Details}.
#'
#' The permitted values for `model` and `species` are as follows. Note
#' that if `model` is not provided, then Fortune (2012) is used
#' for `"N. Atl. Right Whale"`, and Lockyer (1976) is used for all
#' other species.
#'
#' * `"moore2005"` (which only works if `species` is `"N. Atl. Right Whale"`) yields
#' \eqn{242.988 * exp(0.4 * length)}{242.988 * exp(0.4 * L)},
#' which (apart from a unit change on `L`) is the regression equation
#' shown above Figure 1d in Moore et al. (2005) for right whales. A
#' difficult in the Moore et al. (2005) use of a single nonzero digit
#' in the multiplier on `L` is illustrated in \dQuote{Examples}.
#'
#' * `"fortune2012"` with `species="N. Atl. Right Whale"` yields the formula
#' \eqn{exp(-10.095 + 2.825*log(100*L))}{exp(-10.095 + 2.825*log(100*L))}
#' for North Atlantic right whales, according to a corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'
#' * `"fortune2012"` with `species="N. Pac. Right Whale"` yields the formula
#' \eqn{exp(-12.286 + 3.158*log(100*L))}{exp(-12.286 + 3.158*log(100*L))}
#' for North Pacific right whales, according to a corrected version of the
#' erroneous formula given in the caption of Figure 4 in Fortune et al (2012).
#' (The error, an exchange of slope and intercept, was confirmed by
#' S. Fortune in an email to D. Kelley dated June 22, 2018.)
#'
#' * `"lockyer1976"` uses formulae from Table 1 of Lockyer (1976). The
#' permitted `species` and the formulae used are as follows (note that
#' the `"Gray Whale"` formula is in the table's caption, not in the table itself).
#'     * `"Blue Whale"`:       \eqn{2.899 L^{3.25}}{2.899 * L^3.25}
#'     * `"Bryde Whale"`:      \eqn{12.965 L^{2.74}}{12.965 * L^2.74}
#'     * `"Fin Whale"`:        \eqn{7.996 L^{2.90}}{7.996 * L^2.90}
#'     * `"Gray Whale"`:       \eqn{5.4 L^{3.28}}{5.4 * L^3.28}
#'     * `"Humpback Whale"`:   \eqn{16.473 L^{2.95}}{16.473 * L^2.95}
#'     * `"Minke Whale"`:      \eqn{49.574 L^{2.31}}{49.574 * L^2.31}
#'     * `"Pac. Right Whale"`: \eqn{13.200 L^{3.06}}{13.200 * L^3.06}
#'     * `"Sei Whale"`:        \eqn{25.763 L^{2.43}}{25.763 * L^2.43}
#'     * `"Sperm Whale"`:      \eqn{6.648 L^{3.18}}{6.648 * L^3.18}
#'
#' @param L whale length in m.
#'
#' @param species character value specifying the species (see \dQuote{Details}).
#' If  only one value is given, then it will repeated to have the same length
#' as `L`. Otherwise, its length must match the length of `L`.
#'
#' @param model either NULL (the default), to choose a model based on the
#' particular species, or a character value specifying the model (see \dQuote{Details}).
#' If  only one value is given, then it will repeated to have the same length
#' as `L`. Otherwise, its length must match the length of `L`.
#'
#' @return Mass in kg.
#'
#' @examples
#' library(whalestrike)
#' L <- seq(5, 15, length.out = 100)
#' kpt <- 1000 # kg per tonne
#' # Demonstrate (with dashing) the sensitivity involved in the single-digit
#' # parameter in Moore's formula, and (with colour) the difference to the
#' # Fortune et al. (2012) formulae.
#' plot(L, whaleMassFromLength(L, model = "moore2005") / kpt,
#'     type = "l", lwd = 2,
#'     xlab = "Right-whale Length [m]", ylab = "Mass [tonne]"
#' )
#' lines(L, 242.988 * exp(0.35 * L) / kpt, lty = "dotted", lwd = 2)
#' lines(L, 242.988 * exp(0.45 * L) / kpt, lty = "dashed", lwd = 2)
#' lines(L, whaleMassFromLength(L,
#'     species = "N. Atl. Right Whale",
#'     model = "fortune2012"
#' ) / kpt, col = 2, lwd = 2)
#' lines(L, whaleMassFromLength(L,
#'     species = "N. Pac. Right Whale",
#'     model = "fortune2012"
#' ) / kpt, col = 3, lwd = 2)
#' legend("topleft",
#'     lwd = 2, col = 1:3,
#'     legend = c("moore2005", "fortune2012 Atlantic", "fortune2012 Pacific")
#' )
#'
#' # Emulate Figure 1 of Lockyer (1976), with roughly-chosen plot limits.
#' L <- seq(0, 18, 0.5)
#' m <- whaleMassFromLength(L, species = "Pac. Right Whale", model = "lockyer1976") / kpt
#' plot(L, m,
#'     col = 1, xlab = "Length [m]", ylab = "Mass [tonne]", type = "l", lwd = 2,
#'     xaxs = "i", yaxs = "i", xlim = c(3, 30), ylim = c(0, 180)
#' )
#' L <- seq(0, 28, 0.5)
#' m <- whaleMassFromLength(L, species = "Blue Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 2, lwd = 2)
#' L <- seq(0, 24, 0.5)
#' m <- whaleMassFromLength(L, species = "Fin Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 3, lwd = 2)
#' L <- seq(0, 18, 0.5)
#' m <- whaleMassFromLength(L, species = "Sei Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 1, lty = 2, lwd = 2)
#' L <- seq(0, 17, 0.5)
#' m <- whaleMassFromLength(L, species = "Bryde Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 2, lty = 2, lwd = 2)
#' L <- seq(0, 12, 0.5)
#' m <- whaleMassFromLength(L, species = "Minke Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 3, lty = 2, lwd = 2)
#' L <- seq(0, 17, 0.5)
#' m <- whaleMassFromLength(L, species = "Humpback Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 1, lty = 3, lwd = 2)
#' L <- seq(0, 18, 0.5)
#' m <- whaleMassFromLength(L, species = "Sperm Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 2, lty = 3, lwd = 2)
#' L <- seq(0, 15, 0.5)
#' m <- whaleMassFromLength(L, species = "Gray Whale", model = "lockyer1976") / kpt
#' lines(L, m, col = 3, lty = 3, lwd = 2)
#' grid()
#' legend("topleft",
#'     col = c(1:3, 1:3, 1:2), lwd = 2, lty = c(rep(1, 3), rep(2, 3), rep(3, 3)),
#'     legend = c("Right", "Blue", "Fin", "Sei", "Bryde", "Minke", "Humpback", "Sperm", "Gray")
#' )
#'
#' @references
#' * Lockyer, C. "Body Weights of Some Species of Large Whales." J. Cons. Int.
#' Explor. Mer. 36, no. 3 (1976): 259-73.
#'
#' * Moore, M.J., A.R. Knowlton, S.D. Kraus, W.A. McLellan, and R.K. Bonde.
#' "Morphometry, Gross Morphology and Available Histopathology in North Atlantic
#' Right Whale (Eubalaena Glacialis) Mortalities (1970 to 2002)." Journal of
#' Cetacean Research and Management 6, no. 3 (2005): 199-214.
#'
#' * Fortune, Sarah M. E., Andrew W. Trites, Wayne L. Perryman, Michael J. Moore,
#' Heather M. Pettis, and Morgan S. Lynn. "Growth and Rapid Early Development of
#' North Atlantic Right Whales (Eubalaena Glacialis)." Journal of Mammalogy 93,
#' no. 5 (2012): 1342-54. \doi{10.1644/11-MAMM-A-297.1}.
#'
#' @seealso [whaleLengthFromMass()] is the reverse of this.
#'
#' @author Dan Kelley
#'
#' @export
whaleMassFromLength <- function(L, species = "N. Atl. Right Whale", model = NULL) {
    n <- length(species)
    if (length(L) < n) {
        L <- rep(L, n)
    }
    if (length(species) < length(L)) {
        species <- rep(species, length(L))
    }
    n <- length(species)
    if (is.null(model)) {
        model <- unlist(
            sapply(
                species,
                function(s) {
                    switch(s,
                        "N. Atl. Right Whale" = "fortune2012",
                        "Blue Whale" = "lockyer1976",
                        "Bryde Whale" = "lockyer1976",
                        "Fin Whale" = "lockyer1976",
                        "Gray Whale" = "lockyer1976",
                        "Humpback Whale" = "lockyer1976",
                        "Minke Whale" = "lockyer1976",
                        "Pac. Right Whale" = "lockyer1976",
                        "Sei Whale" = "lockyer1976",
                        "Sperm Whale" = "lockyer1976"
                    )
                }
            )
        )
    }
    #print(data.frame(species = species, model = model))
    if (length(model) == 1) {
        model <- rep(model, n)
    }
    if (length(species) == 1) {
        species <- rep(species, n)
    }
    if (n != length(model)) {
        stop("length of species (", n, ") does not equal length of model (", length(model), ")")
    }
    if (n != length(L)) {
        stop("length of species (", n, ") does not equal length of L (", length(L), ")")
    }
    rval <- rep(NA, n)
    for (i in 1:n) {
        if (model[i] == "moore2005") {
            if (species[i] == "N. Atl. Right Whale") {
                rval[i] <- 242.988 * exp(0.4 * L[i])
            } else {
                stop("The 'moore2005' model only works if species is 'N. Atl. Right Whale'")
            }
        } else if (model[i] == "fortune2012") {
            if (species[i] == "N. Atl. Right Whale") {
                rval[i] <- exp(-10.095 + 2.825 * log(100 * L[i]))
            } else if (species[i] == "N. Pac. Right Whale") {
                rval[i] <- exp(-12.286 + 3.158 * log(100 * L[i]))
            } else {
                stop("The 'fortune2012' model only works if species is 'N. Atl. Right Whale' or 'N. Pac. Right Whale'")
            }
        } else if (model[i] == "lockyer1976") {
            if (species[i] == "Blue Whale") {
                rval[i] <- 2.899 * L[i]^(3.25)
            } else if (species[i] == "Bryde Whale") {
                rval[i] <- 12.965 * L[i]^(2.74)
            } else if (species[i] == "Fin Whale") {
                rval[i] <- 7.996 * L[i]^(2.90)
            } else if (species[i] == "Gray Whale") {
                rval[i] <- 5.4 * L[i]^(3.28)
            } else if (species[i] == "Humpback Whale") {
                rval[i] <- 16.473 * L[i]^(2.95)
            } else if (species[i] == "Minke Whale") {
                rval[i] <- 49.574 * L[i]^(2.31)
            } else if (species[i] == "Pac. Right Whale") {
                rval[i] <- 13.200 * L[i]^(3.06)
            } else if (species[i] == "Sei Whale") {
                rval[i] <- 25.763 * L[i]^(2.43)
            } else if (species[i] == "Sperm Whale") {
                rval[i] <- 6.648 * L[i]^(3.18)
            } else {
                stop(
                    "species[", i, "]=\"", species[i], "\" must be one of the following:",
                    "\"Blue Whale\", \"Bryde Whale\", \"Fin Whale\", \"Gray Whale\", \"Humpback Whale\"",
                    "\"Minke Whale\", \"Pac. Right Whale\", \"Sei Whale\", or \"Sperm Whale\""
                )
            }
        } else {
            stop(
                "model[", i, "]=", model[i], "\" is unknown. This must be one of: ",
                "\"moore2005\", \"fortune2012\", or \"lockyer1976]\""
            )
        }
    }
    rval
}

#' Compute whale length from mass
#'
#' This works by inverting [whaleMassFromLength()] using [uniroot()].
#'
#' @param M Whale mass (kg).
#'
#' @param species A string indicating the whale species
#' (see [whaleMassFromLength()] for details).
#'
#' @param model Character string specifying the model
#' (see [whaleMassFromLength()] for details).
#'
#' @return Whale length (m).
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom stats uniroot
#'
#' @seealso [whaleMassFromLength()] is the reverse of this.
whaleLengthFromMass <- function(M, species = "N. Atl. Right Whale", model = "fortune2012") {
    n <- length(M)
    if (length(species) == 1) {
        species <- rep(species, n)
    }
    if (length(model) == 1) {
        model <- rep(model, n)
    }
    if (n != length(model)) {
        stop("length of M (", n, ") does not equal length of model (", length(model), ")")
    }
    if (n != length(species)) {
        stop("length of M (", n, ") does not equal length of species (", length(species), ")")
    }
    rval <- rep(NA, n)
    for (i in seq_along(M)) {
        rval[i] <- uniroot(
            function(x) {
                M[i] - whaleMassFromLength(x, species = species[i], model = model[i])
            },
            c(0.1, 100)
        )$root
    }
    rval
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
shipAreaFromMass <- function(ms) {
    length <- 11.73 # m
    beam <- 4.63 # m
    draft <- 1.58 # m
    displacement <- 20.46e3 # m^3
    factor <- (ms / displacement)^(1 / 3) # lengthscale factor
    length * (beam + 2 * draft) * factor^2
}
