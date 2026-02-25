# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

whaleMeasurementsTable <- structure(list(species = c(
    "Blue Whale", "Bryde Whale", "Fin Whale",
    "Gray Whale", "Humpback Whale", "Minke Whale", "N. Atl. Right Whale",
    "Pac. Right Whale", "Sei Whale", "Sperm Whale"
), properName = c(
    "Balaenoptera musculus",
    "Balaenoptera brydei", "Balaenoptera physalus", "Eschrichtius robustus",
    "Megaptera novaengliae", "Balaenoptera acutorostrata", "Eubalaena glacialis",
    "Eubalaena japonica", "Balaenoptera borealis", "Physeter macrocephalus"
), length = c(
    21.9, 13.4, 16.9, 12.5, 11.6, 6.7, 13.8, 13.8,
    13.4, 11.8
), bone = c(
    0.173, 0.093, 0.144, 0.136, 0.142, 0.078,
    0.143, 0.143, 0.093, 0.164
), sublayer = c(
    1.687, 0.631, 0.714,
    0.847, 1.372, 0.313, 1.325, 1.325, 0.631, 0.81
), blubber = c(
    0.076,
    0.044, 0.058, 0.097, 0.082, 0.033, 0.163, 0.163, 0.044, 0.115
), skin = c(
    0.004, 0.002, 0.005, 0.01, 0.009, 0.003, 0.009, 0.009,
    0.002, 0.004
)), row.names = c(NA, -10L), class = "data.frame")

#' Get values for various whale measurements
#'
#' This uses a data frame containing information about several whale
#' species, compiled by Alexandra Mayette and provided to Dan
#' Kelley as a personal communication on 2026-01-28.
#'
#' There are two species in the table that are not in Mayette's table.
#' These are `"Pac. Right Whale"` and `"Bryde Whale"`. For these, Mayette has suggesting
#' using values for the `"N. Atl. Right Whale"` and `"Sei Whale"` cases, respectively,
#' as conditional estimates for use in this package.
#'
#' @param species either (1) the name of a species or (2) NULL. In the first
#' case, the table is consulted to find a row with the given species
#' name, and that row is returned. In the second case, the whole table
#' is returned.
#'
#' @return The return value contains
#' * `name` species name, as used in e.g. whaleMassFromLength().
#' * `Species` proper species name. (This is not used in this package.)
#' * `length` whale length in metres.
#' * `bone` whale bone thickness in metres, measured from the centre to the sublayer.
#' * `sublayer` thickness of sublayer in meters; this was called `muscle` in Mayette's document.
#' * `blubber` whale blubber thickness in meters.
#' * `skin` whale skin thickness in metres.
#'
#' @examples
#' library(whalestrike)
#' # All species in database
#' whaleMeasurements()
#' # A particular species
#' whaleMeasurements("N. Atl. Right Whale")
#'
#' @export
#'
#'
#' @author Dan Kelley, using data and advice from Alexandra Mayette
#'
#' @references
#' Mayette, A. (2026). Measurements of large whale tissue thickness (Data set).
#' Zenodo. \doi{10.5281/zenodo.18764979}.
whaleMeasurements <- function(species = NULL) {
    if (is.null(species)) {
        whaleMeasurementsTable
    } else {
        i <- which(species == whaleMeasurementsTable$species)
        if (length(i) != 1) {
            stop("cannot find measurements for species '", species, "'")
        }
        whaleMeasurementsTable[i, ]
    }
}

#' Set parameters for a whale strike simulation
#'
#' Assembles control parameters into a list suitable for passing to [strike()]
#' and the functions that it calls. If `file` is provided, then all the other
#' arguments are read from that source. Note that [updateParameters()] may
#' be used to modify the results of `parameters`, e.g. for use in sensitivity
#' tests.
#'
#' @param ms Ship mass (kg).
#'
#' @param Ss Ship wetted area (m^2). This, together with `Cs`, is used by
#' [shipWaterForce()] to estimate ship drag force. If `Ss`
#' is not given, then an estimate is made by calling [shipAreaFromMass()] with
#' the provided value of `ms`.
#'
#' @param Ly Ship impact horizontal extent (m); defaults to 1.15m if not specified,
#' based on an analysis of the shape of the bow of typical coastal fishing boats
#' of the Cape Islander variety.
#'
#' @param Lz Ship impact vertical extent (m); defaults to 1.15m if not specified,
#' based on the same analysis as for Ly.
#'
#' @param lw either (1) whale length in metres or (2) the string `"from_species"`.
#' If the latter, then the length is determined from [whaleMeasurements()].
#' In either case, the length is used by [whaleAreaFromLength()] to
#' calculate area, which is needed for the water drag calculation done by
#' [whaleWaterForce()].
#'
#' @param species a string indicating the whale species. For the permitted values,
#' see [whaleMassFromLength()]. (The `species` value can also set the
#' `lw` and `l` values, as noted in their portions of this documentation.)
#'
#' @param mw either (1) the whale mass in kg or (2) NULL. In the latter case,
#' the mass is calculated from whale length, using [whaleMassFromLength()]
#' with `type="wetted"`.
#'
#' @param Sw either (1) the whale surface area in m^2 or (2) NULL. If the
#' latter case, the area is calculated from whale length using
#' [whaleAreaFromLength()].
#'
#' @param l either (1) a numerical vector of length 4 that indicates
#' the thicknesses in metres of skin, blubber, sublayer and bone; (2) NULL
#' to set these four values to 0.025, 0.16, 1.12, and 0.1; or (3) the
#' string `"from_species"`, in which case these four values are
#' determined by calling [whaleMeasurements()].
#' The default skin thickness of 0.025 m represents the 0.9-1.0 inch value
#' stated in Section 2.2.3 of Raymond (2007).
#' The blubber default of 0.16 m is a rounded average of the values inferred
#' by whale necropsy, reported in Appendix 2 of Daoust et al., 2018.
# > round(mean(c(17,14,18.13,18,21.25,16.75,13.33,7)/100),2)
# [1] 0.16
#' The sublayer default of 1.12 m may be reasonable at some spots on the whale body.
#' The bone default of 0.1 m may be reasonable at some spots on the whale body.
#' The sum of these default values, 1.40 m, is a whale radius that
#' is consistent with a half-circumference of 4.4 m, reported in Table 2.2
#' of Raymond (2007).  Note, however, that these values are not identical
#' to those found in `whaleMeasurements`.
#'
#' @param a,b Numerical vectors of length 4, giving values to use in the
#' stress-strain law `stress=a*(exp(b*strain)-1)`, where `a` is in Pa
#' and `b` is unitless. By construction, `a*b` is the local modulus at
#' low strain (i.e. at low `b*strain` values), and that `b` is the
#' efolding scale for nonlinear increase in stress with strain.
#' This exponential relationship has been mapped out
#' for whale blubber, using a curve fit to Figure 2.13 of Raymond (2007), and
#' these values are used for the second layer (blubber); see
#' the documentation for the [raymond2007] dataset, to see
#' for how that fit was done.
#' If not provided, `a` defaults to
#' `c(17.8e6/0.1, 1.58e5, 1.58e5, 8.54e8/0.1)`
#' and `b` defaults to
#' `c(0.1, 2.54, 2.54, 0.1)`.
#' The skin defaults are set up to give a linear shape (since `b` is small)
#' with the `a*b` product
#' being 17.8e6 Pa, which is the adult-seal value
#' given in Table 3 of Grear et al. (2017).
#' The blubber defaults are from a regression of the stress-strain
#' relationship shown in Figure 2.13 of Raymond (2007).
#' The sublayer defaults are set to match those of blubber, lacking
#' any other information.
#' The bone default for `b` is small, to set up a linear function,
#' and `a*b` is set to equal 8.54e8 Pa,
#' given in Table 2.3 of Raymond (2007) and Table 4.5 of
#' Campbell-Malone (2007).
#'
#' @param s Numerical vector of length 4, giving the ultimate strengths (Pa) of
#' skin, blubber, sublayer, and bone, respectively. If not provided, the
#' value is set to `1e6 * c(19.600,0.255,0.255,22.900)`
#' with reasoning as follows.
#' The skin default of 19.6 MPa
#' is a rounded value from Table 3 of Grear et al. (2018) for adult seal skin strength at
#' an orientation of 0 degrees.  The blubber and sublayer values were chosen
#' as the central point of a logistic fit of whale collision damage
#' to maximal stress during a default impact simulation.
#' (For comparison, a strength of
#' 0.437 MPa may be inferred by
#' multiplying Raymond's (2007) Figure 2.13 elastic modulus of 0.636 MPa
#' by the ratio 0.97/1.41 determined for adult seal strength/modulus, as reported
#' in Table 3 of Grear et al. (2018).)
#' The bone default o 22.9 MPa is from Table 2.3 of Raymond (2007) and
#' Table 4.5 of Campbell-Malone (2007).
#'
#' @param theta Whale skin deformation angle (deg); defaults to 55 degrees,
#' if not supplied, because that angle produces a good match to Raymond's (2007)
#' Figure 6.1 for the total force as a function of vessel speed, for large
#' vessels. Note that the match works almost as well in the range 50 deg
#' to 70 deg.
#'
#' @param Cs Drag coefficient for ship (dimensionless),
#' used by [shipWaterForce()] to estimate ship drag force. Defaults
#' to 1e-2, which is 4 times the frictional coefficient of 2.5e-3
#' inferred from Figure 4 of Manen and van Oossanen (1988), assuming
#' a Reynolds number of 5e7, computed from speed 5m/s, lengthscale 10m
#' and viscosity 1e-6 m^2/s. The factor of 4 is under the assumption
#' that frictional drag is about a quarter of total drag.
#' The drag force is computed with [shipWaterForce()].
#'
#' @param Cw Drag coefficient for whale (dimensionless),
#' used by [whaleWaterForce()] to estimate whale drag force.
#' Defaults to 2.5e-3, for Reynolds number 2e7, computed from speed
#' 2 m/s, lengthscale 5m which is chosen to be between radius and length, and
#' viscosity 1e-6 m^2/s.  The drag force is computed with
#' [whaleWaterForce()].
#'
#' @param logistic a [list] containing `logStressCenter` and `logStressWidth`,
#' which define an empirical logistic fit of an index of whale injury in
#' observed strikes (ranging from 0 for no injury to 1 for fatal injury),
#' as a function of the base-10 logarithm of compressive
#' stress, as well as `tau25`, `tau50` and `tau75`, which are the stresses
#' in that fit that yield index values of 0.25, 0.50 and 0.75, respectively;
#' these values set colour boundaries in [plot.strike()] plots that have
#' `which="threat"`.
#'
#' @param file Optional name a comma-separated file that holds all of the
#' previous values, except `Cs` and `Cw`. If provided,
#' then other parameters except `Cs` and `Cw` are
#' ignored, because values are sought from the file. The purpose of
#' this is in shiny apps that want to save a simulation framework.
#' The file should be saved [write.csv()] with
#' `row.names=FALSE`.
#'
#' @return
#' A named list holding the parameters, with defaults and alternatives reconciled
#' according to the system described above, along with some items used internally,
#' including `lsum`, which is the sum of the values in `l`, and `stressFromStrain()`,
#' a function created by [stressFromStrainFunction()] that computes compression
#' force from engineering strain.
#'
#' @examples
#' parms <- parameters()
#' epsilon <- seq(0, 1, length.out = 100) # strain
#' sigma <- parms$stressFromStrain(epsilon) # stress
#' plot(epsilon, log10(sigma), xlab = "Strain", ylab = "log10(Stress [MPa])", type = "l")
#' mtext("Note sudden increase in stress, when bone compression starts")
#'
#' @author Dan Kelley
#'
#' @references
#'
#' Campbell-Malone, Regina. "Biomechanics of North Atlantic Right Whale Bone:
#' Mandibular Fracture as a Fatal Endpoint for Blunt Vessel-Whale Collision
#' Modeling." PhD Thesis, Massachusetts Institute of Technology and Woods Hole
#' Oceanographic Institution, 2007. \doi{10.1575/1912/1817}.
#'
#' Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' "Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine Animal Response
#' Society, and Fisheries and Oceans Canada, 2018.
#' \url{https://publications.gc.ca/site/eng/9.850838/publication.html}.
#'
#' Grear, Molly E., Michael R. Motley, Stephanie B. Crofts, Amanda E. Witt, Adam
#' P. Summers, and Petra Ditsche. "Mechanical Properties of Harbor Seal Skin and
#' Blubber - a Test of Anisotropy." Zoology 126 (2018): 137-44.
#' \doi{10.1016/j.zool.2017.11.002}.
#'
#' Raymond, J. J. "Development of a Numerical Model to Predict Impact Forces on a
#' North Atlantic Right Whale during Collision with a Vessel." University of New
#' Hampshire, 2007. \url{https://scholars.unh.edu/thesis/309/}.
#'
#' @export
#'
#' @importFrom utils read.csv
parameters <- function(
    ms = 45e3, Ss = NULL, Ly = 1.15, Lz = 1.15,
    species = "N. Atl. Right Whale", lw = 13.7, mw = NULL, Sw = NULL,
    l = NULL, a = NULL, b = NULL, s = NULL,
    theta = 55,
    Cs = 0.01, Cw = 0.0025,
    logistic = list(
        logStressCenter = 5.38, logStressWidth = 0.349,
        tau25 = 0.100e6, tau50 = 0.241e6, tau75 = 0.581e6
    ),
    file = NULL) {
    if (!is.null(file)) {
        rval <- as.list(read.csv(file))
        rval$Ss <- shipAreaFromMass(rval$ms)
        rval$mw <- whaleMassFromLength(rval$lw, species = species)
        rval$Sw <- whaleAreaFromLength(rval$lw, species = species, type = "wetted")
        rval$tmax <- NULL
        rval$vs <- NULL
        rval$Cs <- Cs
        rval$Cw <- Cw
        rval$l <- c(rval$l1, rval$l2, rval$l3, rval$l4)
        rval$l1 <- rval$l2 <- rval$l3 <- rval$l4 <- NULL
        rval$lsum <- sum(rval$l)
        # the next are copied from below. The app doesn't let the user
        # set these things, so we know their values.
        # NOTE: keep in synch with 'BBBB' below!
        rval$a <- c(17.8e6 / 0.1, 1.58e5, 1.58e5, 8.54e8 / 0.1)
        rval$b <- c(0.1, 2.54, 2.54, 0.1)
        rval$s <- c(19.6e6, 0.437e6, 0.437e6, 22.9e6)
        o <- sort(names(rval))
        rval <- rval[o]
    } else {
        # Check some elements, and set some defaults
        if (length(ms) != 1) {
            stop("ms must be a single numeric value")
        }
        if (ms <= 0) {
            stop("ms must be positive, but it is ", ms)
        }
        if (is.null(Ss)) {
            Ss <- shipAreaFromMass(ms)
        }
        if (length(Ss) != 1) {
            stop("Ss must be a single numeric value")
        }
        if (Ss <= 0) {
            stop("Ss must be positive, but it is ", Ss)
        }
        if (length(Ly) != 1) {
            stop("Ly must be a single numeric value")
        }
        if (Ly <= 0) {
            stop("Ly must be positive, but it is ", Ly)
        }
        if (length(Lz) != 1) {
            stop("Lz must be a single numeric value")
        }
        if (Lz <= 0) {
            stop("Lz must be positive, but it is ", Lz)
        }
        if (length(lw) != 1) {
            stop("lw must be a string or a single numeric value")
        }
        if (identical(lw, "from_species")) {
            lw <- whaleMeasurements(species)$length
        }
        if (lw <= 0) {
            stop("lw must be positive, but it is ", lw)
        }
        if (is.null(mw)) {
            mw <- whaleMassFromLength(lw, species = species)
        }
        if (length(mw) != 1) {
            stop("cannot handle more than one 'mw' at a time")
        }
        if (is.null(Sw)) {
            Sw <- whaleAreaFromLength(lw, species = species, type = "wetted")
        }
        if (length(Sw) != 1) {
            stop("cannot handle more than one 'Sw' at a time")
        }
        if (is.null(l)) {
            l <- c(0.025, 0.16, 1.12, 0.1)
        }
        if (is.null(a)) {
            a <- c(17.8e6 / 0.1, 1.58e5, 1.58e5, 8.54e8 / 0.1)
        }
        if (is.null(b)) {
            b <- c(0.1, 2.54, 2.54, 0.1)
        }
        if (is.null(s)) {
            s <- 1e6 * c(19.600, 0.255, 0.255, 22.900)
        }
        if (any(s <= 0) || length(s) != 4) {
            stop("'s' must be a vector with 4 positive numbers")
        }
        if (identical(l, "from_species")) {
            tmp <- whaleMeasurements(species)
            l <- c(tmp$skin, tmp$blubber, tmp$sublayer, tmp$bone)
            # message("l: ", paste(l, collapse = " "))
        }
        if (any(l <= 0) || length(l) != 4) {
            stop("'l' must be a vector with 4 positive numbers (even if 'from_species')")
        }
        if (any(a <= 0) || length(a) != 4) {
            stop("'a' must be a vector with 4 positive numbers")
        }
        if (any(b <= 0) || length(b) != 4) {
            stop("'b' must be a vector with 4 positive numbers")
        }
        if (length(theta) != 1) {
            stop("cannot handle more than one 'theta' at a time")
        }
        if (theta < 0 || theta > 89) {
            stop("whale skin deformation angle (theta) must be between 0 and 89 deg, but it is ", theta)
        }
        if (length(Cs) != 1) {
            stop("cannot handle more than one 'Cs' at a time")
        }
        if (Cs <= 0) {
            stop("ship resistance parameter (Cs) must be positive, but it is ", Cs)
        }
        if (length(Cw) != 1) {
            stop("cannot handle more than one 'Cw' at a time")
        }
        if (Cw <= 0) {
            stop("ship resistance parameter (Cw) must be positive, but it is ", Cw)
        }
        rval <- list(
            ms = ms, Ss = Ss,
            Ly = Ly, Lz = Lz,
            mw = mw, Sw = Sw, lw = lw,
            l = l, lsum = sum(l), a = a, b = b, s = s,
            theta = theta,
            Cs = Cs, Cw = Cw
        )
    }
    # For efficiency, create and store an overall stress-strain function
    rval$stressFromStrain <- stressFromStrainFunction(rval$l, rval$a, rval$b)
    rval$logistic <- logistic
    class(rval) <- "parameters"
    rval
}

#' Summarize a parameters object
#'
#' This provides an overview of the contents of an object
#' created with [parameters()].
#'
#' @param object an object of class `"parameters"`, as created with [parameters()].
#'
#' @param \dots ignored
#'
#' @examples
#' summary(parameters())
#'
#' @author Dan Kelley
#'
#' @return `summary.parameters` returns nothing. It is called for
#' its side effect of printing information about the parameters
#' in a ship-whale collision simulation.
#'
#' @export
summary.parameters <- function(object, ...) {
    cat("Whale and ship properties, as created by parameters():\n")
    cat("  Ship properties:\n")
    cat(sprintf("    ms:   %12g kg  -- mass\n", object$ms))
    cat(sprintf("    Ss:   %12g m^2 -- wetted area\n", object$Ss))
    cat(sprintf("    Ly:   %12g m   -- width of impact area\n", object$Ly))
    cat(sprintf("    Lz:   %12g m   -- height of impact area\n", object$Lz))
    cat(sprintf("    Cs:   %12g     -- drag coefficient\n", object$Cs))
    cat("  Whale properties:\n")
    cat(sprintf("    mw:   %12g kg  -- mass\n", object$mw))
    cat(sprintf("    Sw:   %12g m^2 -- wetted area\n", object$Sw))
    cat(sprintf("    lw:   %12g m   -- length\n", object$Lw))
    cat(sprintf("    l[1]: %12g m   -- thickness of skin\n", object$l[1]))
    cat(sprintf("    l[2]: %12g m   -- thickness of blubber\n", object$l[2]))
    cat(sprintf("    l[3]: %12g m   -- thickness of sublayer\n", object$l[3]))
    cat(sprintf("    l[4]: %12g m   -- half-thickness of bone\n", object$l[4]))
    cat(sprintf("    lsum: %12g m   -- a[1]+a[2]+a[3]+a[4]\n", object$lsum))
    cat(sprintf("    a[1]: %12.3e Pa  -- skin stress factor\n", object$a[1]))
    cat(sprintf("    a[2]: %12.3e Pa  -- blubber stress factor\n", object$a[2]))
    cat(sprintf("    a[3]: %12.3e Pa  -- sublayer stress factor\n", object$a[3]))
    cat(sprintf("    a[4]: %12.3e Pa  -- bone stress factor\n", object$a[4]))
    cat(sprintf("    b[1]: %12g     -- skin stress nonlinearity term\n", object$b[1]))
    cat(sprintf("    b[2]: %12g     -- blubber stress nonlinearity term\n", object$b[2]))
    cat(sprintf("    b[3]: %12g     -- sublayer stress nonlinearity term\n", object$b[3]))
    cat(sprintf("    b[4]: %12g     -- bone stress nonlinearity term\n", object$b[4]))
    cat(sprintf("    s[1]: %12.3e Pa  -- skin strength\n", object$s[1]))
    cat(sprintf("    s[2]: %12.3e Pa  -- blubber strength\n", object$s[2]))
    cat(sprintf("    s[3]: %12.3e Pa  -- sublayer strength\n", object$s[3]))
    cat(sprintf("    s[4]: %12.3e Pa  -- bone strength\n", object$s[4]))
    cat(sprintf("    theta: %11g deg -- impact dimple angle\n", object$theta))
    cat(sprintf("    Cw: %12g       -- drag coefficient\n", object$Cw))
    cat("  Functions:\n")
    cat("    stressFromStrain()     -- function to compute stress\n")
    cat("    logistic()             -- function to compute lethality index\n")
}


#' Update parameters
#'
#' `updateParameters()` is used to alter one or more components of an existing
#' object of type `"parameters"` that was created by [parameters()]. This
#' can be useful for e.g. sensitivity tests (see \dQuote{Details}).
#'
#' Two important differences between argument handling in `updateParameters()`
#' and [parameters()] should be kept in mind.
#'
#' First, `updateParameters()` does not check its arguments
#' for feasible values.  This can lead to bad results when using
#' [strike()], which is e.g. expecting four layer thicknesses to
#' be specified, and also that each thickness is positive.
#'
#' Second, `updateParameters()` does not perform ancillary
#' actions that [parameters()] performs, with regard to certain interlinking
#' argument values.  Such actions are set up for
#' whale length and ship mass, which are easily-observed
#' quantities from other quantities can be estimated
#' using `whalestrike` functions.  If `lw` (whale length) is
#' supplied to [parameters()] without also supplying `mw`
#' (whale mass), then [parameters()] uses [whaleMassFromLength()]
#' to infer `mw` from `lw`. The same procedure is used to infer
#' `Sw` if it is not given, using [whaleAreaFromLength()].
#' Similarly, [parameters()] uses [shipAreaFromMass()] to
#' compute `Ss` (ship area) from `ms` (ship mass), if the `Ss`
#' argument is not given.  Importantly, these three inferences
#' are *not* made by `updateParameters()`, which alters only
#' those values that are supplied explicitly. It is easy
#' to supply those values, however; for example,
#' ```
#' parms <- updateParameters(PARMS, lw=1.01 * PARMS$lw))
#' parms <- updateParameters(parms, mw=whaleMassFromLength(parms$lw))
#' parms <- updateParameters(parms, Sw=whaleAreaFromLength(parms$lw))
#' ```
#' modifies a base state stored in `PARMS`, increasing whale length
#' by 1% and then increasing whale mass and area accordingly.  This
#' code block is excerpted from a sensitivity test of the model, in
#' which
#' ```
#' parms <- updateParameters(PARMS, ms=1.01 * PARMS$ms)
#' parms <- updateParameters(parms, Ss=shipAreaFromMass(parms$ms))
#' ```
#' was also used to perturb ship mass (and inferred area).
#'
#' @param original An object of class `"parameters"`, as created by [parameters()]
#' and perhaps later altered by previous calls to `updateParameters()`.
#'
#' @inheritParams parameters
#'
#' @return A named list holding the items of the same name as those in the list
#' returned by [parameters()].
#'
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#'
#' @references
#'
#' Daoust, Pierre-Yves, Emilie L. Couture, Tonya Wimmer, and Laura Bourque.
#' "Incident Report. North Atlantic Right Whale Mortality Event in the Gulf of St.
#' Lawrence, 2017." Canadian Wildlife Health Cooperative, Marine Animal Response
#' Society, and Fisheries and Oceans Canada, 2018.
#' \url{https://publications.gc.ca/site/eng/9.850838/publication.html}.
#'
#' @author Dan Kelley
#'
#' @export
updateParameters <- function(original,
                             ms, Ss,
                             Ly, Lz,
                             species, lw, mw, Sw,
                             l, a, b, s, theta,
                             Cs, Cw,
                             logistic,
                             debug = 0) {
    rval <- original
    if (!missing(ms)) rval$ms <- ms
    if (!missing(Ss)) rval$Ss <- Ss
    if (!missing(Ly)) rval$Ly <- Ly
    if (!missing(Lz)) rval$Lz <- Lz
    if (!missing(species)) rval$species <- species
    if (!missing(lw)) rval$lw <- lw
    if (!missing(mw)) rval$mw <- mw
    if (!missing(Sw)) rval$Sw <- Sw
    if (!missing(l)) {
        rval$l <- l
        rval$lsum <- sum(l)
    }
    if (!missing(a)) rval$a <- a
    if (!missing(b)) rval$b <- b
    if (!missing(s)) rval$s <- s
    if (!missing(theta)) rval$theta <- theta
    if (!missing(Cs)) rval$Cs <- Cs
    if (!missing(Cw)) rval$Cw <- Cw
    if (!missing(logistic)) rval$logistic <- logistic
    if (debug > 0) {
        cat("at end of updateParameters(): a=", paste(rval$a, collapse = " "), "\n")
        cat("at end of updateParameters(): b=", paste(rval$b, collapse = " "), "\n")
        cat("at end of updateParameters(): l=", paste(rval$l, collapse = " "), "\n")
    }
    rval$stressFromStrain <- stressFromStrainFunction(rval$l, rval$a, rval$b)
    class(rval) <- "parameters"
    rval
}
