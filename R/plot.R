# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

#' Draw polygon between two xy curves
#'
#' This adds to an existing plot by filling the area between the
#' lower=lower(x) and upper=upper(x) curves.  In most cases, as
#' shown in \dQuote{Examples}, it is helpful
#' to use `xaxs="i"` in the preceding plot call, so that the
#' polygon reaches to the edge of the plot area.
#'
#' @param x Coordinate along horizontal axis
#'
#' @param lower Coordinates of the lower curve, of same length as `x`,
#' or a single value that gets repeated to the length of `x`.
#'
#' @param upper Coordinates of the upper curve, or a single value that gets
#' repeated to the length of `x`.
#'
#' @param ... passed to [polygon()]. In most cases, this
#' will contain `col`, the fill colour, and possibly `border`,
#' the border colour, although cross-hatching with `density`
#' and `angle` is also a good choice.
#'
#' @examples
#' # 1. CO2 record
#' plot(co2, xaxs = "i", yaxs = "i")
#' fillplot(time(co2), min(co2), co2, col = "pink")
#'
#' # 2. stack (summed y) plot
#' x <- seq(0, 1, 0.01)
#' lower <- x
#' upper <- 0.5 * (1 + sin(2 * pi * x / 0.2))
#' plot(range(x), range(lower, lower + upper),
#'     type = "n",
#'     xlab = "x", ylab = "y1, y1+y2",
#'     xaxs = "i", yaxs = "i"
#' )
#' fillplot(x, min(lower), lower, col = "darkgray")
#' fillplot(x, lower, lower + upper, col = "lightgray")
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom graphics polygon
fillplot <- function(x, lower, upper, ...) {
    n <- length(x)
    if (length(lower) == 1L) lower <- rep(lower, n)
    if (length(upper) == 1L) upper <- rep(upper, n)
    if (n != length(lower)) stop("lengths of x and lower must match")
    if (n != length(upper)) stop("lengths of x and upper must match")
    xx <- c(x, rev(x), x[1])
    yy <- c(upper, rev(lower), upper[1])
    polygon(xx, yy, ...)
}

fillplot4 <- function(x, y, yOffset = 0, breaks, col, ...) {
    if (missing(x)) stop("must give x")
    if (missing(y)) stop("must give y")
    nx <- length(x)
    if (length(y) != nx) stop("x and must be of equal length")
    if (missing(breaks)) stop("must give breaks")
    if (missing(col)) stop("must give col")
    nbreaks <- length(breaks)
    if (nbreaks != 3) stop("invalid use of non-exported function (programmr error)")
    if (length(col) != nbreaks + 1) stop("must have 1 more col than break")
    xx <- c(x, rev(x), x[1])
    yy <- c(y, rep(0, nx), y[1])
    polygon(xx, yy + yOffset, col = col[4], border = col[4])
    yy <- ifelse(yy < breaks[3], yy, breaks[3])
    # message("breaks[3]=", breaks[3], "; max(y)=", max(y), "; max(yy) after trim=", max(yy), "; col=", col[3])
    polygon(xx, yy + yOffset, col = col[3], border = col[3])
    yy <- ifelse(yy < breaks[2], yy, breaks[2])
    # message("breaks[2]=", breaks[3], "; max(y)=", max(y), "; max(yy) after trim=", max(yy), "; col=", col[2])
    polygon(xx, yy + yOffset, col = col[2], border = col[2])
    yy <- ifelse(yy < breaks[1], yy, breaks[1])
    # message("breaks[1]=", breaks[3], "; max(y)=", max(y), "; max(yy) after trim=", max(yy), "; col=", col[1])
    polygon(xx, yy + yOffset, col = col[1], border = col[1])
    lines(x, y + yOffset)
}

#' Plot a strike object
#'
#' Creates displays of various results of a simulation performed
#' with [strike()].
#'
#' @param x An object created by [strike()].
#'
#' @param which A character vector that indicates what to plot.
#' This choices for its entries are listed below, in no particular order.
#' \itemize{
#'
#' \item `"location"` for a time-series plot of boat location `xw` in
#' dashed black, whale centerline `xs` in solid gray,
#' blubber-interior interface in red, and skin in blue. The maximum
#' acceleration of ship and whale (in "g" units) are indicated in notes
#' placed near the horizontal axes. Those acceleration indications report
#' just a single value for each of ship and whale, but if the blubber
#' and sublayer have been squeezed to their limits, yielding a short and
#' intense force spike as the bone compresses, then the summaries will
#' also report on the spike duration and intensity. The spike is computed
#' based on using [runmed()] on the acceleration data, with a
#' `k` value that is set to correspond to 5 ms, or to k=11, whichever
#' is larger.
#'
#' \item `"section"` to plot skin thickness, blubber thickness and sublayer thickness
#' in one panel, creating a cross-section diagram.
#'
# \item \code{"injury"} a stacked plot showing time-series traces of health
# indicators for skin extension, blubber compression, and sublayer compression.
# These take the form of dotted horizontal lines, with labels above and at the
# left side of the panel.  During times when an injury criterion is halfway met,
# e.g. that blubber stress 1/2 the value of blubber strength, then the dotted
# line is overdrawn with a thick gray line. During times when the criterion
# is exceeded, the colour shifts to black.
#'
#' \item `"threat"` a stacked plot showing time-series traces of an
#' ad-hoc measure of the possible threat to skin, blubber and sublayer.
#' The threat level is computed as the ratio
#' of stress to ultimate strength, e.g. for blubber, it is
#' `x$WCF$stress/x$parms$s[2]`. The same vertical scale is used in
#' each of the subpanels that make up the stack. Any values exceeding
#' 10 are clipped to 10, and in such a case the overall label on the vertical
#' axis will note this clipping, although it should be easy to see, because
#' the way it most often occurs is if the soft layers "bottom out" onto
#' the bone, which yields a short period of very high stress, owing to the
#' very high compression modulus of bone. Each of the curves is filled in with a light
#' gray colour for stress/strength values up to 1, and with black
#' for higher values; this makes it easy to tell at a glance whether the
#' threat level is high.
#'
#' \item `"whale acceleration"` for a time-series plot of whale acceleration.
#'
#' \item `"blubber thickness"` for a time-series plot of blubber thickness.
#'
#' \item `"sublayer thickness"` for a time-series plot of the thickness
#' of the layer interior to the blubber.
#'
#' \item `"reactive forces"` for a time-series showing the reactive
#' forces associated with skin stretching (solid) and the compression of the
#' blubber and sublayer components (dashed).
#'
#' \item `"compression stress"` for a time-series plot of the compression stress on the blubber
#' and the layer to its interior. (These stresses are equal, under an equilibrium assumption.)
#'
#' \item `"skin stress"` for a time-series of skin stress in the along-skin y and z directions.
#'
#' \item `"lethality index"` for a time-series of Lethality Index, computed from compression stress
#' using [lethalityIndexFromStress()]. Values of Lethality Index that exceed 0.5 are highlighted,
#' and a dotted line is drawn at that value.
#'
#' \item `"values"` for a listing of `param` values.
#'
#' \item `"all"` for all of the above.
#'
#' \item `"default"` for a three-element plot showing `"location"`,
#' `"section"`, and `"threat"`.
#' }
#'
#' @param drawEvents Logical, indicating whether to draw lines for some events,
#' such as the moment of closest approach.
#'
#' @param colwcenter Colour used to indicate the whale centre.
#'
#' @param colwinterface Colour used to indicate the interface
#' between whale blubber and sublayer.
#'
#' @param colwskin Colour used to indicate the whale skin.
#'
#' @param cols As `colw`, but the colour to be used for the ship bow location,
#' which is drawn with a dashed line.
#'
# @param colInjury Two-element colour specification used in `"injury"`
# panels. The first colour is used to indicate values that are halfway to the
# injury crition, and the second is used to indicate values that exceed the
# criterion.
#'
#' @param colThreat a 4-element colour specification used in `"threat"` plots.
#' The colour transitions depend on the layer being plotted.
#' For skin and bone (i.e. in the bottom and top subpanels of the plot),
#' the first colour is used up to a stress:strength ratio
#' of 1/4, the second up to 1/2, the third up to 3/4, and the fourth above
#' 3/4.  For blubber and sublayer (i.e. the second and third subpanels,
#' counting from the bottom), the colours are defined by comparing compressive
#' stress to the values of `tau25`, `tau50` and `tau75` (see the [parameters()]
#' documentation), with the first colour for stress under `tau25`, the second
#' for stress up to `tau50`, the third for stress up to `tau75`, and the fourth
#' for higher stresses.
#'
#' @param lwd Line width used in plots for time intervals in which damage
#' criteria are not exceeded.
#'
#' @param D Factor by which to thicken lines during times during which damage
#' criteria are exceeded.
#'
#' @param debug Integer indicating debugging level, 0 for quiet operation and higher values
#' for more verbose monitoring of progress through the function.
#'
#' @param ... Ignored.
#'
#' @references
#' See [whalestrike()] for a list of references.
#'
#' @examples
#' # 1. default 3-panel plot
#' t <- seq(0, 0.7, length.out = 200)
#' state <- c(xs = -2, vs = knot2mps(12), xw = 0, vw = 0) # 12 knot ship
#' parms <- parameters() # default values
#' sol <- strike(t, state, parms)
#' par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0), mfrow = c(1, 3))
#' plot(sol)
#' # 2. all 12 plot types
#' par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0), mfrow = c(4, 3))
#' plot(sol, "all")
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom graphics abline axis box lines legend mtext par plot text
#' @importFrom grDevices hcl
#' @importFrom stats approx runmed
plot.strike <- function(x, which = "default", drawEvents = TRUE,
                        colwcenter = "black", # Slate Gray",
                        colwinterface = "black", # colwinterface="Firebrick",
                        colwskin = "black", # colwskin="Dodger Blue 4",
                        cols = "black",
                        colThreat = c("white", "lightgray", "darkgray", "black"),
                        lwd = 1, D = 3, debug = 0, ...) {
    g <- 9.8 # gravity
    t <- x$t
    xs <- x$xs
    vs <- x$vs
    xw <- x$xw
    vw <- x$vw
    death <- xs >= xw
    if (any(death)) {
        firstDead <- which(death)[1]
        dead <- firstDead:length(t)
        xw[dead] <- xw[firstDead]
        x$WCF$compressed[dead, 1] <- 0
        x$WCF$compressed[dead, 2] <- 0
    }
    showEvents <- function(xs, xw) {
        if (drawEvents) {
            death <- which(xs >= xw)[1]
            tdeath <- if (is.finite(death)) t[death] else NA
            if (is.finite(tdeath)) {
                abline(v = tdeath, lwd = lwd, col = "blue")
                mtext("Fatality", at = tdeath, side = 3, line = 0, col = "blue", font = 2, cex = par("cex"))
            }
        }
    }
    all <- "all" %in% which
    if (length(which) == 1 && which == "default") {
        which <- c("location", "section", "lethality index")
    }
    # Ensure that the plot type is known.
    allowed <- c(
        "all", "location", "section", "threat", "whale acceleration",
        "blubber thickness", "sublayer thickness",
        "whale water force", "reactive forces", "skin stress",
        "compression stress", "lethality index", "values"
    )
    for (w in which) {
        if (!(w %in% allowed) && !length(grep("NEW", w))) {
            stop(
                "which value \"", w, "\" is not allowed; try one of: \"",
                paste(allowed, collapse = "\" \""), "\""
            )
        }
    }
    # x(t) and xw(t)
    if (all || "location" %in% which) {
        ylim <- range(c(xs, xw), na.rm = TRUE)
        plot(t, xs, type = "l", xlab = "Time [s]", ylab = "Location [m]", col = cols, ylim = ylim, lwd = lwd, lty = "84", xaxs = "i")
        lines(t, xw, lwd = lwd, col = colwcenter)
        compressed <- x$WCF$compressed
        y <- xw - compressed[, 4]
        lines(t, y, col = colwinterface, lwd = lwd)
        y <- y - compressed[, 3]
        lines(t, y, col = colwinterface, lwd = lwd)
        y <- y - compressed[, 2]
        lines(t, y, col = colwskin, lwd = lwd)
        y <- y - compressed[, 1]
        lines(t, y, col = colwskin, lwd = lwd)
        # Accelerations (quite complicated; possibly too confusing to viewer)
        k <- round(0.005 / (t[2] - t[1]))
        k <- max(k, 11L)
        if (!(k %% 2)) {
            k <- k + 1
        }
        as <- derivative(vs, t)
        asmax <- max(abs(as))
        asmaxs <- max(abs(runmed(as, k))) # smoothed
        if (asmax > 2 * asmaxs) {
            peakTime <- sum(abs(as) > 0.5 * (asmax + asmaxs)) * (t[2] - t[1])
            labelShip <- sprintf(
                "%.1fg w/ spike to %.0fg for %.1fms (ship)",
                asmaxs / g, peakTime * 1e3, asmax / g
            )
        } else {
            labelShip <- sprintf("%.1fg (ship)", asmax / g)
        }
        aw <- derivative(vw, t)
        awmax <- max(abs(aw))
        awmaxs <- max(abs(runmed(aw, k)))
        if (awmax > 2 * awmaxs) {
            peakTime <- sum(abs(aw) > 0.5 * (awmax + awmaxs)) * (t[2] - t[1])
            labelWhale <- sprintf(
                "%.1fg w/ spike to %.0fg for %.1fms (whale)",
                awmaxs / g, peakTime * 1e3, awmax / g
            )
        } else {
            labelWhale <- sprintf("%.1fg (whale)", awmax / g)
        }
        mtext(paste("Max: ", labelWhale, ",", labelShip), side = 3, line = 0, cex = 0.8 * par("cex"))
        showEvents(xs, xw)
    }
    if (all || "section" %in% which) {
        skin <- x$WCF$compressed[, 1]
        blubber <- x$WCF$compressed[, 2]
        sublayer <- x$WCF$compressed[, 3]
        bone <- x$WCF$compressed[, 4]
        maxy <- max(c(skin + blubber + sublayer + bone))
        ylim <- c(-maxy * 1.2, 0)
        plot(t, -(skin + blubber + sublayer + bone),
            xlab = "Time [s]", ylab = "Whale-centred location [m]",
            type = "l", lwd = lwd, ylim = ylim, xaxs = "i", yaxs = "i", col = colwskin
        ) # outside skin margin
        lines(t, -(blubber + sublayer + bone), lwd = lwd, col = colwskin)
        lines(t, -(sublayer + bone), lwd = lwd, col = colwinterface) # , lty="42")
        lines(t, -bone, lwd = lwd, col = colwinterface) # , lty="42")
        showEvents(xs, xw)
        xusr <- par("usr")[1:2]
        x0 <- xusr[1] - 0.01 * (xusr[2] - xusr[1]) # snuggle up to axis
        text(x0, -0.5 * x$parms$l[4], "bone", pos = 4)
        text(x0, -x$parms$l[4] - 0.5 * x$parms$l[3], "sublayer", pos = 4)
        text(x0, -x$parms$l[4] - x$parms$l[3] - 0.5 * x$parms$l[2], "blubber", pos = 4)
        if (x$parms$l[1] > 0.1 * sum(x$parms$l)) {
            text(x0, -x$parms$l[4] - x$parms$l[3] - x$parms$l[2] - 0.5 * x$parms$l[1], "skin", pos = 4)
        }
        text(x0, 0.5 * (ylim[1] - x$parms$lsum), "", pos = 4)
    }
    if (all || "threat" %in% which) {
        # tcol <- rep(1, 4)
        skinzThreat <- x$WSF$sigmaz / x$parms$s[1]
        skinyThreat <- x$WSF$sigmay / x$parms$s[1]
        skinThreat <- ifelse(skinyThreat > skinzThreat, skinyThreat, skinzThreat)
        blubberThreat <- x$WCF$stress / x$parms$s[2]
        sublayerThreat <- x$WCF$stress / x$parms$s[3]
        boneThreat <- x$WCF$stress / x$parms$s[4]
        worst <- max(c(skinThreat, blubberThreat, sublayerThreat, boneThreat))
        trimThreat <- 10
        trimmed <- worst > trimThreat
        if (trimmed) {
            skinThreat <- pin(skinThreat, upper = trimThreat)
            blubberThreat <- pin(blubberThreat, upper = trimThreat)
            sublayerThreat <- pin(sublayerThreat, upper = trimThreat)
            boneThreat <- pin(boneThreat, upper = trimThreat)
            worst <- pin(worst, upper = trimThreat)
        }
        dy <- round(0.5 + worst)
        ylim <- c(0, 3 * dy + worst)
        plot(range(t), ylim, type = "n", xlab = "Time [s]", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
        axis(1)
        box()
        yTicks <- pretty(c(0, worst))
        mtext(
            paste(
                "Threat, Stress/Strength",
                if (trimmed) paste(" trimmed to", trimThreat) else ""
            ),
            side = 2, line = 2, cex = par("cex")
        )
        # Skin
        mtext("Skin", side = 4, at = 0, cex = 0.8 * par("cex"))
        y0 <- 0 # if (log) -1 else 0
        Y <- skinThreat
        fillplot4(t, skinThreat, yOffset = dy, breaks = c(1 / 4, 1 / 2, 3 / 4), col = colThreat)
        abline(h = 0)
        axis(2, at = y0 + yTicks, labels = yTicks)
        # Blubber
        mtext("Blubber", side = 4, at = dy, cex = 0.8 * par("cex"))
        Y <- x$WCF$stress / x$parms$s[2]
        tau25scaled <- x$parms$logistic$tau25 / x$parms$s[2]
        tau50scaled <- x$parms$logistic$tau50 / x$parms$s[2]
        tau75scaled <- x$parms$logistic$tau75 / x$parms$s[2]
        fillplot4(t, Y, yOffset = dy, breaks = c(tau25scaled, tau50scaled, tau75scaled), col = colThreat)
        abline(h = dy)
        axis(2, at = y0 + dy + yTicks, labels = rep("", length(yTicks)), tcl = 0.5)
        # Sublayer
        Y <- x$WCF$stress / x$parms$s[3]
        tau25scaled <- x$parms$logistic$tau25 / x$parms$s[3]
        tau50scaled <- x$parms$logistic$tau50 / x$parms$s[3]
        tau75scaled <- x$parms$logistic$tau75 / x$parms$s[3]
        mtext("Sublayer", side = 4, at = 2 * dy, cex = 0.8 * par("cex"))
        fillplot4(t, Y, yOffset = 2 * dy, breaks = c(tau25scaled, tau50scaled, tau75scaled), col = colThreat)
        abline(h = 2 * dy)
        axis(2, at = y0 + 2 * dy + yTicks, labels = yTicks)
        # Bone
        mtext("Bone", side = 4, at = 3 * dy, cex = 0.8 * par("cex"))
        fillplot4(t, boneThreat, yOffset = dy, breaks = c(1 / 4, 1 / 2, 3 / 4), col = colThreat)
        axis(2, at = y0 + 3 * dy + yTicks, labels = rep("", length(yTicks)), tcl = 0.5)
        abline(h = 3 * dy)
        showEvents(xs, xw)
    }
    if (all || "whale acceleration" %in% which) {
        a <- derivative(vw, t)
        plot(t, a, xlab = "Time [s]", ylab = "Whale accel. [m/s^2]", type = "l", lwd = lwd, xaxs = "i")
        mtext(sprintf("Max. %.3g m/s^2", max(a, na.rm = TRUE)), side = 3, line = 0, cex = 0.8 * par("cex"))
        showEvents(xs, xw)
    }
    if (all || "blubber thickness" %in% which) {
        WCF <- x$WCF
        y <- WCF$compressed[, 2]
        ylim <- c(min(0, min(y, na.rm = TRUE), na.rm = TRUE), max(y, na.rm = TRUE)) # include 0 if not there by autoscale
        plot(t, y, xlab = "Time [s]", ylab = "Blubber thickness [m]", type = "l", lwd = lwd, ylim = ylim, xaxs = "i")
        showEvents(xs, xw)
    }
    if (all || "sublayer thickness" %in% which) {
        WCF <- x$WCF
        y <- WCF$compressed[, 3]
        ylim <- c(min(0, min(y, na.rm = TRUE), na.rm = TRUE), max(y, na.rm = TRUE)) # include 0 if not there by autoscale
        plot(t, y, xlab = "Time [s]", ylab = "Sublayer thickness [m]", type = "l", lwd = lwd, ylim = ylim, xaxs = "i")
        showEvents(xs, xw)
    }
    if (all || "whale water force" %in% which) {
        y <- whaleWaterForce(vw, x$parms) / 1e6
        plot(t, y, xlab = "Time [s]", ylab = "Water force [MN]", type = "l", lwd = lwd, xaxs = "i")
        mtext(sprintf("Max. %.3g MN", max(y, na.rm = TRUE)), side = 3, line = 0, cex = 0.8 * par("cex"))
        showEvents(xs, xw)
    }
    if (all || "reactive forces" %in% which) {
        SF <- x$WSF$force
        CF <- x$WCF$force
        WWF <- x$WWF
        ylim <- range(c(WWF, SF, CF), na.rm = TRUE) / 1e6
        plot(t, SF / 1e6, type = "l", xlab = "Time [s]", ylab = "Forces [MN]", lwd = lwd, ylim = ylim, xaxs = "i")
        lines(t, CF / 1e6, lty = "dotted", lwd = lwd)
        # nolint start T_and_F_symbol_linter
        mtext(expression(" " * F[E]), side = 3, line = -1.2, adj = 0, cex = par("cex"))
        # nolint end T_and_F_symbol_linter
        mtext(" (solid)", side = 3, line = -2.2, adj = 0, cex = par("cex"))
        # nolint start T_and_F_symbol_linter
        mtext(expression(F[C] * " "), side = 3, line = -1.2, adj = 1, cex = par("cex"))
        # nolint end T_and_F_symbol_linter
        mtext(" (dotted) ", side = 3, line = -2.2, adj = 1, cex = par("cex"))
        mtext(sprintf("Max. %.3g MN", max(c(SF, CF, WWF) / 1e6, na.rm = TRUE)), side = 3, line = 0, cex = 0.8 * par("cex"))
        showEvents(xs, xw)
    }
    if (all || "skin stress" %in% which) {
        Fs <- whaleSkinForce(xs, xw, x$parms)
        ylim <- range(c(Fs$sigmay, Fs$sigmaz) / 1e6)
        plot(t, Fs$sigmay / 1e6, type = "l", xlab = "Time [s]", ylab = "Skin Stress [MPa]", lwd = lwd, ylim = ylim, xaxs = "i")
        lines(t, Fs$sigmaz / 1e6, lty = "dotted", lwd = lwd)
        mtext(" horiz.", side = 3, line = -1.2, adj = 0, cex = par("cex"))
        mtext(" (solid)", side = 3, line = -2.2, adj = 0, cex = par("cex"))
        mtext("vert. ", side = 3, line = -1.2, adj = 1, cex = par("cex"))
        mtext("(dotted) ", side = 3, line = -2.2, adj = 1, cex = par("cex"))
        mtext(sprintf("Max. %.3g MPa", max(c(Fs$sigmay, Fs$sigmaz) / 1e6, na.rm = TRUE)),
            side = 3, line = 0, cex = 0.8 * par("cex")
        )
        showEvents(xs, xw)
    }
    if (all || "compression stress" %in% which) {
        force <- x$WCF$force
        stress <- force / (x$parms$Lz * x$parms$Ly)
        plot(t, stress / 1e6, type = "l", xlab = "Time [s]", ylab = "Compress. Stress [MPa]", lwd = lwd, xaxs = "i")
        mtext(sprintf("Max. %.3g MPa", max(stress / 1e6, na.rm = TRUE)), side = 3, line = 0, cex = 0.8 * par("cex"))
        showEvents(xs, xw)
    }
    if (all || "lethality index" %in% which) {
        stress <- x$WCF$stress
        LI <- lethalityIndexFromStress(stress)
        plot(t, LI, type = "l", xlab = "Time [s]", ylab = "Lethality Index", lwd = lwd, xaxs = "i", ylim = c(0, 1), yaxs = "i")
        nt <- length(t)
        maxLI <- max(LI, na.rm = TRUE)
        # Redraw the supercritical in a thicker line. But, first, refine the grid if it's coarse.
        dangerTime <- t[LI > 0.5]
        dangerTimeInterval <- 0.0
        if (length(dangerTime) > 0L) {
            dangerTimeInterval <- diff(range(dangerTime, na.rm = TRUE))
            if (nt < 2000) {
                t2 <- seq(t[1], t[nt], length.out = 2000)
                LI2 <- approx(t, LI, t2)$y
                LI2[LI2 < 0.5] <- NA
                if (any(is.finite(LI2))) {
                    lines(t2, LI2, lwd = 2 * lwd)
                }
            } else {
                LI[LI < 0.5] <- NA
                if (any(is.finite(LI))) {
                    lines(t, LI, lwd = 2 * lwd)
                }
            }
        }
        abline(h = 0.5, lty = "dotted")
        if (dangerTimeInterval == 0.0) {
            mtext(sprintf("Max. %.2g", maxLI), side = 3, line = 0, cex = 0.8 * par("cex"))
        } else {
            mtext(sprintf("Max. %.2g (>0.5 for %.1gs)", maxLI, dangerTimeInterval), side = 3, line = 0, cex = 0.8 * par("cex"))
        }
        showEvents(xs, xw)
    }
    if (all || "values" %in% which) {
        omar <- par("mar")
        par(mar = rep(0, 4))
        parms <- x$parms[unlist(lapply(x$parms, function(p) is.vector(p)))]
        parms$logistic <- NULL # we have no GUI for this, so do not display
        parms$engineForce <- NULL # inserted during calculation, not user-supplied
        parms$lsum <- NULL # inserted during calculation, not user-supplied
        parms <- lapply(parms, function(x) signif(x, 4))
        parms <- lapply(parms, function(p) deparse(p))
        parms$vs_knots <- mps2knot(x$vs[1])
        names <- names(parms)
        values <- unname(unlist(parms))
        n <- 1 + length(values)
        plot(1:n, 1:n, type = "n", xlab = "", ylab = "", axes = FALSE)
        o <- order(names(parms), decreasing = TRUE)
        for (i in seq_along(values)) {
            text(1, i + 0.5, paste(names[o[i]], "=", values[o[i]]), pos = 4, cex = 1)
        }
        par(mar = omar)
    }
}

