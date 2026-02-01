# vim:textwidth=100:expandtab:shiftwidth=4:softtabstop=4

helpText <- "<p>When this app starts, the sliders and tick-boxes are set
up to model a small (45 tonne) fishing boat moving at 10 knots towards a
North Atlantic Right Whale of length 13.7m.</p>

<p>Clicking on the down-arrow icons labelled <tt>plot</tt>,
<tt>whale</tt> and <tt>ship</tt> opens up controllers for setting
certain key properties of the simulation. Hovering the mouse
over the circled 'i' symbols provides more information on
how to work with these controllers.</p>

<p>The simulation works by calling the <tt>strike()</tt> function within
the package, so a good way to learn more about it is to type
<tt>?strike</tt> in the R console. Also, please note that clicking the
</tt>Code</tt> button will show a pop-up box that displays code that
will mimic the simulation being viewed.</p>

<p>A video demonstration of an older version of this app
is at <a href=\"https://youtu.be/kTMl3nXa5A4\">youtube</a>).
You can learn more about this app, and about the whalestrike
package, by consulting the package documentation.</p>"

#' GUI application for whale simulation
#'
#' Graphical-user-interface tool for exploring whale-strike simulations.
#'
#' Sliders, buttons, and choosers are grouped into panes that appear on
#' the left of the view. When [app()] first opens, all of these panes
#' are closed. To get acquainted with the app, try adjusting the controllers
#' that *are* visible on the initial view.  Then, open the "ship" pane and
#' increase the ship mass.  Do you find that the results make qualitative
#' sense?  Continue this process, exploring all the panes. A
#' half-hour of such exploration should be enough to build enough
#' confidence to start investigating practical applications.
#' To learn more about how the simulations are carried out, and
#' to read more about the underlying goals of this tool,
#' please consult Kelley et al. (2021) and Kelley (2024). Extensive
#' details on the calculations are provided in the help pages
#' for the various functions of the whalestrike package, of
#' which that for [whalestrike()] is a good starting point.
#'
#' More information on [app()] in video form on
#' [youtube](https://youtu.be/kTMl3nXa5A4).
#'
#' Note that an older version of a similar GUI application is
#' still available as [app_2025()], but it is not maintained
#' and is slated for removal in the early months of 2026.
#'
#' @param debug logical value indicating whether to print output to
#' the R console as the computation is done.
#'
#' @importFrom bslib accordion accordion_panel card card_header sidebar tooltip
#'
#' @importFrom shiny actionButton checkboxGroupInput h6 observeEvent plotOutput selectInput shinyApp showModal sliderInput stopApp
#'
#' @export
#'
#' @family interactive apps
#'
#' @references
#'
#' Kelley, Dan E., James P. Vlasic, and Sean W. Brillant. "Assessing the Lethality
#' of Ship Strikes on Whales Using Simple Biophysical Models." Marine Mammal
#' Science 37, no. 1 (2021): 251–67. https://doi.org/10.1111/mms.12745.
#'
#' Kelley, Dan E."“Whalestrike: An R Package for Simulating Ship Strikes on Whales."
#' Journal of Open Source Software 9, no. 97 (2024): 6473.
#' https://doi.org/10.21105/joss.06473.
#'
#' Mayette, Alexandra. "Whale Layer Thickness." December 15, 2025. (Personal
#' communication of a 5-page document.)
#'
#' Mayette, Alexandra, and Sean W. Brillant. "A Regression-Based Method to Estimate
#' Vessel Mass for Use in Whale-Ship Strike Risk Models." PloS One 21, no. 1 (2026):
#' e0339760. https://doi.org/10.1371/journal.pone.0339760.
#'
#' @author Dan Kelley
app <- function(debug = FALSE) {
    dmsg <- function(...) {
        if (debug) message(...)
    }
    if (!requireNamespace("bslib")) {
        stop("must install.packages(\"bslib\") for app() to work")
    }
    if (!requireNamespace("shiny")) {
        stop("must install.packages(\"shiny\") for app() to work")
    }
    ui <- bslib::page_sidebar(
        # title = "app", # a waste of space (the user launched this, and can do help)
        shiny::tags$script(paste0(
            "$(document).on(\"keypress\", function (e) {",
            "Shiny.onInputChange(\"keypress\", e.which);",
            "Shiny.onInputChange(\"keypressTrigger\", Math.random());",
            "});"
        )),
        sidebar = bslib::sidebar(
            # title = "Controls", # a waste of space (who could not figure this out?)
            tooltip(
                "\u24D8",
                paste(
                    "The three 'accordion' menus below are used to select the",
                    "desired plot types, and to set the characteristics of the",
                    "whale and the ship that collides with it."
                )
            ),
            width = 280, # default, 250, too narrow for one of the plot types
            bslib::accordion(
                open = FALSE,
                multiple = TRUE,
                accordion_panel(
                    "plot",
                    tooltip(
                        "\u24D8",
                        paste(
                            "The checkboxes select the plots to be displayed.",
                            "To learn more about the choices, try typing",
                            "'?whalestrike::plot.strike' in an R session."
                        )
                    ),
                    shiny::checkboxGroupInput("plot_panels", "",
                        choices = c(
                            "location",
                            "section",
                            "lethality index",
                            "threat",
                            "whale acceleration",
                            "whale water force",
                            "reactive forces",
                            "skin stress",
                            "compression stress",
                            "values"
                        ),
                        selected = c("location", "section", "lethality index")
                    )
                ),
                accordion_panel(
                    "whale",
                    tooltip(
                        "\u24D8",
                        paste(
                            "These sliders set the whale characteristics.",
                            "Select 'Generic' for the whale model discussed in Kelley et al. (2021)",
                            "and Kelley (2024). Otherwise, select one of the provided list of species.",
                            "In this latter case, the sliders are adjusted according to a tabulation",
                            "of whale properties assembled by Alexandra Mayette."
                        )
                    ),
                    shiny::selectInput("species", "Species",
                        choices = c(
                            paste("Generic"),
                            paste("Blue", "Whale"),
                            paste("Bryde", "Whale"),
                            paste("Fin", "Whale"),
                            paste("Gray", "Whale"),
                            paste("Humpback", "Whale"),
                            paste("Minke", "Whale"),
                            paste("N.", "Atl.", "Right", "Whale"),
                            paste("Pac.", "Right", "Whale"),
                            paste("Sei", "Whale"),
                            paste("Sperm", "Whale")
                        ),
                        selected = "Generic"
                    ),
                    shiny::sliderInput("lw", shiny::h6("Length [m]"),
                        ticks = FALSE,
                        min = 5, max = 25, value = 13.7, step = 0.1
                    ),
                    # bone: range in whaleMeasurements() is 7.8 to 17.3 cm
                    shiny::sliderInput("l4", shiny::h6("Bone thickness [cm]"),
                        ticks = FALSE,
                        min = 1, max = 20, value = 10, step = 1
                    ),
                    # sublayer: range in whaleMeasurements() is 31.3 to 168.7 cm
                    shiny::sliderInput("l3", shiny::h6("Sublayer thickness [cm]"),
                        ticks = FALSE,
                        min = 20, max = 200, value = 112, step = 1
                    ),
                    # blubber: range in whaleMeasurements() is 3.3 to 16.3 cm
                    shiny::sliderInput("l2", shiny::h6("Blubber thickness [cm]"),
                        ticks = FALSE,
                        min = 1, max = 40, value = 16, step = 1
                    ),
                    # skin: range in whaleMeasurements() is 0.2 to 1.0 cm
                    shiny::sliderInput("l1", shiny::h6("Skin thickness [cm]"),
                        ticks = FALSE,
                        min = 0.0, max = 3, value = 2.5, step = 0.1
                    )
                ),
                accordion_panel(
                    "ship",
                    tooltip(
                        "\u24D8",
                        paste(
                            "These sliders set the ship characteristics. There are two choices.",
                            "Choise 1: select 'Generic' for a vessel model discussed in Kelley et al. (2021)",
                            "and Kelley (2024). Choice 2: select one of the specified vessel",
                            "types, and then use a slider to set vessel length, after which the app",
                            "computes ship mass using shipMassFromLength()."
                        )
                    ),
                    shiny::selectInput("vessel", "Vessel",
                        choices = c(
                            "Generic",
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
                    ),
                    shiny::conditionalPanel(
                        condition = "input.vessel == 'Generic'",
                        shiny::sliderInput("ms",
                            shiny::h6("Mass [tonne]"),
                            ticks = FALSE,
                            min = 10, max = 500, value = 45, step = 1
                        ),
                    ),
                    shiny::conditionalPanel(
                        condition = "input.vessel != 'Generic'",
                        shiny::sliderInput("vesselLength",
                            shiny::h6("Length [m]"),
                            ticks = FALSE,
                            min = 10, max = 500, value = 50, step = 5
                        ),
                    ),
                    shiny::sliderInput("Ly", shiny::h6("Impact width [m]"),
                        ticks = FALSE,
                        min = 0.1, max = 2, value = 1.15, step = 0.05
                    ),
                    shiny::sliderInput("Lz", shiny::h6("Impact height [m]"),
                        ticks = FALSE,
                        min = 0.1, max = 2, value = 1.15, step = 0.05
                    )
                ),
            ),
            bslib::tooltip(
                shiny::sliderInput("tmax", shiny::h6("Max time [s]"),
                    ticks = FALSE,
                    min = 0.1, max = 5, value = 1, step = 0.05
                ),
                "Set the time interval to be simulated (in seconds)."
            ),
            bslib::tooltip(
                shiny::sliderInput("vs", shiny::h6("Ship Speed [knot]"),
                    ticks = FALSE,
                    min = 1, max = 30, value = 10, step = 0.1
                ),
                "Set the initial (pre-impact) ship speed, in knots. Note that 1 knot is 0.514 m/s."
            ),
            bslib::tooltip(
                shiny::actionButton("help", "Help"),
                "Open a pop-up window that provides some information on using this app."
            ),
            bslib::tooltip(
                shiny::actionButton("code", "Code"),
                paste(
                    "Open a pop-up window that provides R code that is aligned with",
                    "the present app settings. This can be useful for users who want to set up",
                    "a suite of simulations, e.g. to produce graphs of lethality index versus ship speed."
                )
            ),
            bslib::tooltip(
                shiny::actionButton("quit", "Quit"),
                "Quite this application."
            ),
        ),
        card(
            card_header(
                "", # Results of collision simulation",
                tooltip(
                    "\u24D8",
                    paste(
                        "This diagram represents the results of a simulation",
                        "of a whale being struck by a ship, with parameters",
                        "set by the controllers shown to the left."
                    )
                ),
                plotOutput("plot", width="100%", height="600px")
            ),
        ),
    )
    server <- function(input, output, session) {
        dmsg("in server")
        whaleMass <- function(length, species) {
            if (species == "Generic") {
                whaleMassFromLength(length, species = "N. Atl. Right Whale", model = "fortune2012")
            } else if (species == "N. Atl. Right Whale") {
                whaleMassFromLength(length, species = species, model = "fortune2012")
            } else if (species == "Blue Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Bryde Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Fin Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Gray Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Humpback Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Minke Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Pac. Right Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Sei Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else if (species == "Sperm Whale") {
                whaleMassFromLength(length, species = species, model = "lockyer1976")
            } else {
                stop("programming error: species not handled: ", species)
            }
        }
        shiny::observeEvent(input$vessel, {
            if (input$vessel != "Generic") {
                length <- shipLength(input$vessel)
                shiny::updateSliderInput(session, "vesselLength", shiny::h6("Length [m]"), length)
            }
        })
        shiny::observeEvent(input$species, {
            wm <- if (identical(input$species, "Generic")) {
                p <- parameters()
                with(p, list(length = lw, skin = l[1], blubber = l[2], sublayer = l[3], bone = l[4]))
            } else {
                whaleMeasurements(input$species)
            }
            shiny::updateSliderInput(session, "lw", shiny::h6("Length [m]"), wm$length)
            shiny::updateSliderInput(session, "l4", shiny::h6("Bone thickness [cm]"), 100 * wm$bone)
            shiny::updateSliderInput(session, "l3", shiny::h6("Sublayer thickness [cm]"), 100 * wm$sublayer)
            shiny::updateSliderInput(session, "l2", shiny::h6("Blubber thickness [cm]"), 100 * wm$blubber)
            shiny::updateSliderInput(session, "l1", shiny::h6("skin thickness [cm]"), 100 * wm$skin)
        })
        #<unused> shiny::observeEvent(input$keypressTrigger, {
        #<unused>     key <- intToUtf8(input$keypress)
        #<unused>     # NOTE: this keystroke is not explained. I may delete
        #<unused>     # it. And I might add other keystrokes. One that I
        #<unused>     # think might be good would be to try small-ship and
        #<unused>     # large-ship simulations.
        #<unused>     if (key == "?") {
        #<unused>         shiny::showModal(shiny::modalDialog(shiny::HTML(help), title = "Using this application", size = "xl"))
        #<unused>     }
        #<unused> })
        shiny::observeEvent(input$help, {
            shiny::showModal(shiny::modalDialog(shiny::HTML(helpText),
                title = "Using this application", size = "l"
            ))
            # print(help)
            # shiny::showModal(shiny::modalDialog(shiny::HTML("WTF"), title = "Using this application", size = "l"))
        })
        observeEvent(input$quit, {
            shiny::stopApp()
        })
        shiny::observeEvent(input$code, {
            msg <- "<pre>library(whalestrike)<br>"
            msg <- paste0(
                msg,
                "# Simulation time interval<br>",
                "t <- seq(0.0, ", input$tmax, ", length.out = 2000)<br>"
            )
            # BEGIN bookmark 1 (MUST keep in-synch with bookmark 2 below)
            msg <- paste0(
                msg,
                "# Initial state<br>",
                "state <- list(<br>",
                "    xs = ", -(1 + 0.01 * (input$l1 + input$l2 + input$l3 + input$l4)), ",<br>",
                "    vs = knot2mps(", input$vs, "),<br>",
                "    xw = 0, vw = 0)<br>"
            )
            ms <- if (input$vessel == "Generic") {
                1000 * input$ms
            } else {
                shipMassFromLength(input$vessel, input$vesselLength)
            }
            Ss <- shipAreaFromMass(ms)
            Sw <- whaleAreaFromLength(input$lw, species = "N. Atl. Right Whale", "wetted")
            mw <- whaleMass(input$lw, input$species)
            # END bookmark 1
            msg <- paste0(
                msg,
                "# Whale parameters<br>",
                "parms <- parameters(<br>",
                "    ms = ", round(ms), ",<br",
                "    Ss = ", round(Ss), ",<br>",
                "    Ly = ", input$Ly, ",<br>",
                "    Lz = ", input$Lz, ",<br>",
                "    lw = ", input$lw, ",<br>",
                "    mw = ", round(mw), ",<br>",
                "    Sw = ", round(Sw), ",<br>",
                "    l = c(", input$l1 / 100,
                ", ", input$l2 / 100,
                ", ", input$l3 / 100,
                ", ", input$l4 / 100, "))<br>"
            )
            msg <- paste0(
                msg,
                "# Perform the simulation<br>"
            )
            msg <- paste0(msg, "sol <- strike(t, state, parms)<br>")
            if (length(input$plot_panels)) {
                msg <- paste0(
                    msg,
                    "# Plot results<br>"
                )
                npanels <- length(input$plot_panels)
                nrows <- floor(sqrt(npanels))
                ncols <- ceiling(npanels / nrows)
                msg <- paste0(msg, "par(mfrow = c(", nrows, ", ", ncols, "))<br>")
                msg <- paste0(msg, "par(mar = c(3.2, 3, 2.5, 2))<br>")
                msg <- paste0(msg, "par(mgp = c(1.7, 0.6, 0))<br>")
                for (which in input$plot_panels) {
                    msg <- paste0(msg, "plot(sol, which = \"", which, "\")<br>")
                }
            }
            msg <- paste0(msg, "</pre><br>")
            shiny::showModal(shiny::modalDialog(shiny::HTML(msg), size = "l"))
        })
        output$plot <- renderPlot(
            {
                # BEGIN bookmark 2 (MUST keep in-synch with bookmark 1 above)
                ms <- if (input$vessel == "Generic") {
                    1000 * input$ms
                } else {
                    shipMassFromLength(input$vessel, input$vesselLength)
                }
                Ss <- shipAreaFromMass(ms)
                Sw <- whaleAreaFromLength(input$lw, species = "N. Atl. Right Whale", "wetted")
                mw <- whaleMass(input$lw, input$species)
                # END bookmark 2
                l <- c(input$l1 / 100, input$l2 / 100, input$l3 / 100, input$l4 / 100)
                if (debug) {
                    message("About to run strike with params:")
                    message("  ms=", round(ms), " [kg] ship mass")
                    message("  Ss=", round(Ss), " [m^2] ship area")
                    message("  Ly=", input$Ly, " [m] impact area width")
                    message("  Lz=", input$Lz, " [m] impact area height")
                    message("  mw=", round(mw), " [kg] whale mass")
                    message("  Sw=", round(Sw), " [m^2] whale area")
                    message("  l=", paste(l, collapse=","), " [m] whale layer thicknesses")
                }
                parms <- whalestrike::parameters(
                    ms = ms,
                    Ss = Ss,
                    Ly = input$Ly,
                    Lz = input$Lz,
                    lw = input$lw,
                    mw = mw,
                    Sw = Sw,
                    l = l
                )
                state <- list(xs = -(1 + parms$lsum), vs = whalestrike::knot2mps(input$vs), xw = 0, vw = 0)
                if (debug) {
                    message("and with state:")
                    message("  xs=", state$xs, " [m] ship position")
                    message("  vs=", round(state$vs,2), " [m/s] ship speed")
                    message("  xs=", state$xw, " [m] whale position")
                    message("  vs=", state$vw, " [m/s] whale speed\n")
                }
                t <- seq(0, input$tmax, length.out = 2000)
                sol <- strike(t, state, parms)
                if (sol$refinedGrid) {
                    showNotification("Refined grid for accel. peak")
                }
                npanels <- length(input$plot_panels)
                nrows <- floor(sqrt(npanels))
                ncols <- ceiling(npanels / nrows)
                par(mfrow = c(nrows, ncols), mar = c(3.2, 3, 2.5, 2), mgp = c(1.7, 0.6, 0), cex = 1)
                for (which in input$plot_panels) {
                    plot(sol, which = which)
                }
            },
            pointsize = 18
        )
    }
    shiny::shinyApp(ui = ui, server = server)
}
