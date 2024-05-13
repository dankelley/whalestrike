# vim:textwidth=100:expandtab:shiftwidth=4:softtabstop=4

help <- "
<p>

When this app starts, the sliders and tick-boxes are set up to model a
small (45 tonne) fishing boat moving at 10 knots towards a whale of
length 13.7m.

</p>

<p>

Clicking on the down-arrow icons labelled <tt>plot</tt>,
<tt>whale</tt> and <tt>ship</tt> opens up controllers for setting
certain key properties of the simulation. See the package
documentation to learn more about these things.

</p>

<p>

More information on this app may be retrieved by typing <tt>?app2</tt>
in the R console, and a video demonstration is provided on
<a href=\"https://youtu.be/kTMl3nXa5A4\">youtube</a>).

</p>

<p>

The simulation works by calling the <tt>strike()</tt> function within
the package, so a good way to learn more about it is to type
<tt>?strike</tt> in the R console. Also, please note that clicking the
</ttCode</tt> button, which displays the code used to set up the
simulation under view at any given time.

</p>

"

#' GUI application for whale simulation (version 2)
#'
#' This is similar to [app()], except that it relies on the `bslib` package
#' to provide a cleaner interface, in which sub-windows for controllers
#' can be opened and closed by the user. Since `bslib` is still experimental,
#' [app2()] may start failing in the future, and users are asked to
#' report an issue on the package github page, so the author
#' can make adjustments.
#'
#' Compared with [app()], the present function lacks the ability to save settings
#' and reload them later. This is mainly because it only works with locally-run
#' operations, not from server-run operations.  The latter would require extra
#' coding to set up user's storage space, to prevent against web attacks, etc.,
#' which is beyond the present purpose.  However, there is an addition with
#' [app2()] that might prove more useful: a button to display the code required
#' to reproduce the simulated state.  This may be of help to the those seeking
#' to explore the results of simulations more precisely and with
#' greater reproducibility.
#'
#' Sliders, buttons, and choosers are grouped into panes that appear on
#' the left of the view. When [app2()] first opens, all of these panes
#' are closed. To get acquainted with the app, try adjusting the controllers
#' that *are* visible on the initial view.  Then, open the "ship" pane and
#' increase the ship mass.  Do you find that the results make qualitative
#' sense?  Continue this process, exploring all the panes. It is hoped
#' that a half hour of such exploration will let users start to
#' investigate practical applications.  For more about how the simulations
#' are carried out, as well as comments on some applications that
#' may be of interest, please consult Kelley et al. (2021).
#'
#' More information on [app2()] in video form on
#' [youtube](https://youtu.be/kTMl3nXa5A4).
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
#' Science, 37(1), 2021. \doi{10.1111/mms.12745}.
##'
#' @author Dan Kelley
app2 <- function(debug = FALSE) {
    dmsg <- function(...) {
        if (debug) message(...)
    }
    if (!requireNamespace("bslib")) {
        stop("must install.packages(\"bslib\") for app2() to work")
    }
    if (!requireNamespace("shiny")) {
        stop("must install.packages(\"shiny\") for app2() to work")
    }
    ui <- bslib::page_sidebar(
        # title = "app2", # a waste of space (the user launched this, and can do help)
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
                            "The default values are meant to represent",
                            "an adult Right Whale."
                        )
                    ),
                    shiny::selectInput("species", "Species",
                        choices = c(
                            "N. Atl. Right", "Blue", "Bryde", "Fin", "Gray", "Humpback", "Minke",
                            "Pac. Right", "Sei", "Sperm"
                        ),
                        selected = "N. Atl. Right"
                    ),
                    shiny::sliderInput("lw", shiny::h6("Length [m]"),
                        ticks = FALSE,
                        min = 5, max = 20, value = 13.7, step = 0.1
                    ),
                    shiny::sliderInput("l3", shiny::h6("Sublayer thickness [cm]"),
                        ticks = FALSE,
                        min = 5, max = 200, value = 112, step = 1
                    ),
                    shiny::sliderInput("l4", shiny::h6("Bone thickness [cm]"),
                        ticks = FALSE,
                        min = 5, max = 30, value = 10, step = 1
                    ),
                    shiny::sliderInput("l1", shiny::h6("Skin thickness [cm]"),
                        ticks = FALSE,
                        min = 1, max = 3, value = 2.5, step = 0.1
                    ),
                    shiny::sliderInput("l2", shiny::h6("Blubber thickness [cm]"),
                        ticks = FALSE,
                        min = 5, max = 40, value = 16, step = 1
                    )
                ),
                accordion_panel(
                    "ship",
                    tooltip(
                        "\u24D8",
                        paste(
                            "These sliders set the mass of the ship and the",
                            "geometry of the impact zone. For the latter,",
                            "consider the area a few decimetres aft of the first",
                            "contact point. The default values correspond to a small",
                            "fishing vessel with a prow shape similar to that of a",
                            "Cape Islander."
                        )
                    ),
                    shiny::sliderInput("ms",
                        shiny::h6("Mass [tonne]"),
                        ticks = FALSE,
                        min = 10, max = 500, value = 45, step = 1
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
                plotOutput("plot")
            ),
        ),
    )
    server <- function(input, output, session) {
        dmsg("in server")
        whaleMass <- function(length, species) {
            if (species == "N. Atl. Right") {
                whaleMassFromLength(length, species = "N. Atl. Right Whale", model = "fortune2012")
            } else if (species == "Blue") {
                whaleMassFromLength(length, species = "Blue Whale", model = "lockyer1976")
            } else if (species == "Bryde") {
                whaleMassFromLength(length, species = "Bryde Whale", model = "lockyer1976")
            } else if (species == "Fin") {
                whaleMassFromLength(length, species = "Fin Whale", model = "lockyer1976")
            } else if (species == "Gray") {
                whaleMassFromLength(length, species = "Gray Whale", model = "lockyer1976")
            } else if (species == "Humpback") {
                whaleMassFromLength(length, species = "Humpback Whale", model = "lockyer1976")
            } else if (species == "Minke") {
                whaleMassFromLength(length, species = "Minke Whale", model = "lockyer1976")
            } else if (species == "Pac. Right") {
                whaleMassFromLength(length, species = "Pac. Right Whale", model = "lockyer1976")
            } else if (species == "Sei") {
                whaleMassFromLength(length, species = "Sei Whale", model = "lockyer1976")
            } else if (species == "Sperm") {
                whaleMassFromLength(length, species = "Sperm Whale", model = "lockyer1976")
            } else {
                stop("programming error: species not handled: ", species)
            }
        }
        shiny::observeEvent(input$keypressTrigger, {
            key <- intToUtf8(input$keypress)
            # NOTE: this keystroke is not explained. I may delete
            # it. And I might add other keystrokes. One that I
            # think might be good would be to try small-ship and
            # large-ship simulations.
            if (key == "?") {
                shiny::showModal(shiny::modalDialog(shiny::HTML(help), title = "Using this application", size = "l"))
            }
        })
        shiny::observeEvent(input$help, {
            shiny::showModal(shiny::modalDialog(shiny::HTML(help), title = "Using this application", size = "l"))
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
            msg <- paste0(
                msg,
                "# Initial state<br>",
                "state <- list(<br>",
                "    xs = ", -(1 + 0.01 * (input$l1 + input$l2 + input$l3 + input$l4)), ",<br>",
                "    vs = knot2mps(", input$vs, "),<br>",
                "    xw = 0, vw = 0)<br>"
            )
            msg <- paste0(
                msg,
                "# Whale parameters<br>",
                "parms <- parameters(<br>",
                "    ms = ", 1000 * input$ms, ",<br",
                "    Ss = shipAreaFromMass(1000 * ", input$ms, "),<br>",
                "    Ly = ", input$Ly, ",<br>",
                "    Lz = ", input$Lz, ",<br>",
                "    lw = ", input$lw, ",<br>",
                "    mw = whaleMassFromLength(", input$lw, ", \"", input$species, " Whale\"),<br>",
                "    Sw = whaleAreaFromLength(", input$lw, ", \"", input$species, " Whale\"),<br>",
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
                mw <- whaleMass(input$lw, input$species)
                dmsg("mw=", mw)
                parms <- whalestrike::parameters(
                    ms = 1000 * input$ms,
                    Ss = shipAreaFromMass(1000 * input$ms),
                    Ly = input$Ly,
                    Lz = input$Lz,
                    lw = input$lw,
                    mw = mw,
                    Sw = whalestrike::whaleAreaFromLength(input$lw, species = "N. Atl. Right Whale", "wetted"),
                    l = c(input$l1 / 100, input$l2 / 100, input$l3 / 100, input$l4 / 100)
                )
                if (debug) print(parms)
                state <- list(xs = -(1 + parms$lsum), vs = whalestrike::knot2mps(input$vs), xw = 0, vw = 0)
                if (debug) print(state)
                t <- seq(0, input$tmax, length.out = 2000)
                if (debug) print(range(t))
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
