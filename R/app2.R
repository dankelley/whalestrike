# vim:textwidth=100:expandtab:shiftwidth=4:softtabstop=4

help <- "
<p>
When this app starts, the sliders and tick-boxes are set up to model a small
fishing boat, of mass 45 tonnes, moving at 10 knots towards a whale of
length 13.7m. (The whale length is used to compute its mass, using a formula
that is described by the output of typing
<tt>help(\"whaleMassFromLength\",\"whalestrike\")</tt>
in an R console).
</p>

<p>
Sliders are provided for setting certain key properties of the ship and the
whale. The details of these and the other parameters are revealed by
typing <tt>help(\"parameters\",\"whalestrike\")</tt>
and
<tt>help(\"strike\",\"whalestrike\")</tt>
in an R console.
</p>

Checkboxes are provided to select the desired plot types.
For the details of these type, type
<tt>help(\"plot.strike\",\"plot\")</tt>
in an
R console.)

Typing
<tt>\"?app2\"</tt>
will yield some more information about this app."

#' GUI application for whale simulation (version 2)
#'
#' This is similar to [app()], except that it relies on the `bslib` package
#' to provide a cleaner interface, in which sub-windows for controllers
#' can be opened and closed by the user. Since `bslib` is still experimental,
#' [app2()] may not work in the future.
#'
#' Compared with [app()], the present function lacks the ability to save settings
#' and reload them later. This is mainly because it only works with locally-run
#' operations, not from server-run operations.  The latter would require extra
#' coding to set up user's storage space, to prevent against web attacks, etc.,
#' which is beyond the present purpose.  However, there is an addition with
#' [app2()] that might prove more useful: a button to display the code required
#' to reproduce the simulated state.  This may be of help to the those seeking
#' to explore the simulator in programmatically, for precision and reproducibility.
#'
#' Sliders, buttons, and choosers are grouped into panes that appear on
#' the left of the view. When [app2()] first opens, the only pane that is
#' open is the "ship" pane.  To start, try adjusting
#' the "Ship Speed" slider in that pane, to see what happens with
#' the plots.  In particular, notice how the "Lethality Index" plot indicates
#' more threatening conditions, as ship speed is increased.  Then, at a fixed
#' speed, try lowering the "Impact width" controller.  You will see that
#' the lethality index. If that does not make sense, think of the difference
#' about the penetration depth for a hammer hitting a piece of wood, versus
#' a hammer hitting a nail held to that wood.  Proceed in this way, adjusting
#' controllers and trying to build intuition.  As you explore the controllers
#' in this way, you may find it helpful to examine the code produced by the
#' "Code" button.  Note that the package provides documentation on all
#' the functions that are used in that code.
#'
#' The default setup, which is shown before the user alters any of the sliders,
#' etc., is a simulation of a small fishing boat, of mass 45 tonnes, moving at
#' speed 10 knots towards a whale of length 13.7m. (The whale length is used to
#' compute its mass, using a formula that is described by the output of typing
#' `help("whaleMassFromLength","whalestrike")` in an R console).
#'
#' Much can be learned by adjusting the sliders and examining the plotted output.
#' As an exercise, try setting to a particular ship mass of interest, and then to
#' slide the ship speed to higher and lower values, whilst monitoring the "threat"
#' panel for black regions. This will reveal a critical speed for conditions that
#' threaten the whale.  Next, try altering the sublayer thickness, which is a
#' surrogate for location along the whale body, because e.g. the sublayer is
#' thinner near the mandible.
#'
#' Advanced users are likely to want to alter the values of impact width and
#' height. The default setting are intended to mimic a small fishing boat, such as
#' a Cape Islander. Try lowering the width, to simulate a strike by a daggerboard
#' or keel of a sailing boat.
#'
#' Note that the pulldown menu for setting the whale species affects *only* whale
#' area and mass. It does not affect the thicknesses of blubber or sublayer, the values of
#' which have been set up to represent a midsection strike on a North Atlantic
#' Right Whale.
#'
#' @param debug logical value indicating whether to print output to
#' the R console as the computation is done.
#'
#' @importFrom bslib accordion accordion_panel card sidebar
#' @importFrom shiny actionButton checkboxGroupInput h6 observeEvent plotOutput selectInput shinyApp showModal sliderInput stopApp
#' @export
#'
#' @family interactive apps
#'
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
        sidebar = bslib::sidebar(
            # title = "Controls", # a waste of space (who could not figure this out?)
            width = 280, # default, 250, too narrow for one of the plot types
            bslib::accordion(
                open = "ship",
                multiple = TRUE,
                bslib::accordion_panel(
                    "inspect",
                    shiny::actionButton("help", "Help"),
                    shiny::actionButton("quit", "Quit"),
                    shiny::actionButton("code", "Code"),
                    shiny::sliderInput("tmax", shiny::h6("Max time [s]"),
                        ticks = FALSE,
                        min = 0.1, max = 5, value = 1, step = 0.05
                    )
                ),
                accordion_panel(
                    "plot",
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
                    shiny::sliderInput("ms", shiny::h6("Mass [tonne]"),
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
                    ),
                    shiny::sliderInput("vs", shiny::h6("Speed [knot]"),
                        ticks = FALSE,
                        min = 1, max = 30, value = 10, step = 0.1
                    )
                )
            )
        ),
        bslib::card(
            shiny::plotOutput("plot")
        )
    )
    server <- function(input, output, session) {
        dmsg("in server")
        whaleMass <- function(length, species)
        {
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
        shiny::observeEvent(input$help, {
            shiny::showModal(shiny::modalDialog(shiny::HTML(help), title = "Using this application", size = "l"))
        })
        observeEvent(input$quit, {
            shiny::stopApp()
        })
        shiny::observeEvent(input$code, {
            msg <- "<pre>library(whalestrike)<br>"
            msg <- paste0(msg,
                "# Simulation time interval<br>",
                "t <- seq(0.0, ", input$tmax, ", length.out = 2000)<br>")
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
                "    mw = whaleMassFromLength(length = ", input$lw, ", species = ", input$species, "),<br>",
                "    Sw = whaleAreaFromLength(", input$lw, ",<br>",
                "        species = \"", input$species, " Whale\")<br>",
                "    l = c(", input$l1 / 100,
                ", ", input$l2 / 100,
                ", ", input$l3 / 100,
                ", ", input$l4 / 100, "))<br>"
            )
            msg <- paste0(msg,
                "# Perform the simulation<br>")
            msg <- paste0(msg, "sol <- strike(t, state, parms)<br>")
            if (length(input$plot_panels)) {
                msg <- paste0(msg,
                    "# Plot results individually<br>")
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
                    showNotification("Auto-refined grid to capture acceleration peak")
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