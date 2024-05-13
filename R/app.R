# vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4


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
whale, with italic labels used for properties likely to be adjusted during
simulations.  The details of these and the other parameters are revealed by
typing <tt>help(\"parameters\",\"whalestrike\")</tt>
and
<tt>help(\"strike\",\"whalestrike\")</tt>
in an R console.
</p>

To the right of the sliders is a column of checkboxes that control the plotted
output. (For the details of the plots, type
<tt>help(\"plot.strike\",\"plot\")</tt>
in an
R console.) At startup, three of these boxes are ticked, yielding a display
with three panels showing the time history of the simulation.  The left-hand
plot panel shows whale and boat location, the former with an indication of the
interfaces between skin, blubber, sublayer, and bone. Peak ship and whale
accelerations are indicated with labels inside this panel.  The middle panel
shows the same information as the one to its left, but with a whale-centred
coordinate system, and with labels for the components. The right panel is an
indication of the estimated threat to the four layers of the whale, with curves
that are filled with colours that darken with the degree of threat. The
dividing lines are the quantiles of a logistic fit of published reports of
whale injury (with 0 meaning no injury or minor injury and 1 meaning severe or
fatal injury) to the base-10 logarithm of compressive stress.

Much can be learned by adjusting the sliders and examining the plotted output.
As an exercise, try setting to a particular ship mass of interest, and then to
slide the ship speed to higher and lower values, whilst monitoring the \"threat\"
panel. Next, try altering the blubber and sublayer thicknesses, e.g. using
smaller values to represent a strike at the whale mandible.  Having built some
intuition with these experiments, move on to altering the properties of the
ship, exploring the effect of changing ship speed, ship mass, and the width and
height of the impact zone."

# app for simulating the effect of a vessel strike on a whale
ui <- shiny::fluidPage(
    shiny::tags$style(shiny::HTML("body {font-family: 'Arial'; font-size: 12px; margin-left:1ex}")),
    shiny::fluidRow(shiny::radioButtons("instructions", "Instructions", choices = c("Hide", "Show"), selected = "Show", inline = TRUE)),
    shiny::conditionalPanel(
        condition = "input.instructions=='Show'",
        shiny::fluidRow(shiny::HTML(help))
    ),
    shiny::fluidRow(
        shiny::column(
            2,
            shiny::sliderInput("tmax", h6("Max time [s]"),
                ticks = FALSE,
                min = 0.1, max = 5, value = 1, step = 0.05
            ),
            shiny::sliderInput("ms", h6(tags$i("Ship mass [tonne]")),
                ticks = FALSE,
                min = 10, max = 500, value = 45, step = 1
            ),
            shiny::sliderInput("vs", h6(tags$i("Ship speed [knot]")),
                ticks = FALSE,
                min = 1, max = 30, value = 10, step = 0.1
            )
        ),
        shiny::column(
            2,
            shiny::sliderInput("Ly", h6("Impact width [m]"),
                ticks = FALSE,
                min = 0.1, max = 2, value = 1.15, step = 0.05
            ),
            shiny::sliderInput("Lz", h6("Impact height [m]"),
                ticks = FALSE,
                min = 0.1, max = 2, value = 1.15, step = 0.05
            )
        ),
        shiny::column(
            2,
            shiny::sliderInput("lw", h6("Whale length [m]"),
                ticks = FALSE,
                min = 5, max = 20, value = 13.7, step = 0.1
            ),
            shiny::sliderInput("l1", h6("Skin thickness [cm]"),
                ticks = FALSE,
                min = 1, max = 3, value = 2.5, step = 0.1
            ),
            shiny::sliderInput("l2", h6(tags$i("Blubber thickness [cm]")),
                ticks = FALSE,
                min = 5, max = 40, value = 16, step = 1
            )
        ),
        shiny::column(
            2,
            shiny::selectInput("species", "",
                choices = c(
                    "N. Atl. Right", "Blue", "Bryde", "Fin", "Gray", "Humpback", "Minke",
                    "Pac. Right", "Sei", "Sperm"
                ),
                selected = "N. Atl. Right"
            ),
            shiny::sliderInput("l3", h6(tags$i("Sublayer thickness [cm]")),
                ticks = FALSE,
                min = 5, max = 200, value = 112, step = 1
            ),
            shiny::sliderInput("l4", h6("Bone thickness [cm]"),
                ticks = FALSE,
                min = 5, max = 30, value = 10, step = 1
            )
        ),
        shiny::column(
            2,
            shiny::fileInput("loadFile", "Configuration", multiple = FALSE, accept = c("text/csv", ".csv")),
            shiny::actionButton("saveFile", "Save"),
            shiny::hr(),
            shiny::actionButton("quit", "Quit")
        ),
        shiny::column(
            2,
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
        )
    ),
    shiny::fluidRow(shiny::plotOutput("plot"))
)


#' Server for app, with standard arguments.
#' @param input A list created by the shiny server, with entries for slider settings, etc.
#' @param output A list of output entries, for plotting, etc.
#' @param session A list used for various purposes.
#' @importFrom utils write.csv
#' @importFrom shiny observeEvent reactiveValuesToList renderPlot showNotification updateSliderInput
server <- function(input, output, session) {
    observeEvent(input$quit, {
        shiny::stopApp()
    })

    observeEvent(input$saveFile, {
        file <- "boat_whale.csv"
        # Remove the load and save file items from the list of input items to save
        home <- normalizePath("~")
        fullfile <- paste0(home, .Platform$file.sep, file)
        config <- reactiveValuesToList(input)
        configNames <- names(config)
        w <- which("loadFile" == configNames | "saveFile" == configNames | "plot_panels" == configNames)
        config <- config[-w]
        # Convert ship speed from to m/s, from knots in the GUI
        config$vs <- whalestrike::knot2mps(config$vs)
        # Convert ship mass to kg, from tonne in the GUI
        config$ms <- 1e3 * config$ms
        # save in alphabetical order
        o <- order(names(config))
        write.csv(config[o], row.names = FALSE, file = fullfile)
        showNotification(paste0('Saved configuration to "', fullfile, '"'), type = "message")
    })

    observeEvent(input$loadFile, {
        config <- read.csv(input$loadFile$datapath)
        # Convert ship speed from m/s in the file, to knots in the GUI
        config$vs <- mps2knot(config$vs)
        # Convert ship mass from kg in file, to tonne in the GUI
        config$ms <- 1e-3 * config$ms
        # Insert individual thickness entries (one slider each)
        config$l1 <- config$l[1]
        config$l2 <- config$l[2]
        config$l3 <- config$l[3]
        config$l4 <- config$l[4]
        # FIXME: l1, l2 etc
        for (s in c("tmax", "ms", "lw", "vs", "Ly", "Lz", "l1", "l2", "l3", "l4")) {
            updateSliderInput(session, s, value = config[[s]])
        }
    })

    output$plot <- renderPlot(
        {
            # message("species: ", input$species)
            # aDefault <- whalestrike::parameters()$a
            # cat(file=stderr(), "input$a23=", input$a23, "\n")
            mw <- if (input$species == "N. Atl. Right") {
                whaleMassFromLength(input$lw, species = "N. Atl. Right Whale", model = "fortune2012")
            } else if (input$species == "Blue") {
                whaleMassFromLength(input$lw, species = "Blue Whale", model = "lockyer1976")
            } else if (input$species == "Bryde") {
                whaleMassFromLength(input$lw, species = "Bryde Whale", model = "lockyer1976")
            } else if (input$species == "Fin") {
                whaleMassFromLength(input$lw, species = "Fin Whale", model = "lockyer1976")
            } else if (input$species == "Gray") {
                whaleMassFromLength(input$lw, species = "Gray Whale", model = "lockyer1976")
            } else if (input$species == "Humpback") {
                whaleMassFromLength(input$lw, species = "Humpback Whale", model = "lockyer1976")
            } else if (input$species == "Minke") {
                whaleMassFromLength(input$lw, species = "Minke Whale", model = "lockyer1976")
            } else if (input$species == "Pac. Right") {
                whaleMassFromLength(input$lw, species = "Pac. Right Whale", model = "lockyer1976")
            } else if (input$species == "Sei") {
                whaleMassFromLength(input$lw, species = "Sei Whale", model = "lockyer1976")
            } else if (input$species == "Sperm") {
                whaleMassFromLength(input$lw, species = "Sperm Whale", model = "lockyer1976")
            } else {
                stop("programming error: cannot compute mass from length, for species '", input$species)
            }
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
            state <- list(xs = -(1 + parms$lsum), vs = whalestrike::knot2mps(input$vs), xw = 0, vw = 0)
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
        pointsize = 14
    ) # , height=500)
}

#' GUI app for interactive whale-strike simulations
#'
#' The `app()` function starts a GUI application that makes it easy to
#' run simple simulations and see the results in graphical form. Sliders
#' and buttons permit a fair degree of customization.  The application
#' has some build-in documentation, which supplements what can be found
#' in the \sQuote{Details} section of the present documentation.
#'
#' When `app()` is run, a window will appear within a few moments. At the top of that
#' is a textual introduction to the system, with a button to hide that information. Below is
#' a user-interaction area, with buttons and sliders that control the simulation
#' and the plotted output. Below that is a plotting area, the contents of which
#' depend on the configuration of the simulation as well as the user's selection
#' of items to display.
#'
#' The default setup, which is shown before the user alters any of the sliders,
#' etc., is a simulation of a small fishing boat, of mass 45 tonnes, moving at
#' speed 10 knots towards a whale of length 13.7m. (The whale length is used to
#' compute its mass, using a formula that is described by the output of typing
#' `help("whaleMassFromLength","whalestrike")` in an R console).
#'
#' Sliders are provided for setting certain key properties of the ship and the
#' whale, with italic labels for those properties that are deemed most likely to
#' be adjusted during simulations.  The details of these and the other parameters
#' are revealed by typing `help("parameters","whalestrike")` and
#' `help("strike","whalestrike")` in an R console.
#'
#' To the right of the sliders is a column of checkboxes that control the plotted
#' output. At startup, three of these boxes are ticked, yielding a display with
#' three panels showing the time history of the simulation
#' (`help("plot.strike","whalestrike")` provides details of the plots):
#' * The left-hand plot panel shows whale and boat location, the former with an
#'   indication of the interfaces between skin, blubber, sublayer, and bone.
#' * The middle panel shows the same information as the left one, but with a
#'   whale-centred coordinate system, and with labels for the components. This
#' makes it easier to see the degree to which the layers are compressed during the
#' impact.
#' * The right panel is an indication of the estimated threat to the four layers
#'   of the whale, with curves that are filled with grey for time intervals when
#' the impact stress (force/area) is less than the strength of the material in the
#' layer, and black for times when that threshold is exceeded.
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
#' mass. It does not affect the thicknesses of blubber or sublayer, the values of
#' which have been set up to represent a midsection strike on a North Atlantic
#' Right Whale.
#'
#' @param mode character value specifying the style to use.  Only
#' the value \code{"simple"} is permitted at present. This yields a
#' 3-panel plot, constructed by \code{\link{plot.strike}},
#' called with default arguments.
#'
#' @param options list containing options that are provided
#' to \code{\link[shiny]{shinyApp}}, which creates the GUI app.
#'
#' @export
#'
#' @importFrom shiny shinyApp
#'
#' @family interactive apps
#'
#' @author Dan Kelley
app <- function(mode = "simple", options = list(height = 500)) # NOTE: height has no effect
{
    if (mode == "simple") {
        shinyApp(ui = ui, server = server, options = options)
    } else {
        stop("unknown mode; only \"simple\" is permitted")
    }
}

shinyApp(ui, server)
