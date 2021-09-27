## vim:textwidth=128:expandtab:shiftwidth=4:softtabstop=4

## app for simulating the effect of a vessel strike on a whale
ui <- fluidPage(tags$style(HTML("body {font-family: 'Arial'; font-size: 12px; margin-left:1ex}")),
  fluidRow(radioButtons("instructions", "Instructions", choices=c("Hide", "Show"), selected="Show", inline=TRUE)),
  conditionalPanel(condition="input.instructions=='Show'",
    fluidRow(includeMarkdown(system.file("extdata", "app_help.md", package="whalestrike")))),
  fluidRow(column(2,
      sliderInput("tmax",  h6("Max time [s]"), ticks=FALSE,
        min=0.1,  max=5, value=1, step=0.05),
      sliderInput("ms",  h6(tags$i("Ship mass [tonne]")), ticks=FALSE,
        min=10, max=500,  value=45, step=1),
      sliderInput("vs", h6(tags$i("Ship speed [knot]")), ticks=FALSE,
        min=1,  max=30,  value=10, step=0.1)),
    column(2,
      sliderInput("Ly",  h6("Impact width [m]"), ticks=FALSE,
        min=0.1, max=2,  value=1.15, step=0.05),
      sliderInput("Lz",  h6("Impact height [m]"), ticks=FALSE,
        min=0.1, max=2,  value=1.15, step=0.05)),
    column(2,
      sliderInput("lw",  h6("Whale length [m]"), ticks=FALSE,
        min=5,  max=20, value=13.7, step=0.1),
      sliderInput("l1", h6("Skin thickness [cm]"), ticks=FALSE,
        min=1, max=3, value=2.5, step=0.1),
      sliderInput("l2", h6(tags$i("Blubber thickness [cm]")), ticks=FALSE,
        min=5, max=40, value=16, step=1)),
    column(2,
      selectInput("species", "",
        choices= c("N. Atl. Right", "Blue", "Bryde", "Fin", "Gray", "Humpback", "Minke",
          "Pac. Right", "Sei", "Sperm"),
        selected="N. Atl. Right"),
      sliderInput("l3", h6(tags$i("Sublayer thickness [cm]")), ticks=FALSE,
        min=5, max=200, value=112, step=1),
      sliderInput("l4", h6("Bone thickness [cm]"), ticks=FALSE,
        min=5, max=30, value=10, step=1)),
    column(2,
      fileInput("loadFile", "Configuration", multiple=FALSE, accept=c("text/csv", ".csv")),
      actionButton("saveFile", "Save"),
      hr(),
      actionButton("quit", "Quit")),
    column(2,
      checkboxGroupInput("plot_panels", "",
        choices=c("location",
          "section",
          "lethality index",
          "threat",
          "whale acceleration",
          "whale water force",
          "reactive forces",
          "skin stress",
          "compression stress",
          "values"),
        selected=c("location", "section", "lethality index")))),
  fluidRow(plotOutput("plot")))


#' Server for app, with standard arguments.
#' @param input A list created by the shiny server, with entries for slider settings, etc.
#' @param output A list of output entries, for plotting, etc.
#' @param session A list used for various purposes.
#' @importFrom utils write.csv
#' @importFrom shiny observeEvent reactiveValuesToList renderPlot showNotification updateSliderInput
#' @md
server <- function(input, output, session)
{
    observeEvent(input$quit, {
        shiny::stopApp()
  })

    observeEvent(input$saveFile, {
        file <- "boat_whale.csv"
        ## Remove the load and save file items from the list of input items to save
        home <- normalizePath("~")
        fullfile <- paste0(home, .Platform$file.sep, file)
        config <- reactiveValuesToList(input)
        configNames <- names(config)
        w <- which("loadFile" == configNames | "saveFile" == configNames | "plot_panels" == configNames)
        config <- config[-w]
        ## Convert ship speed from to m/s, from knots in the GUI
        config$vs <- whalestrike::knot2mps(config$vs)
        ## Convert ship mass to kg, from tonne in the GUI
        config$ms <- 1e3 * config$ms
        ## save in alphabetical order
        o <- order(names(config))
        write.csv(config[o], row.names=FALSE, file=fullfile)
        showNotification(paste0('Saved configuration to "', fullfile, '"'), type="message")
  })

    observeEvent(input$loadFile, {
        config <- read.csv(input$loadFile$datapath)
        ## Convert ship speed from m/s in the file, to knots in the GUI
        config$vs <- mps2knot(config$vs)
        ## Convert ship mass from kg in file, to tonne in the GUI
        config$ms <- 1e-3 * config$ms
        ## Insert individual thickness entries (one slider each)
        config$l1 <- config$l[1]
        config$l2 <- config$l[2]
        config$l3 <- config$l[3]
        config$l4 <- config$l[4]
        ## FIXME: l1, l2 etc
        for (s in c("tmax", "ms", "lw", "vs", "Ly", "Lz",
                "l1", "l2", "l3", "l4"))
            updateSliderInput(session, s, value=config[[s]])
  })

    output$plot <- renderPlot({
        ##message("species: ", input$species)
        ## aDefault <- whalestrike::parameters()$a
        ##cat(file=stderr(), "input$a23=", input$a23, "\n")
        mw <- if (input$species == "N. Atl. Right") {
            whaleMassFromLength(input$lw, species="N. Atl. Right Whale", model="fortune2012")
        } else if (input$species == "Blue") {
            whaleMassFromLength(input$lw, species="Blue Whale", model="lockyer1976")
        } else if (input$species == "Bryde") {
            whaleMassFromLength(input$lw, species="Bryde Whale", model="lockyer1976")
        } else if (input$species == "Fin") {
            whaleMassFromLength(input$lw, species="Fin Whale", model="lockyer1976")
        } else if (input$species == "Gray") {
            whaleMassFromLength(input$lw, species="Gray Whale", model="lockyer1976")
        } else if (input$species == "Humpback") {
            whaleMassFromLength(input$lw, species="Humpback Whale", model="lockyer1976")
        } else if (input$species == "Minke") {
            whaleMassFromLength(input$lw, species="Minke Whale", model="lockyer1976")
        } else if (input$species == "Pac. Right") {
            whaleMassFromLength(input$lw, species="Pac. Right Whale", model="lockyer1976")
        } else if (input$species == "Sei") {
            whaleMassFromLength(input$lw, species="Sei Whale", model="lockyer1976")
        } else if (input$species == "Sperm") {
            whaleMassFromLength(input$lw, species="Sperm Whale", model="lockyer1976")
        } else {
            stop("programming error: cannot compute mass from length, for species '", input$species)
        }
        parms <- whalestrike::parameters(ms=1000*input$ms,
            Ss=shipAreaFromMass(1000*input$ms),
            Ly=input$Ly,
            Lz=input$Lz,
            lw=input$lw,
            mw=mw,
            Sw=whalestrike::whaleAreaFromLength(input$lw, species="N. Atl. Right Whale", "wetted"),
            l=c(input$l1/100, input$l2/100, input$l3/100, input$l4/100))
        state <- list(xs=-(1 + parms$lsum), vs=whalestrike::knot2mps(input$vs), xw=0, vw=0)
        t <- seq(0, input$tmax, length.out=2000)
        sol <- strike(t, state, parms)
        if (sol$refinedGrid)
            showNotification("Auto-refined grid to capture acceleration peak")
        npanels <- length(input$plot_panels)
        nrows <- floor(sqrt(npanels))
        ncols <- ceiling(npanels / nrows)
        par(mfrow=c(nrows,ncols), mar=c(3,3,2.5,2), mgp=c(2,0.7,0), cex=1)
        for (which in input$plot_panels)
            plot(sol, which=which)
    }, pointsize=12)#, height=500)
}

#' Run a GUI app for interactive simulations
#'
#' @param mode Character string specifying the style to use.  Only
#' the value \code{"simple"} is permitted at present. This yields a
#' 3-panel plot, constructed by \code{\link{plot.strike}},
#' called with default arguments.
#'
#' @param options List containing options that are provided
#' to \code{\link[shiny]{shinyApp}}, which creates the GUI app.
#'
#' @author Dan Kelley
#'
#' @export
#'
#' @importFrom shiny shinyApp
#' @md
app <- function(mode="simple", options=list(height=500)) # NOTE: height has no effect
{
    if (mode == "simple")
        shinyApp(ui=ui, server=server, options=options)
    else
        stop("unknown mode; only \"simple\" is permitted")
}

shinyApp(ui, server)
