#' User interface for app
ui <- fluidPage(tags$style(HTML("body {font-family: 'Arial'; font-size: 12px; margin-left:1ex}")),
                fluidRow(radioButtons("instructions", "Instructions", choices=c("Hide", "Show"), selected="Show", inline=TRUE)),
                conditionalPanel(condition="input.instructions=='Show'", fluidRow(includeMarkdown("app_help.md"))),
                fluidRow(column(2,
                                sliderInput("tmax",  h6("Max time [s]"), ticks=FALSE,
                                            min=0.1,  max=5, value=1, step=0.05),
                                sliderInput("ms",  HTML("<font color=\"FF0000\">Ship mass [tonne]</font>"), ticks=FALSE,
                                            min=10, max=500,  value=45, step=0.5),
                                sliderInput("vs", HTML("<font color=\"FF0000\">Ship speed [knot]</font>"), ticks=FALSE,
                                            min=0.5,  max=30,  value=10, step=0.5)),
                         column(2,
                                sliderInput("Ly",  h6("Impact width [m]"), ticks=FALSE,
                                            min=0.1, max=2,  value=1.15, step=0.05),
                                sliderInput("Lz",  h6("Impact height [m]"), ticks=FALSE,
                                            min=0.1, max=2,  value=1.15, step=0.05)),
                                ##selectInput("species", "Whale species:",
                                ##            choices=c("N. Atl. Right Whale",
                                ##                      "NOTHING ELSE CODED YET"))),
                         column(2,
                                sliderInput("lw",  h6("Right whale length [m]"), ticks=FALSE,
                                            min=5,  max=20, value=13.7, step=0.1),
                                ##sliderInput("theta", h6("Skin theta [deg]"), ticks=FALSE,
                                ##            min=30, max=70, value=55, step=1),
                                sliderInput("l1", h6("Skin thickness [m]"), ticks=FALSE,
                                            min=0.01, max=0.03, value=0.025, step=0.001),
                                sliderInput("l2", h6("Blubber thickness [m]"), ticks=FALSE,
                                            min=0.05, max=.4, value=0.16, step=0.01)),
                         column(2,
                                ## default: a[2]=a[3]=1.58e5 Pa
                                sliderInput("a23", HTML("<font color=\"FF0000\">Blubber/Sublayer compressibility [MPa]</font>"), ticks=FALSE,
                                            min=0.100, max=0.200, value=0.158, step=0.01),
                                sliderInput("l3", HTML("<font color=\"FF0000\">Sublayer thickness[m]</font>"), ticks=FALSE,
                                            min=0.05, max=2, value=1.12, step=0.01),
                                sliderInput("l4", h6("Bone thickness [m]"), ticks=FALSE,
                                            min=0.05, max=.3, value=0.1, step=0.01)),
                         column(2,
                                fileInput("loadFile", "Configuration", multiple=FALSE, accept=c("text/csv", ".csv")),
                                actionButton("saveFile", "Save")),
                         column(2,
                                checkboxGroupInput("plot_panels", "",
                                                   choices=c("location", "section", "threat", "whale acceleration",
                                                             "blubber thickness", "sublayer thickness",
                                                             "whale water force", "reactive forces", "skin stress",
                                                             "compression stress", "values"),
                                                   selected=c("location", "section", "threat")))),
                fluidRow(plotOutput("plot")))# , height='200px'))) # NOTE: height has no effect


#' Server for app, with standard arguments.
#' @param input A list created by the shiny server, with entries for slider settings, etc.
#' @param output A list of output entries, for plotting, etc.
#' @param session A list used for various purposes.
server <- function(input, output, session)
{
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
                 config$vs <- knot2mps(config$vs)
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
                 config$vs <- (1/knot2mps(1)) * config$vs
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
        aDefault <- parameters()$a
        cat(file=stderr(), "input$a23=", input$a23, "\n")
        parms <- parameters(ms=1000*input$ms, Ss=shipAreaFromMass(1000*input$ms),
                            Ly=input$Ly, Lz=input$Lz,
                            lw=input$lw,
                            mw=whaleMassFromLength(input$lw, species="N. Atl. Right Whale"),
                            Sw=whaleAreaFromLength(input$lw, species="N. Atl. Right Whale", "wetted"),
                            l=c(input$l1, input$l2, input$l3, input$l4),
                            a=c(aDefault[1], 1e6*input$a23, 1e6*input$a23, aDefault[4]))
        state <- list(xs=-(1 + parms$lsum), vs=knot2mps(input$vs), xw=0, vw=0)
        t <- seq(0, input$tmax, length.out=2000)
        sol <- strike(t, state, parms)
        if (sol$refinedGrid)
            showNotification("Auto-refined grid to capture acceleration peak")
        npanels <- length(input$plot_panels)
        nrows <- floor(sqrt(npanels))
        ncols <- ceiling(npanels / nrows)
        par(mfrow=c(nrows,ncols), mar=c(3,3,1,2), mgp=c(2,0.7,0), cex=1)
        for (which in input$plot_panels)
            plot(sol, which=which)
    }, pointsize=12)#, height=500)
}

#' Run a GUI app for interactive simulations
#' @param mode Character string specifying the style to use.  Only
#' the value \code{"simple"} is permitted at present. This yields a
#' 3-panel plot, constructed by \code{\link{plot.strike}},
#' called with default arguments.
#' @param options List containing options that are provided
#' to \code{\link[shiny]{shinyApp}}, which creates the GUI app.
app <- function(mode="simple", options=list(height=500)) # NOTE: height has no effect
{
    if (mode == "simple")
        shinyApp(ui=ui, server=server, options=options)
    else
        stop("unknown mode; only \"simple\" is permitted")
}
