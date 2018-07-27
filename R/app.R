#' User interface for app
ui <- fluidPage(tags$style(HTML("body {font-family: 'Arial'; font-size: 12px;}")),
                fluidRow(column(2,
                                sliderInput("ltmax",  h6("log10(interval [s])"), tick=FALSE,
                                            min=-1,  max=1, value=log10(1), step=0.01),
                                sliderInput("ms",  h6("Ship mass [tonne]"), tick=FALSE,
                                            min=10, max=311,  value=20, step=0.5),
                                sliderInput("vs", h6("Ship speed [knot]"), tick=FALSE,
                                            min=1,  max=20,  value=10, step=0.1)),
                         column(2,
                                sliderInput("Ly",  h6("Impact width [m]"), tick=FALSE,
                                            min=0.1, max=2,  value=0.5, step=0.05),
                                sliderInput("Lz",  h6("Impact height [m]"), tick=FALSE,
                                            min=0.1, max=2,  value=1.0, step=0.05),
                                sliderInput("lw",  h6("Whale length [m]"), tick=FALSE,
                                            min=5,  max=15, value=13.7, step=0.1)),
                         column(2,
                                sliderInput("theta", h6("Skin theta [deg]"), tick=FALSE,
                                            min=40, max=70, value=55, step=1),
                                sliderInput("l1", h6("Skin thickness [m]"), tick=FALSE,
                                            min=0.01, max=0.03, value=0.025, step=0.001),
                                sliderInput("l2", h6("Blubber thickness [m]"), tick=FALSE,
                                            min=0.05, max=.4, value=0.16, step=0.01)),
                         column(2,
                                sliderInput("l3", h6("Sub-layer thickness [m]"), tick=FALSE,
                                            min=0.05, max=2, value=1.12, step=0.01),
                                sliderInput("l4", h6("Bone thickness [m]"), tick=FALSE,
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
                fluidRow(plotOutput("plot", height=500)))


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
                 config$vs <- 0.514444 * config$vs
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
                 config$vs <- (1/0.514444) * config$vs
                 ## Convert ship mass from kg in file, to tonne in the GUI
                 config$ms <- 1e-3 * config$ms
                 ## Insert individual thickness entries (one slider each)
                 config$l1 <- config$l[1]
                 config$l2 <- config$l[2]
                 config$l3 <- config$l[3]
                 config$l4 <- config$l[4]
                 ## FIXME: l1, l2 etc
                 for (s in c("ltmax", "ms", "lw", "vs", "Ly", "Lz", "theta",
                             "l1", "l2", "l3", "l4"))
                     updateSliderInput(session, s, value=config[[s]])
                })
    output$plot <- renderPlot({
        lm <- 10
        parms <- parameters(ms=1000*input$ms, Ss=shipAreaFromMass(1000*input$ms),
                            Ly=input$Ly, Lz=input$Lz,
                            lw=input$lw,
                            mw=whaleMassFromLength(input$lw),
                            Sw=whaleAreaFromLength(input$lw, "wetted"),
                            l=c(input$l1, input$l2, input$l3, input$l4),
                            theta=input$theta) # in degrees; 0 means no bevel
        state <- c(xs=-(1 + parms$lsum), vs=input$vs * 0.514444, xw=0, vw=0)
        t <- seq(0, 10^input$ltmax, length.out=500)
        sol <- strike(t, state, parms)
        npanels <- length(input$plot_panels)
        nrows <- floor(sqrt(npanels))
        ncols <- ceiling(npanels / nrows)
        par(mfrow=c(nrows,ncols), mar=c(3,3,1,2), mgp=c(2,0.7,0), cex=1)
        for (which in input$plot_panels)
            plot(sol, which=which)
    }, pointsize=12, height=500)
}

#' Run a GUI app for interactive simulations
#' @param mode Character string specifying the style to use.  Only
#' the value \code{"simple"} is permitted at present. This yields a
#' 3-panel plot, constructed by \code{\link{plot.strike}},
#' called with default arguments.
#' @param options List containing options that are provided
#' to \code{\link[shiny]{shinyApp}}, which creates the GUI app.
app <- function(mode="simple", options=list(height=800))
{
    if (mode == "simple")
        shinyApp(ui=ui, server=server, options=options)
    else
        stop("unknown mode; only \"simple\" is permitted")
}
