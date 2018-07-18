#' User interface for app
ui <- fluidPage(tags$style(HTML("body {font-family: 'Arial'; font-size: 12px;}")),
                fluidRow(column(2,
                                radioButtons("gammaType", h6("Sublayer type"),
                                             c("Tissue"=1, "Bone"=2), inline=TRUE)),
                         column(2,
                                fileInput("loadFile", "Configuration",
                                          multiple=FALSE,
                                          accept=c("text/csv",
                                        #"text/comma-separated-values,text/plain",
                                                   ".csv"))),
                         column(2,
                                actionButton("saveFile", "Save"))),
                ##actionButton("saveFile", "Save", icon("save")))),
                fluidRow(column(2,
                                sliderInput("tmax",  h6("log10(interval [s])"), tick=FALSE,
                                            min=-1,  max=1, value=log10(1), step=0.01),
                                sliderInput("ms",  h6("Ship mass [tonne]"), tick=FALSE,
                                            min=10, max=40,  value=20, step=0.5),
                                sliderInput("vs", h6("Ship speed [knot]"), tick=FALSE,
                                            min=1,  max=20,  value=10, step=0.1)),
                         column(2,
                                sliderInput("Ly",  h6("Impact width [m]"), tick=FALSE,
                                            min=0.1, max=2,  value=0.5, step=0.05),
                                sliderInput("Lz",  h6("Impact height [m]"), tick=FALSE,
                                            min=0.1, max=2,  value=1.5, step=0.05),
                                sliderInput("lw",  h6("Whale length [m]"), tick=FALSE,
                                            min=5,  max=15, value=13.7, step=0.1)),
                         column(2,
                                sliderInput("alpha", h6("Skin thickness [m]"), tick=FALSE,
                                            min=0.01, max=0.03, value=0.025, step=0.001),
                                sliderInput("Ealpha", h6("Skin modulus [MPa]"), tick=FALSE,
                                            min=10, max=30, value=17.80, step=0.5),
                                sliderInput("UTSalpha", h6("Skin strength [MPa]"), tick=FALSE,
                                            min=10, max=30, value=19.56, step=0.5)),
                         column(2,
                                sliderInput("theta", h6("Skin theta [deg]"), tick=FALSE,
                                            min=0, max=85, value=45)),
                         column(2,
                                sliderInput("beta", h6("Blubber thickness [m]"), tick=FALSE,
                                            min=0.05, max=.4, value=0.2, step=0.01),
                                sliderInput("Ebeta", h6("Blub. modulus [MPa]"), tick=FALSE,
                                            min=0.1, max=3, value=0.6, step=0.1),
                                sliderInput("UTSbeta", h6("Blub. strength [MPa]"), tick=FALSE,
                                            min=0.1, max=3, value=(0.8/1.2)*0.6), step=0.01),
                         conditionalPanel(condition="input.gammaType == 1",
                                          column(2,
                                                 sliderInput("gamma1", h6("Tissue thickness [m]"),
                                                             tick=FALSE,
                                                             min=0.2, max=1.5, value=0.5),
                                                 sliderInput("Egamma1", h6("Tissue modulus [MPa]"),
                                                             tick=FALSE,
                                                             min=0.2, max=0.8, value=0.425294), # Raymond (2007) Table 2.6 for value
                                                 sliderInput("UTSgamma1", h6("Tissue strength [MPa]"),
                                                             tick=FALSE,
                                                             min=0.2, max=1.5, value=(0.8/1.2)*0.4))),
                         conditionalPanel(condition="input.gammaType == 2",
                                          column(2,
                                                 sliderInput("gamma2", h6("Bone thickness [m]"),
                                                             tick=FALSE,
                                                             min=0.01, max=0.20, value=0.12),
                                                 sliderInput("Egamma2", h6("Bone modulus [MPa]"),
                                                             tick=FALSE,
                                                             min=557.9, max=1360.4, value=854.2),# Raymond (2007) Table 2.3
                                                 sliderInput("UTSgamma2", h6("Bone strength [MPa]"),
                                                             tick=FALSE,
                                                             min=12.9, max=51.7, value=22.9)))), # Raymond (2007) Table 2.3
                fluidRow(plotOutput("plot", height=400)))


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
                 w <- which("loadFile" == configNames | "saveFile" == configNames)
                 config <- config[-w]
                 ## Convert ship speed from to m/s, from knots in the GUI
                 config$vs <- 0.514444 * config$vs
                 ## Convert ship mass to kg, from tonne in the GUI
                 config$ms <- 1e3 * config$ms
                 ## Convert some things from Pa in the file, to MPa in the GUI
                 config$Ealpha <- 1e6 * config$Ealpha
                 config$Ebeta <- 1e6 * config$Ebeta
                 config$Egamma1 <- 1e6 * config$Egamma1
                 config$Egamma2 <- 1e6 * config$Egamma2
                 config$UTSalpha <- 1e6 * config$UTSalpha
                 config$UTSbeta <- 1e6 * config$UTSbeta
                 config$UTSgamma1 <- 1e6 * config$UTSgamma1
                 config$UTSgamma2 <- 1e6 * config$UTSgamma2
                 ## Add in un-suffixed items for use by whalestrike::parameters()
                 if (config$gammaType == 1) {
                     config$gamma <- config$gamma1 # convert from MPa to Pa
                     config$Egamma <- config$Egamma1 # convert from MPa to Pa
                     config$UTSgamma <- config$UTSgamma1 # convert from MPa to Pa
                 } else {
                     config$gamma <- config$gamma2 # convert from MPa to Pa
                     config$Egamma <- config$Egamma2 # convert from MPa to Pa
                     config$UTSgamma <- config$UTSgamma2 # convert from MPa to Pa
                 }
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
                 ## Convert some things from Pa in the file, to MPa in the GUI
                 config$Ealpha <- 1e-6 * config$Ealpha
                 config$Ebeta <- 1e-6 * config$Ebeta
                 config$Egamma1 <- 1e-6 * config$Egamma1
                 config$Egamma2 <- 1e-6 * config$Egamma2
                 config$UTSalpha <- 1e-6 * config$UTSalpha
                 config$UTSbeta <- 1e-6 * config$UTSbeta
                 config$UTSgamma1 <- 1e-6 * config$UTSgamma1
                 config$UTSgamma2 <- 1e-6 * config$UTSgamma2
                 for (s in c("tmax", "ms", "lw", "vs", "Ly", "Lz",
                             "alpha", "Ealpha", "UTSalpha",
                             "theta",
                             "beta", "Ebeta", "UTSbeta",
                             "gammaType",
                             "gamma1", "Egamma1", "UTSgamma1",
                             "gamma2", "Egamma2", "UTSgamma2"))
                     updateSliderInput(session, s, value=config[[s]])
                })
    output$plot <- renderPlot({
        lm <- 10
        parms <- parameters(ms=1000*input$ms, Ss=shipAreaFromMass(1000*input$ms),
                            Ly=input$Ly, Lz=input$Lz,
                            lw=input$lw,
                            mw=whaleMassFromLength(input$lw),
                            Sw=whaleAreaFromLength(input$lw, "wetted"),
                            alpha=input$alpha, # already in m
                            Ealpha=input$Ealpha*1e6, # MPa to Pa
                            UTSalpha=input$UTSalpha*1e6, # MPa to Pa
                            theta=input$theta, # angle retained in degree by whalestrike
                            beta=input$beta,   # already in m
                            Ebeta=input$Ebeta*1e6, # MPa to Pa
                            UTSbeta=input$UTSbeta*1e6, # MPa to Pa
                            gamma=if (input$gammaType == 1) input$gamma1 else input$gamma2,
                            Egamma=1e6*if (input$gammaType == 1) input$Egamma1 else input$Egamma2,
                            UTSgamma=1e6*if (input$gammaType == 1) input$UTSgamma1 else input$UTSgamma2)
        state <- c(xs=-(1 + parms$alpha + parms$beta + parms$gamma), vs=input$vs * 0.514444, xw=0, vw=0)
        t <- seq(0, 10^input$tmax, length.out=500)
        sol <- strike(t, state, parms)
        par(mfcol=c(1, 3), mar=c(3, 3, 1, 2), mgp=c(2, 0.7, 0), cex=1)
        plot(sol) # c("location", "section", "threat"))
    }, pointsize=12, height=350)
}

#' Run a GUI app for interactive simulations
#' @param mode Character string specifying the style to use.  Only
#' the value \code{"simple"} is permitted at present. This yields a
#' 3-panel plot, constructed by \code{\link{plot.strike}},
#' called with default arguments.
app <- function(mode="simple")
{
    if (mode == "simple")
        shinyApp(ui, server)
    else
        stop("unknown mode; only \"simple\" is permitted")
}
