
# Setup things ----

library(shiny)

source("appFunctions.R")
cat("Reading in raw data... ")
raw_data_frame <- readRDS("raw_data_table")
cat("Done\n")
cat("Reading in MSMS data... ")
raw_msms_data <- readRDS("raw_msms_table")
cat("Done\n")
cat("Creating TIC... ")
tic <- raw_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))
cat("Done\n")
cat("Reading in metadata... ")
falkor_metadata <- read.csv("falkor_metadata.csv")
cat("Done\n")

# Functions ----

get_EIC <- function(raw_data_frame, mass, ppm=5, 
                    mdframe=falkor_metadata){
  raw_data_frame %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
}

get_Spectrum <- function(raw_data_frame, scan, ret_window=1){
  raw_data_frame %>% 
    filter(rt>scan-ret_window/2&rt<scan+ret_window/2) %>% 
    mutate(mz=round(mz*1000)/1000) %>% # Round to group nearby scans
    group_by(mz) %>% 
    summarize(TIS=sum(int))
}

plotGivenEIC <- function(eic, plotby="depth", plottic=TRUE, tic=NULL){
  if(plottic){
    tic$int <- (tic$int/max(tic$int))*max(eic$int)
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green"), 
                                  unique(eic[[plotby]]))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", 
                line=list(color="black"),
                name="TIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("red", "blue", "green"), unique(eic[[plotby]]))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  }
}

plotGivenSpectrum <- function(spectrum){
  plot_ly(source = "TIS") %>%
    add_trace(data=spectrum, x=~mz, y=~TIS, type = "bar", 
              marker=list(color="black"), hoverinfo="none") %>%
    add_trace(data=spectrum, x=~mz, y=~TIS, marker=list(color="black"),
              type = "scatter", mode="markers") %>%
    layout(xaxis = list(title = "m/z"),
           yaxis = list(title = "Intensity",
                        fixedrange = TRUE))
}


# UI ----
ui <- fluidPage(
    titlePanel("mavenMimic"),
    sidebarLayout(
        sidebarPanel(
            numericInput(inputId = "given_mz",
                        "m/z:",
                        value = 118.0868),
            numericInput(inputId = "given_ppm",
                         "+/- (ppm)",
                         value = 5),
            radioButtons(inputId = "treatment", 
                         label = "Color by which?", 
                         selected = "depth", 
                         choiceNames = c("Depth", "Spin direction", "Time"),
                         choiceValues = c("depth", "spindir", "time")),
            checkboxInput(inputId = "user_plottic",
                          label = "Plot a TIC on the graph?",
                          value = TRUE)
        ),

        mainPanel(
          plotlyOutput("chrom", height = "80%"),
          # verbatimTextOutput("debug"),
          # tableOutput("debug"),
          plotlyOutput("TIS", height = "80%")
        )
    )
)

# Server ----
server <- function(input, output, session) {
  current_mass <- reactiveVal(NULL)
  observe({
    current_mass(event_data(event = "plotly_click", source = "TIS")$x)
  })
  observeEvent(input$given_mz, {
    current_mass(input$given_mz)
  })
  observeEvent(current_mass(), {
    updateNumericInput(session, "given_mz", value = current_mass())
  })
  
  
  given_EIC <- reactive({
    get_EIC(raw_data_frame = raw_data_frame, 
            mass = current_mass(), 
            ppm = input$given_ppm, 
            mdframe = falkor_metadata)
  })
  
  given_Spectrum <- reactive({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    get_Spectrum(raw_data_frame = raw_data_frame, scan=EIC_data$x)
  })
  
  # output$debug <- renderTable({head(given_Spectrum())})
  
  output$chrom <- renderPlotly({
    plotGivenEIC(eic = given_EIC(), 
                 plotby = input$treatment,
                 plottic = input$user_plottic, 
                 tic = tic)
  })

  output$TIS <- renderPlotly({
    plotGivenSpectrum(given_Spectrum())
  })
}

# Run the application ----
shinyApp(ui = ui, server = server)
