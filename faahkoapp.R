
# Setup things ----

library(shiny)
library(dplyr)
library(plotly)

cat("Reading in raw data... ")
alldata <- readRDS("Data/faahko_MS1_data_frame")
cat("Done\n")
cat("Creating TIC... ")
tic <- alldata %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))
cat("Done\n")
cat("Reading in metadata... ")
metadata <- read.csv("Data/faahko_metadata.csv")
cat("Done\n")

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

get_EIC <- function(alldata, mass, ppm=100, 
                    mdframe=metadata){
  alldata %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
}

get_Spectrum <- function(alldata, scan, ret_window=3){
  alldata %>% 
    filter(rt>scan-ret_window/2&rt<scan+ret_window/2) %>% 
    mutate(mz=round(mz*1000)/1000) %>% # Round to group nearby scans
    group_by(mz) %>% 
    summarize(TIS=sum(int))
}

plotGivenEIC <- function(eic, plotby="treatment", plottic=TRUE, tic=NULL,
                         current_mass=335, ppm=100){
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
             yaxis = list(title = "Intensity"),
             title = paste(round(pmppm(current_mass, ppm = ppm), digits = 4), 
                           collapse = " - "))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("red", "blue", "green"), unique(eic[[plotby]]))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = paste(round(pmppm(current_mass, ppm = ppm), digits = 4), 
                           collapse = " - "))
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
                   value = 335),
      numericInput(inputId = "given_ppm",
                   "+/- (ppm)",
                   value = 100),
      checkboxInput(inputId = "user_plottic",
                    label = "Plot a TIC on the graph?",
                    value = FALSE)
    ),
    
    mainPanel(
      h2("Extracted ion chromatogram: "),
      plotlyOutput(outputId = "chrom", height = "300px"),
      h2("Ions identified in chosen scan: "),
      plotlyOutput(outputId = "TIS", height = "300px")
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
  # observeEvent(current_mass(), {
  #   updateNumericInput(session, "given_mz", value = current_mass())
  # })
  
  
  given_EIC <- reactive({
    get_EIC(alldata = alldata, 
            mass = current_mass(), 
            ppm = input$given_ppm, 
            mdframe = metadata)
  })
  
  given_Spectrum <- reactive({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    get_Spectrum(alldata = alldata, scan=EIC_data$x)
  })
  
  output$chrom <- renderPlotly({
    plotGivenEIC(eic = given_EIC(), 
                 plotby = "treatment",
                 plottic = input$user_plottic, 
                 tic = tic,
                 current_mass = current_mass(),
                 ppm = input$given_ppm)
  })
  
  output$TIS <- renderPlotly({
    plotGivenSpectrum(given_Spectrum())
  })
}


# Run the application ----
shinyApp(ui = ui, server = server)
