# Setup things ----

library(shiny)
library(plotly)
library(dplyr)
library(RSQLite)

#load("raw_data_frame")
tic <- raw_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

plotGivenEIC <- function(mass, ppm=5, df=raw_data_frame, plottic=TRUE){
  eic <- df %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(file, rt) %>% 
    summarize(int=sum(int)) %>%
    mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
    mutate(spindir=c("Cyclone", "Anticyclone")[(1-ceiling(file/12)%%2)+1])
  if(plottic){
    tic$int <- (tic$int/max(tic$int))*max(eic$int)
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~spindir, opacity = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone"))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity",
                          fixedrange = TRUE))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~spindir, alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone"))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity",
                          fixedrange = TRUE))
  }
}

plotGivenScan <- function(ret, window=1, df=raw_data_frame){
  scandata <- df %>% 
    filter(rt>ret-window/2&rt<ret+window/2) %>% 
    mutate(mz=round(mz*100)/100) %>%
    group_by(mz) %>% 
    summarize(TIS=sum(int))
  plot_ly(source = "TIS") %>%
    add_trace(data=scandata, x=~mz, y=~TIS, type = "bar", 
              marker=list(color="black"), hoverinfo="none") %>%
    add_trace(data=scandata, x=~mz, y=~TIS, marker=list(color="black"),
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
                         value = 50)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("chrom"),
           #verbatimTextOutput("debug"),
           plotlyOutput("TIS")
        )
    )
)

# Server ----
server <- function(input, output) {

    output$chrom <- renderPlotly({
        suppressWarnings(TIS_data <- event_data(event = "plotly_click", source="TIS"))
        if(is.null(TIS_data)){
            plotGivenEIC(input$given_mz, ppm=input$given_ppm)
        } else {
            plotGivenEIC(TIS_data$x, ppm=input$given_ppm)
        }
    })
    
    output$debug <- renderPrint({
      event_data(event = "plotly_click", source = "EIC")
    })
    
    output$TIS <- renderPlotly({
        EIC_data <- event_data(event = "plotly_click", source = "EIC")
        if (is.null(EIC_data)) {
            plotly_empty(type="scatter", mode="markers")
        } else {
            plotGivenScan(ret = EIC_data$x, window = 1)
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
