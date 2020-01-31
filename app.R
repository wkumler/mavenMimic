
# Setup things ----

library(shiny)

source("appFunctions.R")

raw_data_frame <- readRDS("raw_data_frame")
tic <- raw_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))

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
           plotlyOutput("chrom", height = "90%"),
           #verbatimTextOutput("debug"),
           plotlyOutput("TIS", height = "90%")
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
