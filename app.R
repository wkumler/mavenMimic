
# Setup things ----

library(shiny)

source("appFunctions.R")

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
                         value = 50),
            radioButtons(inputId = "treatment", 
                         label = "Color by which?", 
                         selected = "depth", 
                         choiceNames = c("Depth", "Spin direction"),
                         choiceValues = c("depth", "spindir")),
            checkboxInput(inputId = "user_plottic",
                          label = "Plot a TIC on the graph?",
                          value = TRUE)
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
  observeEvent(input$given_mz, {TIS_data <- NULL})

  output$chrom <- renderPlotly({
      suppressWarnings(TIS_data <- event_data(event = "plotly_click", source="TIS"))
      if(is.null(TIS_data)){
          plotGivenEIC(input$given_mz, ppm=input$given_ppm, 
                       plotby = input$treatment, plottic = input$user_plottic)
      } else {
          plotGivenEIC(TIS_data$x, ppm=input$given_ppm, 
                       plotby = input$treatment, plottic = input$user_plottic)
      }
  })
  
  output$debug <- renderPrint({
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
