# Setup things ----

library(shiny)
library(plotly)
library(dplyr)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

plotGivenEIC <- function(mass, ppm=5, df=raw_data_frame, plotTIC=FALSE){
    eic <- df %>% 
        filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
        group_by(file, rt) %>% 
        summarize(EIC_int=sum(int)) %>%
        mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
        mutate(spindir=c("Cyclone", "Anticyclone")[(1-ceiling(file/12)%%2)+1])
    if(plotTIC){
        TIC_df <- df %>% 
            mutate(rt=round(rt)) %>% 
            group_by(rt) %>% 
            summarize(TIC=sum(int)) %>%
            mutate(TIC=(TIC/max(TIC))*max(eic$EIC_int))
        plot_ly(source = "EIC") %>%
            add_trace(data = eic, x = ~rt, y = ~EIC_int, color = ~spindir, opacity = 0.5,
                      mode="lines", type="scatter",
                      colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone"))) %>%
            add_trace(data = TIC_df, x=~rt, y=~TIC, 
                      mode="lines", type="scatter", line=list(color="black"),
                      hoverinfo="none") %>%
            layout(xaxis = list(title = "Retention time (s)"),
                   yaxis = list(title = "Intensity",
                                fixedrange = TRUE))
    } else {
        plot_ly(data = eic, x = ~rt, y = ~EIC_int, color = ~spindir, alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone")),
                source = "EIC") %>%
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
    plot_ly(data=scandata, x=~mz, y=~TIS, type = "bar", source = "TIS") %>%
        layout(xaxis = list(title = "m/z"),
               yaxis = list(title = "Intensity",
                            fixedrange = TRUE)) %>%
        event_register("plotly_click")
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
                         value = 5)
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
        TIS_data <- event_data(event = "plotly_click", source="TIS")
        if(is.null(TIS_data)){
            plotGivenEIC(input$given_mz, ppm=input$given_ppm)
        } else {
            plotGivenEIC(TIS_data$x, ppm=input$given_ppm)
        }
    })
    
    output$debug <- renderPrint({
        #EIC_data <- event_data(event = "plotly_click", source = "EIC")
        #print(dev.list())
        #print(dev.cur())
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
