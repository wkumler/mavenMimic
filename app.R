
# Setup things ----

library(shiny)
library(dplyr)
library(plotly)

# cat("Reading in raw data... ")
# alldata <- readRDS("raw_data_table")
# cat("Done\n")
# cat("Reading in MSMS data... ")
# raw_msms_data <- readRDS("raw_msms_table")
# cat("Done\n")
# cat("Creating TIC... ")
# tic <- alldata %>% mutate(rt=round(rt)) %>%
#   group_by(rt) %>% summarize(int=sum(int))
# cat("Done\n")
# cat("Reading in metadata... ")
# falkor_metadata <- read.csv("falkor_metadata.csv")
# cat("Done\n")

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

get_EIC <- function(alldata, mass, ppm=5, 
                    mdframe=falkor_metadata){
  alldata %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
}

get_Spectrum <- function(alldata, scan, ret_window=1){
  alldata %>% 
    filter(rt>scan-ret_window/2&rt<scan+ret_window/2) %>% 
    mutate(mz=round(mz*1000)/1000) %>% # Round to group nearby scans
    group_by(mz) %>% 
    summarize(TIS=sum(int))
}

plotGivenEIC <- function(eic, plotby="depth", plottic=TRUE, tic=NULL,
                         current_mass=118.0865, ppm=5){
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

plotMSMS <- function(mass, ret_time = 1, ppm=5, ret_win=20, dataframe=raw_msms_data){
  frags <- dataframe %>% 
    filter(premz>min(pmppm(mass, ppm = ppm))) %>%
    filter(premz<max(pmppm(mass, ppm = ppm))) %>%
    filter(rt>ret_time-ret_win&rt<ret_time+ret_win)
  if(!nrow(frags)){
    empty_plot_ly <- plot_ly(x=ret_time, y=mass, 
                             text="No data for any voltage here",
                             mode="text", type="scatter") %>%
      layout(title = paste(round(pmppm(mass, ppm), digits = 4), 
                           collapse = " - "),
             xaxis = list(range=c(ret_time-ret_win, ret_time+ret_win)),
             yaxis = list(range=pmppm(mass, ppm)))
    return(empty_plot_ly)
  }
  pls <- lapply(split(frags, frags$nrg), function(given_frags){
    frags2plot <- given_frags %>%
      filter(int>max(int/10000)) %>%
      mutate(int=(int/max(int))*100)
    frag_lines <- lapply(seq_len(nrow(frags2plot)), function(x){
      line_obj <- list(type="line", xref="x", yref="y")
      line_obj$x0 <- frags2plot[x, "fragmz"]
      line_obj$x1 <- frags2plot[x, "fragmz"]
      line_obj$y0 <- 0
      line_obj$y1 <- frags2plot[x, "int"]
      return(line_obj)
    })
    pl <- plot_ly(source = "MSMS") %>%
      add_trace(data = frags2plot, x=~fragmz, y=~int,
                mode="markers", type="scatter", marker=list(color="#00000000")) %>%
      layout(xaxis = list(title = "m/z", range=c(0, max(frags2plot$fragmz)+5)),
             yaxis = list(title = "Intensity", range=c(0, 120)),
             shapes = frag_lines,
             title=paste0("Fragments for ", 
                          round(frags2plot[1, "premz"], digits = 4), 
                          " Da @ ~",
                          round(mean(frags2plot$rt)), "seconds"),
             showlegend = FALSE) %>%
      add_annotations(
        text = paste("Voltage:", unique(frags2plot$nrg)),
        x = 0.5,
        y = 1,
        yref = "paper",
        xref = "paper",
        xanchor = "middle",
        yanchor = "top",
        showarrow = FALSE,
        font = list(size = 15)
      )
    return(pl)
  })
  subplot(pls, shareX = TRUE, nrows = length(pls))
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
                          value = FALSE)
        ),

        mainPanel(
          h2("Extracted ion chromatogram: "),
          plotlyOutput(outputId = "chrom", height = "300px"),
          h2("Ions identified in chosen scan: "),
          plotlyOutput(outputId = "TIS", height = "300px"),
          h2("Fragments of the selected ion: "),
          plotlyOutput(outputId = "MSMS", height = "300px")
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
            mdframe = falkor_metadata)
  })
  
  given_Spectrum <- reactive({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    get_Spectrum(alldata = alldata, scan=EIC_data$x)
  })
  
  output$chrom <- renderPlotly({
    plotGivenEIC(eic = given_EIC(), 
                 plotby = input$treatment,
                 plottic = input$user_plottic, 
                 tic = tic,
                 current_mass = current_mass(),
                 ppm = input$given_ppm)
  })

  output$TIS <- renderPlotly({
    plotGivenSpectrum(given_Spectrum())
  })
  
  output$MSMS <- renderPlotly({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    plotMSMS(mass = current_mass(), ret_time = EIC_data$x, 
             ppm = input$given_ppm, ret_win = 20,
             dataframe = raw_msms_data)
  })
}

  
# Run the application ----
shinyApp(ui = ui, server = server)
