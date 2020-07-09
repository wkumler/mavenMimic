
# Setup things ----

library(shiny)
library(dplyr)
library(plotly)
library(data.table)

addisos <- c(H=1.007276, Na=22.98977, K=38.963158, NH4=18.033823,
             C13=1.003355, N15=0.997035, O18=2.004244, S34=1.995796)

cat("Reading in raw data... ")
MS1_data_frame <- as.data.table(readRDS("Data/MS1_data_frame"))
cat("Done\n")
cat("Reading in MSMS data... ")
MS2_data_frame <- as.data.table(readRDS("Data/MS2_data_frame"))
cat("Done\n")
cat("Creating TIC... ")
tic <- MS1_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))
cat("Done\n")
cat("Reading in metadata... ")
falkor_metadata <- read.csv("Data/falkor_metadata.csv")
cat("Done\n")
cat("Reading in standards list... ")
stans_csv <- read.csv("Data/falkor_stans.csv")
stans_namelist <- split(stans_csv, stans_csv$Fraction1) %>%
  lapply(FUN = `[[`, "Compound.Name") %>%
  lapply(FUN = sort) %>%
  `names<-`(c("Negative mode", "Positive mode")) %>%
  `[`(c(2,1))
cat("Done\n")
browseURL("http://127.0.0.1:1337/")


# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

get_EIC <- function(MS1_data_frame, mass, ppm=5, mdframe=falkor_metadata){
  req(mass)
  MS1_data_frame[mz%between%pmppm(mass, ppm = ppm)] %>%
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
}

get_Spectrum <- function(MS1_data_frame, scan, ret_window=1){
  req(scan)
  MS1_data_frame[rt%between%c(scan-ret_window/2, scan+ret_window/2)] %>%
    #mutate(mz=round(mz*1000)/1000) %>% # Round to group nearby scans
    group_by(mz) %>% 
    summarize(TIS=sum(int))
}

plotGivenEIC <- function(eic, plotby="depth", plottic=TRUE, tic=NULL,
                         current_mass=118.0865, ppm=5){
  if(plottic){
    tic$int <- (tic$int/max(tic$int))*max(eic$int)
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, text=~fileid,
                color = ~get(plotby), alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green", "black"), 
                                  unique(eic[[plotby]]))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", 
                line=list(color="black"),
                name="TIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = list(text=paste(round(pmppm(current_mass, ppm = ppm), digits = 4), 
                                     collapse = " - "),
                          x=0.07))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, text=~fileid,
            color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("red", "blue", "green"), 
                              unique(eic[[plotby]]))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = list(text=paste(round(pmppm(current_mass, ppm = ppm), digits = 4), 
                                     collapse = " - "),
                          x=0.07))
  }
}

plotGivenSpectrum <- function(spectrum){
  plot_ly(source = "TIS") %>%
    add_trace(data=spectrum, x=~mz, y=~TIS, type = "scatter", 
              mode="lines", line=list(color="black")) %>%
    layout(xaxis = list(title = "m/z"),
           yaxis = list(title = "Intensity"),
           title = list(
             text="TIS",
             x = 0.1
           ))
}

plotMSMS <- function(mass, ret_time = 1, ppm=5, ret_win=20, dataframe=MS2_data_frame){
  frags <- dataframe %>% 
    filter(premz%between%pmppm(mass, ppm = ppm)) %>%
    filter(rt>ret_time-ret_win&rt<ret_time+ret_win) %>%
    arrange(fragmz)
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
                mode="lines+markers", type="scatter", 
                marker=list(color="#000000FF"),
                line=list(color="#000000FF")) %>%
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

data_MSMS <- function(mass, ret_time = 1, ppm=5, ret_win=20, dataframe=MS2_data_frame){
  dataframe %>% 
    filter(premz>min(pmppm(mass, ppm = ppm))) %>%
    filter(premz<max(pmppm(mass, ppm = ppm))) %>%
    filter(rt>ret_time-ret_win&rt<ret_time+ret_win) %>%
    arrange(fragmz)
}

makeAddTable <- function(named_vec){
  header <- paste0("<th>", names(named_vec), "</th>", collapse = "")
  header <- paste0("<tr>", header, "</tr>")
  addrow <- paste0("<a href='#' onclick='detect_click(this)'>+", named_vec)
  addrow <- paste0("<td>", addrow, "</td>", collapse = "")
  addrow <- paste0("<tr>", addrow, "</tr>")
  rmrow <- paste0("<a href='#' onclick='detect_click(this)'>", named_vec*-1)
  rmrow <- paste0("<td>", rmrow, "</td>", collapse = "")
  rmrow <- paste0("<tr>", rmrow, "</tr>")
  paste0('<table style="width:100%">', 
         header, addrow, rmrow,
         "</table>")
}



# UI ----
ui <- navbarPage(
  title = "Maven Mimic",
  tabPanel(title = "Exploration",
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
                             value = TRUE),
               selectInput(inputId = "given_stan",
                           label = "Choose a standard:",
                           choices = stans_namelist,
                           #selectize = TRUE,
                           selected = "Betaine")
             ),
             
             mainPanel(
               htmlOutput("adduct_table"),
               h3(),
               plotlyOutput(outputId = "chrom", height = "300px"),
               h3(),
               fluidRow(
                 splitLayout(
                   cellWidths = c("49%", "49%"),
                   plotlyOutput(outputId = "TIS", height = "250px"),
                   plotlyOutput(outputId = "MSMS", height = "250px")
                 )
               ),
               includeScript("detect_click.js")
             )
           )
         ),
  tabPanel(
    title = "Export",
    verticalLayout(
      h3("Chromatogram data"),
      verbatimTextOutput("export_chrom"),
      h3("TIS data"),
      verbatimTextOutput("export_spectrum"),
      h3("MS2 data"),
      verbatimTextOutput("export_MSMS")
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
  observeEvent(input$given_stan, {
    stan_mass <- stans_csv$m.z[stans_csv$Compound.Name==input$given_stan]
    current_mass(as.numeric(as.character(stan_mass)))
  })
  observeEvent(input$clicked_text, {
    current_mass(current_mass()+as.numeric(input$clicked_text))
  })
  observeEvent(current_mass(), {
    updateNumericInput(session, "given_mz", value = current_mass())
  })
  
  
  
  given_EIC <- reactive({
    get_EIC(MS1_data_frame = MS1_data_frame, 
            mass = current_mass(), 
            ppm = input$given_ppm, 
            mdframe = falkor_metadata)
  })
  
  given_Spectrum <- reactive({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    get_Spectrum(MS1_data_frame = MS1_data_frame, scan=EIC_data$x)
  })
  
  
  
  output$adduct_table <- renderUI(HTML(makeAddTable(addisos)))
  
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
             dataframe = MS2_data_frame)
  })
  
  output$export_chrom <- renderPrint({
    dput(get_EIC(MS1_data_frame = MS1_data_frame, 
                 mass = current_mass(), 
                 ppm = input$given_ppm, 
                 mdframe = falkor_metadata))
  })
  output$export_spectrum <- renderPrint({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    dput(get_Spectrum(MS1_data_frame = MS1_data_frame, scan=EIC_data$x))
  })
  output$export_MSMS <- renderPrint({
    EIC_data <- event_data(event = "plotly_click", source = "EIC")
    req(EIC_data)
    dput(data_MSMS(mass = current_mass(), ret_time = EIC_data$x, 
              ppm = input$given_ppm, ret_win = 20,
              dataframe = MS2_data_frame))
  })
}

# Run the application ----
shinyApp(ui = ui, server = server, options = list(port=1337))
