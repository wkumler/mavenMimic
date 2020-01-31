# Database examples
# Setup things ----
library(RSQLite)
library(dplyr)
library(plotly)
library(ggplot2)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

plotGivenEIC <- function(mass, ppm=5, db="falkor.db", plotTIC=FALSE,
                         ms_data_dir="G:/My Drive/FalkorFactor/mzMLs"){
  falkor_db <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
  df <- dbGetQuery(falkor_db, "SELECT * FROM raw_data WHERE mz>? AND mz<?", 
                       params=pmppm(mass, ppm = ppm))
  eic <- df %>% group_by(fileid, rt) %>% summarize(int=sum(int))
  eic$names <- list.files(ms_data_dir, pattern = "Smp|Blk")[eic$fileid]
  eic$depth <- c("Cyclone", "Anticyclone")[1-grepl("62|64", eic$names)+1]
  
  tic <- dbGetQuery(falkor_db, "SELECT * FROM tic") %>% 
    mutate(int=`SUM(int)`) %>% select(-`SUM(int)`) %>% 
    mutate(rt=round(rt)) %>% group_by(rt) %>% summarise(int=sum(int))
  tic$int <- (tic$int/max(tic$int))*max(eic$int)
  dbDisconnect(falkor_db)
  if(plotTIC){
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~depth, opacity = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone"))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~EIC_int, color = ~spindir, alpha = 0.5,
            mode="lines+markers", type="scatter",
            colors = setNames(c("blue", "green"), c("Cyclone", "Anticyclone")),
            source = "EIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  }
}

plotGivenScan <- function(ret, window=1, df=raw_data_frame){
  scandata <- df %>% 
    filter(rt>ret-window/2&rt<ret+window/2) %>% 
    mutate(mz=round(mz*1000)/1000) %>%
    group_by(mz) %>% 
    summarize(TIS=sum(int))
  plot_ly(source = "TIS") %>%
    add_trace(data=scandata, x=~mz, y=~TIS, 
              type = "bar", marker=list(color="black"),
              hoverinfo="none") %>%
    add_trace(data=scandata, x=~mz, y=~TIS, 
              type="scatter", mode="markers", 
              marker=list(color="black")) %>%
    layout(xaxis = list(title = "m/z"),
           yaxis = list(title = "Intensity",
                        fixedrange = TRUE)) %>%
    event_register("plotly_click")
}



# Applications ----
plotGivenEIC(mass = 118.086804, plotTIC = FALSE)