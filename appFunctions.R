
# Libraries required ----
library(dplyr)
library(plotly)
library(RSQLite)

# Functions and metadata ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
metadframe <- data.frame(
  fileid=1:25, 
  filenames=list.files("G:/My Drive/FalkorFactor/mzMLs", pattern = "Smp|Blk"),
  depth=c("Blank", "DCM", "25m")[c(1, ceiling(1:24/3)%%2+2)],
  spindir=c("Blank", "Cyclone", "Anticyclone")[c(1, (1-ceiling(1:24/12)%%2)+2)]
  )

# Raw file functions ----
raw_data_frame <- readRDS("raw_data_table")
tic <- raw_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))

plotGivenEIC <- function(mass, ppm=5, df=raw_data_frame, 
                         mdframe = metadframe, plottic=TRUE){
  eic <- df %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
  if(plottic){
    tic$int <- (tic$int/max(tic$int))*max(eic$int)
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~depth, alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green"), unique(eic$depth))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~depth, alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("red", "blue", "green"), unique(depth))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
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


# Database functions ----
plotGivenEIC_DB <- function(mass, ppm=5, db="falkor.db", plotTIC=TRUE,
                            ms_data_dir="G:/My Drive/FalkorFactor/mzMLs"){
  falkor_db <- dbConnect(drv = RSQLite::SQLite(), db)
  df <- dbGetQuery(falkor_db, "SELECT * FROM raw_data WHERE mz>? AND mz<?", 
                   params=pmppm(mass, ppm = ppm))
  eic <- df %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
  
  tic <- dbGetQuery(falkor_db, "SELECT * FROM tic") %>% 
    mutate(int=`SUM(int)`) %>% select(-`SUM(int)`) %>% 
    mutate(rt=round(rt)) %>% group_by(rt) %>% summarise(int=sum(int))
  tic$int <- (tic$int/max(tic$int))*max(eic$int)
  dbDisconnect(falkor_db)
  if(plotTIC){
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~depth, alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green"), unique(eic$depth))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~depth, alpha = 0.5,
            mode="lines", type="scatter",
            colors = setNames(c("red", "blue", "green"), unique(eic$depth)),
            source = "EIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  }
}

plotGivenScan_DB <- function(ret, window=1, db="falkor.db"){
  falkor_db <- dbConnect(drv = RSQLite::SQLite(), db)
  df <- dbGetQuery(falkor_db, "SELECT * FROM raw_data WHERE rt>? AND rt<?", 
                   params=c(ret-window, ret+window))
  dbDisconnect(falkor_db)
  tis <- df %>% group_by(mz=round(mz, digits = 2)) %>% summarize(int=sum(int))
  
  plot_ly(source = "TIS") %>%
    add_trace(data=tis, x=~mz, y=~int, 
              type = "bar", marker=list(color="black"),
              hoverinfo="none") %>%
    add_trace(data=tis, x=~mz, y=~int, 
              type="scatter", mode="markers", 
              marker=list(color="black")) %>%
    layout(xaxis = list(title = "m/z"),
           yaxis = list(title = "Intensity",
                        fixedrange = TRUE)) %>%
    event_register("plotly_click")
}

plotGivenEIC(118.0865)
plotGivenScan(345)
plotGivenEIC_DB(118.0865)
plotGivenScan_DB(345)