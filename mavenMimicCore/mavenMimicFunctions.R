library(shiny)
library(plotly)
library(dplyr)

# Read in raw_data
load("raw_data_frame")

# Functions
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
    plot_ly() %>%
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
  plot_ly(data=scandata, x=~mz, y=~TIS, type = "bar") %>%
    layout(xaxis = list(title = "m/z"),
           yaxis = list(title = "Intensity",
                        fixedrange = TRUE))
}