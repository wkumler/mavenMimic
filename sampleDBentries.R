# Database examples
library(RSQLite)
library(dplyr)
library(plotly)
library(ggplot2)

# Functions ----


pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

plotGivenEIC <- function(mass, ppm=5, db="falkor.db", plotTIC=FALSE){
  db_conn <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
  df <- dbGetQuery(db_conn, "SELECT * FROM raw_data_smol WHERE mz>? AND mz<?", 
                       params=pmppm(mass, ppm = ppm))
  dbDisconnect(db_conn)
  eic <- df %>% group_by(file, rt) %>% summarize(EIC_int=sum(int))
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
                mode="lines+markers", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity",
                          fixedrange = TRUE))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~EIC_int, color = ~spindir, alpha = 0.5,
            mode="lines+markers", type="scatter",
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
falkorDb <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
dbListTables(falkorDb)
bet_df <- dbGetQuery(falkorDb, "SELECT * FROM raw_data_smol WHERE mz>? AND mz<?", 
                     params=pmppm(117.078979+1.007276))
dbDisconnect(falkorDb)

plotGivenEIC(117.078979+1.007276)

db_conn <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
TIC <- dbGetQuery(db_conn, "SELECT rt,mz,SUM(int) FROM raw_data_smol GROUP BY rt;")
dbDisconnect(db_conn)
