# NOTE: This code is no longer operational with the app because
# the new app has been reactivitized. Consider instead adding
# a checkbox to pull from database rather than memory and re-write
# the reactive() functions in the app.

# Code to create a database from a dataframe ----
library(RSQLite)
if(file.exists("falkor.db"))file.remove("falkor.db")
if(!file.exists("falkor.db"))file.create("falkor.db")

falkor_db <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
dbWriteTable(conn = falkor_db, name = "raw_data", value = raw_data, overwrite=TRUE)

# Create TIC
tic_query <- "SELECT fileid,rt,mz,SUM(int) FROM raw_data GROUP BY rt;"
dbWriteTable(falkor_db, "TIC", dbGetQuery(falkor_db, tic_query))


# Check on data export and disconnect
dbListTables(falkor_db)
dbDisconnect(falkor_db)



# Code to extract and plot data from a dataframe ----

plotGivenEIC_DB <- function(mass, ppm=5, db="falkor.db", plotTIC=TRUE,
                            ms_data_dir="G:/My Drive/FalkorFactor/mzMLs",
                            mdframe=metadframe, plotby="depth"){
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
      add_trace(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green"), unique(eic[["plotby"]]))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", line=list(color="black"),
                hoverinfo="none") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter",
            colors = setNames(c("red", "blue", "green"), unique(eic[["plotby"]])),
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