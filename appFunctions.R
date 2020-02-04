
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
  spindir=c("Blank", "Cyclone", "Anticyclone")[c(1, (1-ceiling(1:24/12)%%2)+2)],
  time=c("Blank", "Morning", "Afternoon")[c(1, ceiling(1:24/6)%%2+2)]
  )

# Raw file functions ----
cat("Reading in raw data... ")
raw_data_frame <- readRDS("raw_data_table")
cat("Done\n")
cat("Reading in MSMS data... ")
raw_msms_data <- readRDS("raw_msms_table")
cat("Done\n")
cat("Creating TIC... ")
tic <- raw_data_frame %>% mutate(rt=round(rt)) %>%
  group_by(rt) %>% summarize(int=sum(int))
cat("Done")

plotGivenEIC <- function(mass=118.0865, ppm=5, df=raw_data_frame, 
                         mdframe = metadframe, plottic=TRUE,
                         plotby="depth"){
  eic <- df %>% 
    filter(mz>min(pmppm(mass, ppm = ppm))&mz<max(pmppm(mass, ppm = ppm))) %>% 
    group_by(fileid, rt) %>% 
    summarize(int=sum(int)) %>%
    left_join(mdframe, by="fileid")
  if(plottic){
    tic$int <- (tic$int/max(tic$int))*max(eic$int)
    plot_ly(source = "EIC") %>%
      add_trace(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
                mode="lines", type="scatter",
                colors = setNames(c("red", "blue", "green"), unique(eic[[plotby]]))) %>%
      add_trace(data = tic, x=~rt, y=~int, 
                mode="lines", type="scatter", 
                line=list(color="black"),
                name="TIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = paste("Mass range:", min(pmppm(mass, ppm = ppm)),
                           "-", max(pmppm(mass, ppm = ppm))))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter", source="EIC",
            colors = setNames(c("red", "blue", "green"), unique(eic[[plotby]]))) %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = paste("Mass range:", min(pmppm(mass, ppm = ppm)),
                           "-", max(pmppm(mass, ppm = ppm))))
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
             yaxis = list(title = "Intensity"),
             title = paste("Mass range:", min(pmppm(mass, ppm = ppm)),
                           "-", max(pmppm(mass, ppm = ppm))))
  } else {
    plot_ly(data = eic, x = ~rt, y = ~int, color = ~get(plotby), alpha = 0.5,
            mode="lines", type="scatter",
            colors = setNames(c("red", "blue", "green"), unique(eic[["plotby"]])),
            source = "EIC") %>%
      layout(xaxis = list(title = "Retention time (s)"),
             yaxis = list(title = "Intensity"),
             title = paste("Mass range:", min(pmppm(mass, ppm = ppm)),
                           "-", max(pmppm(mass, ppm = ppm))))
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

# plotGivenEIC(118.0865)
# plotGivenEIC(118.0865, plotby = "spindir")
# plotGivenScan(345)
# plotGivenEIC_DB(118.0865)
# plotGivenScan_DB(345)

# MSMS functions ----
plotMSMS <- function(mass, ppm=5, retwin=50, dataframe=raw_msms_data){
  frags <- dataframe %>% 
    filter(premz>min(pmppm(mass, ppm = ppm))&premz<max(pmppm(mass, ppm = ppm)))
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
             yaxis = list(title = "Intensity"),
             shapes = frag_lines,
             title=paste0("Fragments for ", 
                          round(frags2plot[1, "premz"], digits = 4), 
                          " Da @ ~",
                         round(mean(frags2plot$rt)), "seconds"))
    return(pl)
  })
  subplot(pls, shareX = TRUE, nrows = length(pls))
}
# plotMSMS(132.1019)
# 
# frags20 <- frags[frags$nrg==50,]
# split_msms <- split(frags20, frags20$rt)
# for(i in seq_along(split_msms)){
#   given_frags <- split_msms[[i]]
#   print(plot(given_frags$fragmz, given_frags$int/max(given_frags$int), 
#              type="h", xlim=c(50, 140), ylim = c(0,1), main = paste0(
#                "Retention time: ", round(unique(given_frags$rt), digits = 4), 
#                " Pre m/z: ", round(unique(given_frags$premz), digits = 4)
#              )))
# }
# 
# v <- cbind(frags20$fragmz[frags20$rt>250&frags20$rt<350], 
#            frags20$int[frags20$rt>250&frags20$rt<350])
# apply(v, 1, function(x){cat(x); cat("\n")})
