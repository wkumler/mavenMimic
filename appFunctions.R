
# Libraries required ----
library(dplyr)
library(plotly)
library(RSQLite)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

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
plotMSMS(214.0905, ret_time = 1300)
