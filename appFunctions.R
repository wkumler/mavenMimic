
# Libraries required ----
library(dplyr)
library(plotly)
library(RSQLite)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

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
# plotMSMS(245.11)
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
