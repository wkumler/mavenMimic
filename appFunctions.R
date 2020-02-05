
# Libraries required ----
library(dplyr)
library(plotly)
library(RSQLite)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}


plotMSMS(214.0905, ret_time = 1300)
