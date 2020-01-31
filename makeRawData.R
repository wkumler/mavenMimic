# Code to read MS data out of mzML files and into a more database format

# Setup things ----
library(pbapply)
library(RSQLite)
init <- Sys.time()

# Functions ----
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}

# Metadata ----
ms_data_dir <- "G:/My Drive/FalkorFactor/mzMLs"
sample_files <- list.files(ms_data_dir, pattern = "Smp|Blk", full.names = TRUE)
sample_spins <- c("Cyclone", "Anticyclone")[1-grepl("62|64", sample_files)+1]
sample_depth <- c("25m", "DCM")[grepl("DCM", sample_files)+1]

# Grab the actual data and clean up a little ----
raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data <- lapply(seq_along(raw_data), function(x){
  cbind(fileid=x, raw_data[[x]])
})
raw_data <- do.call(rbind, raw_data)

# Connect to database and write out data ----
# NOTE: Make sure you've manually created falkor.db first
if(file.exists("falkor.db"))file.remove("falkor.db")
if(!file.exists("falkor.db"))file.create("falkor.db")

falkor_db <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
dbWriteTable(conn = falkor_db, name = "raw_data", value = raw_data, overwrite=TRUE)

# Create TIC ----
tic_query <- "SELECT fileid,rt,mz,SUM(int) FROM raw_data GROUP BY rt;"
dbWriteTable(falkor_db, "TIC", dbGetQuery(falkor_db, tic_query))


# Check on data export and disconnect ----
dbListTables(falkor_db)
dbDisconnect(falkor_db)

print(Sys.time()-init)