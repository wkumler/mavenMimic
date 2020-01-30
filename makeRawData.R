# Code to read MS data out of mzML files and into a more friendly, long-form csv

library(pbapply)
library(RSQLite)

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

ms_data_dir <- "G:/My Drive/FalkorFactor/mzMLs"
sample_files <- list.files(ms_data_dir, pattern = "Smp", full.names = TRUE)
sample_spins <- c("Cyclone", "Anticyclone")[1-grepl("62|64", sample_files)+1]
sample_depth <- c("25m", "DCM")[grepl("DCM", sample_files)+1]

raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data <- lapply(seq_along(raw_data), function(x){
  cbind(fileid=x, filename=basename(sample_files[x]), 
        depth=sample_depth[x], spindir=sample_spins[x], raw_data[[x]])
})
raw_data_frame <- as.data.frame(do.call(rbind, raw_data))

falkorDb <- dbConnect(drv = RSQLite::SQLite(), "falkor.db")
dbWriteTable(falkorDb, "raw_data", raw_data_frame, overwrite=TRUE)
dbListTables(falkorDb)
dbDisconnect(falkorDb)

dbGetQuery(falkorDb, "SELECT file,rt FROM raw_data_smol WHERE rt<121")
