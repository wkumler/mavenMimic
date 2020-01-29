# Code to read MS data out of mzML files and into a more friendly, long-form csv

library(pbapply)

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

raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data_frame <- as.data.frame(do.call(rbind, raw_data))

save(raw_data_frame, file = "raw_data_frame")
