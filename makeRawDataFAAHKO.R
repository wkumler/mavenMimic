# Code to read MS data out of mzML files and into a database/dataframe format

# Setup things ----
library(pbapply)

if(!dir.exists("Data")){
  dir.create("Data")
}

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

grabSingleFileMS2 <- function(filename){
  msdata <- mzR::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  ms2rows <- seq_len(nrow(fullhd))[fullhd$msLevel>1]
  spectra_list <- lapply(ms2rows, function(x){
    rtime <- fullhd[x, "retentionTime"]
    premz <- fullhd[x, "precursorMZ"]
    fragments <- mzR::peaks(msdata, x)
    return(cbind(rtime, premz, fragments))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "premz", "fragmz", "int"))
  return(all_data)
}

# Metadata ----
sample_files <- list.files("C:/Users/willi/Documents/R/win-library/3.6/faahKO/cdf/", 
                           full.names = TRUE, recursive = TRUE)


metadframe <- data.frame(
  fileid=1:12, 
  filenames=basename(sample_files),
  treatment=c(rep("KO", 6), rep("WT", 6))
)
write.csv(metadframe, "Data/faahko_metadata.csv", row.names = FALSE)

# Grab the actual data and clean up a little ----
raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data <- lapply(seq_along(raw_data), function(x){
  cbind(fileid=x, raw_data[[x]])
})
raw_data <- do.call(rbind, raw_data)
saveRDS(raw_data, file = "Data/faahko_MS1_data_frame")

# # Grab MSMS data and clean up a little ----
# msms_data_dir <- "mzMLs/MSMS"
# msms_files <- list.files(msms_data_dir, pattern = "DDApos", full.names = TRUE)
# raw_msmsdata <- pblapply(msms_files, grabSingleFileMS2)
# nrgs <- as.numeric(gsub(".*neg|.*pos|\\.mzML", "", msms_files))
# raw_msmsdata <- lapply(seq_along(raw_msmsdata), function(x){
#   cbind(nrg=nrgs[x], raw_msmsdata[[x]])
# })
# raw_msmsdata <- do.call(rbind, raw_msmsdata)
# saveRDS(raw_msmsdata, file = "Data/faahko_MS2_data_frame")