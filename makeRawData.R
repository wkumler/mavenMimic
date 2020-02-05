# Code to read MS data out of mzML files and into a database/dataframe format

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
ms_data_dir <- "G:/My Drive/FalkorFactor/mzMLs"
sample_files <- list.files(ms_data_dir, pattern = "Smp|Blk", full.names = TRUE)

metadframe <- data.frame(
  fileid=1:25, 
  filenames=list.files("G:/My Drive/FalkorFactor/mzMLs", pattern = "Smp|Blk"),
  depth=c("Blank", "DCM", "25m")[c(1, ceiling(1:24/3)%%2+2)],
  spindir=c("Blank", "Cyclone", "Anticyclone")[c(1, (1-ceiling(1:24/12)%%2)+2)],
  time=c("Blank", "Morning", "Afternoon")[c(1, ceiling(1:24/6)%%2+2)]
)
write.csv(metadframe, "falkor_metadata.csv", row.names = FALSE)

# Grab the actual data and clean up a little ----
raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data <- lapply(seq_along(raw_data), function(x){
  cbind(fileid=x, raw_data[[x]])
})
raw_data <- do.call(rbind, raw_data)
raw_data <- raw_data[raw_data$rt>60&raw_data$rt<1100,]
saveRDS(raw_data, file = "raw_data_table")

# Grab MSMS data and clean up a little ----
msms_data_dir <- "G:/My Drive/FalkorFactor/mzMLs/MSMS"
msms_files <- list.files(msms_data_dir, pattern = "DDApos", full.names = TRUE)
raw_msmsdata <- pblapply(msms_files, grabSingleFileMS2)
nrgs <- gsub(".*neg|.*pos|\\.mzML", "", msms_files)
raw_msmsdata <- lapply(seq_along(raw_msmsdata), function(x){
  cbind(nrg=nrgs[x], raw_msmsdata[[x]])
})
raw_msmsdata <- do.call(rbind, raw_msmsdata)
saveRDS(raw_msmsdata, file = "raw_msms_table")

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