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
sample_files <- normalizePath(list.files("falkor_mzMLs", pattern = "Smp|Blk|Std", 
                                         full.names = TRUE))

metadframe <- data.frame(
  fileid=basename(sample_files),
  depth=c("Blank", "DCM", "25m", "Std")[c(1, ceiling(1:24/3)%%2+2)],
  spindir=c("Blank", "Cyclone", "Anticyclone", "Std")[c(1, (1-ceiling(1:24/12)%%2)+2)],
  time=c("Blank", "Morning", "Afternoon", "Std")[c(1, ceiling(1:24/6)%%2+2)]
)
write.csv(metadframe, "Data/falkor_metadata.csv", row.names = FALSE)

# Grab the actual data and clean up a little ----
raw_data <- pblapply(sample_files, grabSingleFileData)
raw_data <- pblapply(seq_along(raw_data), function(x){
  cbind(fileid=basename(sample_files)[x], raw_data[[x]])
})
raw_data <- do.call(rbind, raw_data)
raw_data <- raw_data[raw_data$rt>60&raw_data$rt<1100,]
saveRDS(raw_data, file = "Data/MS1_data_frame")

# Grab MSMS data and clean up a little ----
msms_data_dir <- "falkor_mzMLs/MSMS"
msms_files <- list.files(msms_data_dir, pattern = "DDApos", full.names = TRUE)
raw_msmsdata <- pblapply(msms_files, grabSingleFileMS2)
nrgs <- as.numeric(gsub(".*neg|.*pos|\\.mzML", "", msms_files))
raw_msmsdata <- lapply(seq_along(raw_msmsdata), function(x){
  cbind(nrg=nrgs[x], raw_msmsdata[[x]])
})
raw_msmsdata <- do.call(rbind, raw_msmsdata)
saveRDS(raw_msmsdata, file = "Data/MS2_data_frame")
# Grab standards file and clean up a little ----
stans <- read.csv(paste0("https://raw.githubusercontent.com/kheal/Example_Unta",
                         "rgeted_Metabolomics_Workflow/master/Ingalls_Lab_Stan",
                         "dards.csv"))
clean_stans <- subset(stans, stans$Column=="HILIC")
clean_stans <- clean_stans[c("Compound.Type", "Compound.Name", "Emperical.Formula",
                             "RT..min.", "m.z", "ionization_form", "Fraction1")]
write.csv(x = clean_stans, file = "Data/falkor_stans.csv")
