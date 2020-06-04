# Code to read MS data out of mzML files and into a database/dataframe format

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
  cbind(fileid=basename(filename), all_data)
}

# Metadata ----
sample_files <- list.files("../Data/", pattern = "\\.mzML", full.names = TRUE)
sample_files <- normalizePath(sample_files)

metadframe <- data.frame(
  fileid=basename(sample_files),
  resolution=regmatches(regexpr(pattern = "_1|_2", sample_files), x = sample_files)
)
metadframe$resolution <- c(`_1`="lowres", `_2`="highres")[metadframe$resolution]
write.csv(metadframe, "../Data/metadata.csv", row.names = FALSE)

# Grab the actual data and clean up a little ----
raw_data <- lapply(sample_files, grabSingleFileData)
raw_data <- do.call(rbind, raw_data)
saveRDS(raw_data, file = "../Data/MS1_data_frame.rds")

# Grab standards file and clean up a little ----
stans <- read.csv(paste0("https://raw.githubusercontent.com/kheal/Example_Unta",
                         "rgeted_Metabolomics_Workflow/master/Ingalls_Lab_Stan",
                         "dards.csv"))
clean_stans <- subset(stans, stans$Column=="HILIC")
clean_stans <- clean_stans[c("Compound.Type", "Compound.Name", "Emperical.Formula",
                             "RT..min.", "m.z", "ionization_form", "Fraction1")]
write.csv(x = clean_stans, file = "../Data/stans.csv")
