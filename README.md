# mavenMimic Introduction

I really like how interactive and exploratory [MAVEN](http://genomics-pubs.princeton.edu/mzroll/index.php) is, with its (relatively) intuitive interaction between chromatograms (single mass, many scans) and spectra (single scan, many masses). Combined with a TIC, it's easy to poke quickly at the data and get an initial sense of interesting compounds.

However, MAVEN also had a few drawbacks. For one, it's no longer being developed and is technically defunct. Still works great, but has a few shortcomings. Specifically, MAVEN only accepts the old .mzXML file format, doesn't handle MS<sup>2+</sup> data, and is hard (impossible?) to color-code by treatment condition.

This app mimics MAVEN's functionality while resolving the above issues. It's built on the `shiny` and `plotly` libraries in R and critically exploits the `event_data()` function to pass information back and forth between the chromatograms and the spectra. The app leverages my own code snippets for getting MS data out of weird file types (thanks to the `mzR` package) and into a simple data frame in R, which can then be `filter`ed, `mutate`d, and generally easily manipulated.

Some considerations! This app is nowhere near as fast as MAVEN - it often takes a half-second or so to render the plots after a new interaction takes place. It's also quite memory-intensive, as the 25 HILIC positive files I've been looking at from the R/V Falkor cruise occupy about a gigabyte of RAM as a `data frame` object in the R environment. One way to get around this is to use SQL and the `RSQLite` package to create a database and handle the filtering. This code exists as "oldMavenDBcode.R" but the I/O times made it almost unusable.

## Making the data

Hopefully you've already got some mass-spec data to play with. If so, convert it to open-source formats such as .mzML and modify `makeRawData.R` to point to it. (i.e. modify the lines below to output the files you want to process).

```r
ms_data_dir <- "G:/My Drive/FalkorFactor/mzMLs"
sample_files <- list.files(ms_data_dir, pattern = "Smp|Blk", full.names = TRUE)
```

There's an additional metadata frame to add treatment information, so if you're
sampling different conditions you'll want to append columns to that to include
the relevant information. Metadata is attached to the raw data after it's been
filtered to reduce sample size, and is done via `dplyr`'s `left_join` function, merging by fileid number. For the Falkor cruise data that I'm working on, my metadata dataframe looks like this:

| fileid  | filenames | depth | spindir | time |
| ------- | --------- | ----- | ------- | ---- |
| 1 | 190715_Blk_KM1906U14-Blk_C.mzML | Blank | Blank | Blank |
| 2 | 190715_Smp_FK180310S62C1-25m_A.mzML | 25m | Cyclone | Afternoon |
| 5 | 190715_Smp_FK180310S62C1-25m_A.mzML | DCM | Cyclone | Afternoon |
| 8 | 190715_Smp_FK180310S62C1-25m_A.mzML | 25m | Cyclone | Morning |
| 11 | 190715_Smp_FK180310S62C1-25m_A.mzML | DCM | Cyclone | Morning |
| 14 | 190715_Smp_FK180310S62C1-25m_A.mzML | 25m | Anticyclone | Afternoon |
| 17 | 190715_Smp_FK180310S62C1-25m_A.mzML | DCM | Anticyclone | Afternoon |
| 20 | 190715_Smp_FK180310S62C1-25m_A.mzML | 25m | Anticyclone | Morning |
| 23 | 190715_Smp_FK180310S62C1-25m_A.mzML | DCM | Anticyclone | Morning |

with triplicate samples for each one explaining the missing fileid numbers. Additional columns will appear as coloration options in the app. Three colors are established by default; you'll need to define more if you're getting crazy.

MSMS data is handled separately in the next section of the script, where you'll again need to point to your MSMS data files in .mzML format

`makeRawData` will produce 3 files from your .mzMLs. The first is simply a .csv file containing your metadata - it's usually a good idea to open this in a spreadsheet GUI and make sure everything looks okay. The next two are files produced by the calls to `saveRDS()`, which compresses the data to about a quarter of its size in memory. Reading and writing RDS is faster and safer than using the more general `save()` function.

