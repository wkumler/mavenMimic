# mavenMimic Introduction

I really like how interactive and exploratory [MAVEN](http://genomics-pubs.princeton.edu/mzroll/index.php) is, with its (relatively) intuitive interaction between chromatograms (single mass, many scans) and spectra (single scan, many masses). Combined with a TIC, it's easy to poke quickly at the data and get an initial sense of interesting compounds.

However, MAVEN also had a few drawbacks. For one, it's no longer being developed and is technically defunct. Still works great, but has a few shortcomings. Specifically, MAVEN only accepts the old .mzXML file format, doesn't handle MS<sup>2+</sup> data, and is hard (impossible?) to color-code by treatment condition.

This app mimics MAVEN's functionality while resolving the above issues. It's built on the `shiny` and `plotly` libraries in R and critically exploits the `event_data()` function to pass information back and forth between the chromatograms and the spectra. The app leverages my own code snippets for getting MS data out of weird file types (thanks to the `mzR` package) and into a simple data frame in R, which can then be `filter`ed, `mutate`d, and generally easily manipulated.

Some considerations! This app is nowhere near as fast as MAVEN - it often takes a half-second or so to render the plots after a new interaction takes place. It's also quite memory-intensive, as the 25 HILIC positive files I've been looking at from the R/V Falkor cruise occupy about a gigabyte of RAM as a `data frame` object in the R environment. One way to get around this is to use SQL and the `RSQLite` package to create a database and handle the filtering. This code exists as "oldMavenDBcode.R" but the I/O times made it almost unusable.
