

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak grouping

input <- args[1]
output <- args[2]

   
msobject <- readr::read_rds(input)
pdp <- xcms::PeakDensityParam(sampleGroups = msobject@phenoData@data$group, 
			maxFeatures  = 50,
                        bw           = 5, 
                        minFraction  = 0.7,
                        binSize = 0.01)
msobject <- groupChromPeaks(msobject, pdp)


readr::write_rds(msobject, output)