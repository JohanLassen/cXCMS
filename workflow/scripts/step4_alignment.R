

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak alignment

input <- args[1]
output <- args[2]

msobject <- readr::read_rds(input)
   
pdp <- xcms::PeakGroupsParam(smooth = "loess",
                         minFraction	= 0.9,
                         span = 0.6,
                         extraPeaks = 5,
                         family = "gaussian")

msobject <- adjustRtime(msobject, pdp)


readr::write_rds(msobject, output)
