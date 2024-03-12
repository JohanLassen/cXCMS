
library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak picking parameters
cwp <- CentWaveParam(peakwidth=c(4,20), 
                    snthresh=6, 
                    ppm=12, 
                    prefilter=c(3,1000), 
                    mzdiff=0.01)

input  <- args[1]
output <- args[2]
groups <- args[3]

group_df <- readr::read_delim(file = groups, delim = ";")
group <- group_df$sample_group[sapply(group_df$sample_name, function(x){grepl(x, input)})]
print(group)

cXCMS::cFindChromPeaksStep1(input = input, output = output, cwp = cwp, groups = group)