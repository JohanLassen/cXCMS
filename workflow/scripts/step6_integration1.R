

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 1

input  <- args[1]
output <- args[2]


fcp           <- FillChromPeaksParam(expandMz = 0, expandRt = 0)
cXCMS::cFillChromPeaksStep1(input = input, fcp = fcp, output = output)