

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")


args <- commandArgs(trailingOnly=TRUE)


# Set the peak picking parameters

inputs <- list.files(args[1], full.names=T)
output <- args[2]

   
cXCMS::cFindChromPeaksStep2(input = inputs, output = output)
