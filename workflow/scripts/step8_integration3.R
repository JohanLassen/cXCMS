

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 1

inputs <- list.files(args[1], full.names=T)
output <- args[2]
step1_input <- args[3]

cXCMS::cFillChromPeaksStep3(inputFromStep1 = step1_input, inputFromStep2 = inputs, output = output)
