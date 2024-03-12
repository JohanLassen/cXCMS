

library(tidyverse)
library(xcms)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 2

input  <- args[1]
output <- args[2]
index  <- as.numeric(args[3])

print(input)
cXCMS::cFillChromPeaksStep2(input = input, output = output, index = index)

# Debugging
# > output <- "./tmp/peak_integrated/test1"
# > output <- "./tmp/peak_integrated/test1.rds"
# > index <- 1
# > 
#   cXCMS::cFillChromPeaksStep2(input = input, output = output, index = index)
# 
