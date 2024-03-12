
library(tidyverse)
library(xcms)
library(CAMERA)
devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# ------------> Brug output fra step 8 som input
input  <- args[1]
output <- args[2]

object <- readr::read_rds(input)
rules     <- read.csv("./scripts/CAMERA_rules/rules_pos.csv")

# Converting new XCMS object to old class
xset_integrated            <- as(object, "xcmsSet")
sampclass(xset_integrated) <- pData(object)$group
sampnames(xset_integrated) <- pData(object)$sample_name

xset_annota <- diffreport(xset_integrated, filebase="./",value="into",sortpval=FALSE) 
xsa      <- xsAnnotate(xset_integrated)
xsaF     <- groupFWHM(xsa, perfwhm=0.3,intval="into")
xsaI     <- findIsotopes(xsaF,mzabs=0.001,ppm=3,intval="into",maxcharge=8,maxiso=5,minfrac=0.2)
peaklist <- getPeaklist(xsaI)


# -----------------> følgende filer skal gemmes det korrekte sted. Måske i en folder der hedder results. Kan du huske hvordan vi lavede tmp folderen i workflow.py?
write.csv(peaklist,file="./results/peaklist.csv")
save(xsaI, file = "./results/xsaI.Rdata")

xsaIC     <- groupCorr(xsaI,cor_eic=0.75)

# -----------------> Har du uploadet en rules_pos file og ligger den under CAMERA_rules? ---> Ja

xsaFA     <- findAdducts(xsaIC,polarity="positive",rules=rules, multiplier=2,ppm=3,mzabs=0.001)
peaklist  <- getPeaklist(xsaFA)
diffrep   <- cbind(xset_annota,peaklist[,c("isotopes","adduct","pcgroup")])

readr::write_csv(diffrep,file = output)

