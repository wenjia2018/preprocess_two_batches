########################################################
# PACKAGES
########################################################

# Packages
library(tidyverse)
# tidyverse_update()

########################################################
# BIOCONDUCTOR
########################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite()
pkg = c(
  "DESeq2",
  "edgeR",
  "limma",
  "Biobase",
  "recount",
  "globaltest",
  "gplots",
  "recount", # recount data repository
  "GEOquery", # geo data repository
  "GEOmetadb", # query geo metadata
  "biobroom",
  "compositions"
)
pkg %>%
  map(library, character.only = T)

########################################################
# REPRODUCE STEVES HOUSEKEEPING
########################################################
client = T
source("/Volumes/Share/preprocess/init_data.R")

cpm = function(x){
  # x a matrix
  10^6 * (x / colSums(x))
}

tcounts1 = exprs(dat_ref) # steve
tcounts2 = log2(pmax(cpm(exprs(dat)), 1)) # library size
all.equal(tcounts1, tcounts2) 

hk = colMeans(exprs(dat[housekeepers11, ])) # mean housekeeper per subject
tcounts3 = log2(pmax(cpm(exprs(dat) / hk), 1)) # me, trying to reproduce steve
tcounts4 = log2(pmax(cpm(exprs(dat)) / hk , 1)) # me, trying to reproduce steve


plot(tcounts1[20,], tcounts3[20,])
plot(tcounts1[200,], tcounts3[200,])
plot(tcounts1[2000,], tcounts3[2000,])

plot(tcounts1[20,], tcounts4[20,])
