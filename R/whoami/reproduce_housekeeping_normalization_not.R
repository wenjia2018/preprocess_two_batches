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
  "gplots",
  "recount", # recount data repository
  "GEOquery", # geo data repository
  "GEOmetadb", # query geo metadata
  "biobroom"
)
pkg %>%
  map(library, character.only = T)

########################################################
# REPRODUCE STEVES HOUSEKEEPING
########################################################
client = T

#source("/Volumes/Share/preprocess/R/init_data.R")
dat_ref <- readRDS("/Volumes/Share/preprocess/data/dt_steve.rds")
dat <- readRDS("/Volumes/Share/preprocess/data/dt.rds")
# cpm = function(x){
#   # x a matrix
#   10^6 * (x / colSums(x))
# }
# 
 tcounts1 = exprs(dat_ref) # steve
# tcounts2 = log2(pmax(cpm(exprs(dat)), 1)) # library size
# 
# hk = colMeans(exprs(dat[housekeepers11, ])) # mean housekeeper per subject
# tcounts3 = log2(pmax(cpm(exprs(dat) / hk), 1)) # me, trying to reproduce steve
# tcounts4 = log2(pmax(cpm(exprs(dat)) / hk , 1)) # me, trying to reproduce steve
 HG=c("C1orf43","CHMP2A","EMC7","GPI","PSMB2","PSMB4","RAB7A","REEP5","SNRPD3","VCP","VPS29")

norm.hg=function(dat){
  
  raw_counts=dat %>% exprs() 
  
   hk = colMeans(exprs(dat[HG, ])) # mean housekeeper per subject
  #hk = psych::geometric.mean (exprs(dat[HG, ])) #geometric mean
  counts_hg = (t(t(raw_counts)/hk)*10^6) %>% pmax(1) %>% log2 # me, trying to reproduce steve
  
  return(counts_hg)
}

counts_hg=norm.hg(dat)



# plot a random gene
plot(tcounts1[200,], counts_hg[200,]) 

# 
# all.equal(tcounts1, tcounts2) 
# sum(tcounts1-tcounts2 < .000001)
# raw_counts=dat %>% exprs()
# tcounts4 = (t(t(raw_counts)/hk)*10^6) %>% pmax(1) %>% log2 # me, trying to reproduce steve

N=dim(dat)[1]
out=NULL
for (i in 1:N)
out[[i]] = lm(tcounts1[i,] ~ counts_hg[i,])

a=out %>% purrr::map_df(~ .x$coefficients[2] %>% as.matrix %>% t %>% as.data.frame)



M=dim(dat)[2]
out2=NULL
for (i in 1:N)
  out2[[i]] = lm(tcounts1[,i] ~ counts_hg[,i])

b=out2 %>% purrr::map_df(~ .x$coefficients[2] %>% as.matrix %>% t %>% as.data.frame)