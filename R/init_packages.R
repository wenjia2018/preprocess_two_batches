
########################################################
# THIS SCRIPT SIMPLY LOADS SOME USEFUL R PACKAGES
# PACKAGES SHOULD BE INSTALLED BY ALTERING init_install_packages.R
# BEFORE BUILDING DOCKER IMAGE.
########################################################

library(tidyverse)

########################################################
# BIOCONDUCTOR
########################################################

pkg <- c(
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
  "compositions",
  "naniar",
  "foreign",
  "lubridate",
  "fastDummies",
  "readxl",
  "haven"
)
pkg %>%
  map(library, character.only = T)

########################################################
# CRAN
########################################################

pkg_cran <- c(
  "ggformula",
  "modelr",
  "DescTools",
  "car"
)
pkg_cran %>%
  map(library, character.only = TRUE)

########################################################
# GITHUB
########################################################

library(dbr)

########################################################
# OTHER
########################################################

# model.matrix() FUSSY ABOUT NA
old_na_action <- options("na.action")
options(na.action = "na.pass")
# options(na.action = old_na_action)
