install.packages("BiocManager")
BiocManager::install(version = "3.10")

########################################################
# BIOCONDUCTOR
########################################################
pkg_bioc <- c(
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

purrr::map(pkg_bioc, BiocManager::install, ask = FALSE)

########################################################
# CRAN
########################################################

pkg_cran <- c(
  "ggformula",
  "modelr",
  "DescTools",
  "car",
  "naniar",
  "fastDummies",
  "foreign",
  "lubridate",
  "readxl",
  "haven"
)

purrr::map(pkg_cran, install.packages)


########################################################
# GITHUB
########################################################

install.packages("devtools")
devtools::install_github("chumbleycode/dbr")
