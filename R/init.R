# THESE FUNCTIONS TAKE NO EXPLICIT ARGUMENTS: THEY READ DATA FILES FROM DISK,
# AND WRITE CLEANED VERSIONS.

# 1. tidy expressionSet object of rna counts: dt.rds
# 2. tidy expressionSet object of ref gene normalized counts:dt_steve.rds
# 3. tidy phenotype variable (SES, BMI, etc): waves.rds
# 4. filted versions of the above expressionSet objects (no haemoglobin, etc): dat.rds, data_ref.rds


# IN ORDER FOR THE BELOW TO EXECUTE, YOUR DOCKER DAEMON MUST HAVE SUFFICIENT
# MEMORY. SEE PREFERENCE
#
# # PACKAGES
library(tidyverse)
# library(naniar)
# library(foreign)
# library(lubridate)
# library(fastDummies)
# library(readxl)
# library(tidyverse)
# library(Biobase)

# SOURCE THE CODE AND LOCATE THE INPUT DATA
if (inside_docker <- FALSE) {
  # root_work <- "/home/rstudio/workspace/"
  data_input <- "/home/rstudio/data_input/"
  data_output <- "/home/rstudio/data_output/"
  root_code = "/home/rstudio/legacy/"
} else if(on_jacobs_mac <- FALSE) {
  # just a convenience if we are working outside the docker container
  # root_work <- str_c(getwd(), "/")
  data_input <- "/Volumes/Share/data_input/"
  data_output <- "/Volumes/Share/preprocessed_two_batches"
  root_code <- str_c(getwd(), "/")
} else if(home_office <-TRUE) { 
  # just a convenience if we are working outside the docker container
  # root_work <- str_c(getwd(), "/")
  data_input <- "/home/share/data_input/"
  data_output <- "/home/share/preprocessed_two_batches"
  root_code <- str_c(getwd(), "/")
}

source(str_c(root_code, "R/init_packages.R")) # having done: source("R/init_install_packages.R")
source(str_c(root_code, "R/init_rna.R"))
source(str_c(root_code, "R/init_pheno.R"))
source(str_c(root_code, "R/init_data.R"))

# POSSIBLY COPY ALL INPUT DATA FILES INTO A SINGLE FOLDER
if (0) system("./move_data")

# *****IMPORTANT: EXECUTION ORDER OF THE FOLLOWING FUNCTIONS MATTERS: init_rna() must precede init_pheno() *******
# write data into nice expressionSet object >> dt_two.rds
init_rna(rna_input_fname = "gene.expression.batch1.batch2.notnormalizedorfiltered.030320",
         output_fname = "dt_batches1_2_raw.rds") # raw
init_rna(rna_input_fname = "gene.expression.batch1.batch2.030320", 
         output_fname = "dt_batches1_2_quantile.rds") # brandts quantile normalized
init_rna(rna_input_fname = "gene.expression.cole.batch1.030320",
         output_fname = "dt_steve_batch1.rds")  # steve ref gene normalized batch1
init_rna(rna_input_fname = "gene.expression.cole.batch2.030320",
         output_fname = "dt_steve_batch2.rds")   # steve ref gene normalized batch1  

dt_steve_batch1 <- readRDS("/home/share/preprocessed_two_batches/dt_steve_batch1.rds")
dt_steve_batch2 <- readRDS("/home/share/preprocessed_two_batches/dt_steve_batch2.rds")
dt_steve_batch1_2 <- combine(dt_steve_batch1, dt_steve_batch2)
saveRDS(object = dt_steve_batch1_2, file = str_c(data_output, str_c("/", "dt_steve_batch1_2.rds")))

init_pheno() # writes new phenotype variables >> waves.rds
init_data() # writes filtered expression set objects >> dat.rds, dat_ref.rds



