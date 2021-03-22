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
source(str_c(root_code, "R/merge_pheno_geno.R"))
source(str_c(root_code, "R/init_define_gene_sets_new.R"))
# POSSIBLY COPY ALL INPUT DATA FILES INTO A SINGLE FOLDER
if (0) system("./move_data")

# *****IMPORTANT: EXECUTION ORDER OF THE FOLLOWING FUNCTIONS MATTERS: init_rna() must precede init_pheno() *******
# write data into nice expressionSet object >> dt_two.rds

latest_phenotype_file = "waves_10.09.2020.rds" # fname to write to 

# Phenotype preprocessing
init_pheno(output_fname = latest_phenotype_file) # writes new phenotype variables >> waves_date.rds

# Create expressionSet objects from RNA data recieved from Brandt (col 1), our name for this data (col 2), our phenotype data to merge in (col 3)
different_normalizations = 
  tribble(
    ~rna_input_fname,                                                ~geno_fname,                
    "gene.expression.batch1.batch2.notnormalizedorfiltered.030320", "dt_batches1_2_raw.rds",       
    "gene.expression.batch1.batch2.030320",                         "dt_batches1_2_quantile.rds",
    "gene.expression.genenormalized.batch1.batch2.040820",          "dt_batches1_2_steve.rds"
  )  

# make expressionSet objects
pwalk(different_normalizations, init_rna) 

############################################################
# merge phenotype data into expressionSet objects
############################################################

pwalk(different_normalizations, 
      merge_pheno_geno,
      pheno_fname = latest_phenotype_file) # writes filtered expression set objects  

############################################################
# decide which genes to include: gold standard uses raw counts ...
############################################################
data_2020 = 
  withr::with_dir(
    data_output,
    dir(pattern = str_c("dt.*", latest_phenotype_file)) %>% 
      set_names() %>%
      map(readRDS) 
  )

if(reconciled <- FALSE){
  good_genes = 
    exprs(readRDS(file.path(data_output, "dt_batches1_2_raw.rds"))) %>% 
    edgeR::filterByExpr() %>% 
    keep( ~ .x == TRUE) %>% 
    names()
  
  # find the common intersection of all expressionSet objects, and the good genes

  intersection = 
    data_2020 %>% 
    map(featureNames) %>%
    purrr::reduce(intersect) %>% 
    intersect(good_genes)   
  # overwrite
  data_2020 %>% iwalk(~ saveRDS(.x[intersection], file.path(data_output, .y)))
}else{
  # for unreconciled data using Steve's 1st batch feature name
  dat_ref <- readRDS("/home/share/preprocessed/dat_ref.rds")
  good_genes = dat_ref %>% featureNames()
  intersection = 
    data_2020 %>% 
    map(featureNames) %>%
    map(~ intersect(.x, good_genes))
  
  data_2020 = map2(data_2020, intersection, ~ .x[.y])
  
  data_2020 %>% iwalk(~ saveRDS(.x, file.path(data_output, .y)))
  
}

############################################################
# CREATE SIGNATURES (WENJIA)
############################################################

# CONSTRUCT SIGNATURES
signatures <- 
  names(data_2020) %>% 
  set_names() %>% 
  map(init_define_gene_sets_new)

# SAVE
iwalk(signatures, ~saveRDS(.x, file.path(data_output, str_replace(.y, ".rds", "_signature.rds"))))

############################################################
# INSPECTIONS
############################################################
# CHECK SIG SETS ARE SAME IN EACH DATASET
sigs = withr::with_dir(data_output, dir(pattern = str_c("dt.*", str_remove(latest_phenotype_file, ".rds"), "_signature"), full.names = TRUE) %>% map(readRDS))   
sigs %>% map(c("outcome_set"))  %>% map("IMMAGE_mRNA")


