merge_pheno_geno = function(geno_fname, pheno_fname, rna_input_fname = NULL){
  
  ########################################################
  # FILTERING
  ########################################################
  
  # PHENOTYPE
  waves  <- readRDS(str_c(data_output, "/", pheno_fname))  
  
  # RNA 
  dat = readRDS(str_c(data_output, "/", geno_fname)) # reference gene normalization
  
  # replace Lauren with our variables 
  pData(dat) <- select(pData(dat), AID) %>% left_join(waves, by = "AID")
  
  # SAVE 
  saveRDS(
    dat, 
    str_c(data_output, 
          "/", 
          str_flatten(
            str_remove_all(c(geno_fname, pheno_fname), ".rds"), "_"), 
          ".rds"
    )
  )  
}