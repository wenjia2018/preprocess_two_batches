init_data = function(){
  
  ########################################################
  # FILTERING
  ########################################################
  
  # PHENOTYPE
  waves  <- readRDS(str_c(data_output, "data/waves.rds"))  
  
  # TO AVOID LATER MERGE CONFLICT
  waves$Plate = NULL
  waves$AvgCorrelogram100 = NULL
  
  # RNA
  dat       = readRDS(str_c(data_output, "data/dt.rds")) # expression set
  dat_ref   = readRDS(str_c(data_output, "data/dt_steve.rds")) # reference gene normalization
  
  # add other variables
  pData(dat)     <- pData(dat) %>% left_join(waves %>% dplyr::select(-matches("bmi")), by = "AID")
  pData(dat_ref) <- pData(dat_ref) %>% left_join(waves %>% dplyr::select(-matches("bmi")), by = "AID")
  
  # LOOSE ALL DARK MATTER
  dark    = str_which(featureNames(dat), "^ENSG.*")
  dat     = dat[-dark, ] 
  dat_ref = dat_ref[-dark, ] 
  
  # LOOSE SUMMARIES AND HAEMOGLOBIN
  baddies = c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped",
              "ENSG00000244734", "ENSG00000188536", "ENSG00000206172")
  dat     = dat[setdiff(featureNames(dat), baddies), ]
  dat_ref = dat_ref[setdiff(featureNames(dat), baddies),]
  
  # FIND GENES WITHOUT ENOUGH VARIATION
  e_genes = edgeR::filterByExpr(Biobase::exprs(dat))
  e_genes = names(e_genes[e_genes == TRUE])
  
  # REMOVE THOSE GENES
  dat     = dat[e_genes]
  dat_ref = dat_ref[e_genes]
  
  # SAVE
  saveRDS(dat, str_c(data_output, "data/dat.rds"))
  saveRDS(dat_ref, str_c(data_output, "data/dat_ref.rds"))
  
}