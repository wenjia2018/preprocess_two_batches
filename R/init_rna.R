init_rna <- function(rna_input_fname, geno_fname) {

  ########################################################
  # Same data files, different paths/servers
  ########################################################
  
  # PHENOTYPE FROM LAUREN, COVARIATES: SEE EMAIL.
  ph <- 
    haven::read_dta(str_c(data_input, "180626_RNAControls.dta"))  %>% 
    dplyr::rename(AID = aid) 
  
  # RAW EXPRESSION
  ex <- read.table(str_c(data_input, rna_input_fname), stringsAsFactors = FALSE, header = TRUE)
   
  # Gene nomenclature: if in ENSG, then convert to HUGO (Steve's data is already in HUGO)
  if(str_detect(rownames(ex[1, ]), "ENSG")) { 
     
    if(full_reproducibility <- FALSE) {
      ## The biomart query used to generate G_list below
      library('biomaRt')
      glist  = rownames(ex) # gene names
      mart   = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
      G_list = getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values = glist, mart= mart)
      saveRDS(G_list, "/home/share/data_input/genename_15062020.rds")
    }
    
    G_list = readRDS("/home/share/data_input/genename_15062020.rds") # from Wenjia Xu
    G_list = G_list  %>% 
      filter(hgnc_symbol != "") %>% 
      filter(!duplicated(ensembl_gene_id) & !duplicated(hgnc_symbol)) %>% 
      filter(ensembl_gene_id %in% rownames(ex)) %>% 
      tidyr::unite("hgnc_ensg", c("hgnc_symbol", "ensembl_gene_id"), remove = FALSE) # create a hybrid HUGO ENSG name
    ex = ex[G_list$ensembl_gene_id, ] # ensure correct row order
    rownames(ex) = G_list$hgnc_symbol # alternatively could take the value G_list$hgnc_ensg, for a hybrid name
    
  }
  
  if(0){ 
    # This is the gene set in Brandts
    brandts_genes = readRDS("/home/share/data_input/brandts-genes.rds")
    ex = ex[rownames(ex) %in% brandts_genes, ]
  }
  
  # QUALITY
  quality = readRDS(str_c(data_input, "quality.rds")) 
  
  # SELECT THOSE AUGMENTED PHENOTYPE MATRIX WITH CORRESPONDING SAMPLES
  ph <- ph %>% 
    left_join(quality, by = "AID") %>% 
    mutate(XAID = str_c("X", AID)) %>% 
    as.data.frame
  rownames(ph) <- str_c("X", ph$AID) # needed by ExpressionSet() function
  ph <- ph[colnames(ex), ]
  
  # MAKE AN EXPRESSION SET
  all.equal(colnames(ex), rownames(ph))
  phenData <- new("AnnotatedDataFrame", data = ph)
  dat <- ExpressionSet(
    assayData = ex %>% as.matrix(),
    phenoData = phenData
  )
  
  saveRDS(dat, file = str_c(data_output, str_c("/", geno_fname)))
}
