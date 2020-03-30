init_rna <- function(rna_input_fname = "gene.expression.batch1.batch2.notnormalizedorfiltered.030320", output_fname = "dt_two.rds") {
  
  ########################################################
  # Same data files, different paths/servers
  ########################################################
  
  # PHENOTYPE FROM LAUREN, COVARIATES: SEE EMAIL.
  ph <- 
    haven::read_dta(str_c(data_input, "180626_RNAControls.dta"))  %>% 
    dplyr::rename(AID = aid) 
  
  # RAW EXPRESSION
  ex <- read.table(str_c(data_input, rna_input_fname), stringsAsFactors = FALSE, header = TRUE)
  
  # get brandt's row and colnames: these have carefully excluded dodgy genes and subjects.
  brandt = read.table(str_c(data_input, "gene.expression.batch1.batch2.030320"), stringsAsFactors = FALSE, header = TRUE)

  good_subjects = colnames(brandt)
  good_genes    = rownames(brandt)  
  
  ex = ex[rownames(ex) %in% good_genes, colnames(ex) %in% good_subjects]
  
  # Gene nomenclature
  if(full_reproducibility <- FALSE) {  
    ## The biomart query used to generate G_list below
    library('biomaRt')
    glist  = rownames(ex) # gene names
    mart   = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
    G_list = getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values = glist, mart= mart)
  }
  G_list = readRDS("/home/share/scratch/batch2/genename.rds") # from Wenjia Xu
  G_list = G_list  %>% 
    filter(hgnc_symbol != "")%>% 
    filter(gene %in% rownames(ex)) %>% 
    tidyr::unite("hgnc_ensg", c("hgnc_symbol", "gene"), remove = FALSE) # create a hybrid HUGO ENSG name
  ex = ex[G_list$gene, ] # ensure correct row order
  rownames(ex) = G_list$hgnc_symbol # alternatively could take the value G_list$hgnc_ensg, for a hybrid name
  
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
  
  saveRDS(dat, file = str_c(data_output, str_c("/", output_fname)))
}
