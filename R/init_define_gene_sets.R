init_define_gene_sets <- function(expression_fname){ 
  
  method = str_extract(expression_fname, "steve|quantile|raw") 
  dat = readRDS(file.path(data_output, expression_fname ))
  gene_signature <- featureNames(dat)
  
  # specify which normalization method
    path <- str_c(gsub("(/.+?)/.*", "\\1", data_input)) 

  
  ########################################################
  # HOUSEKEEPERS
  ########################################################
  # Eisenberg, E., & Levanon, E. Y. (2013). Human housekeeping genes, revisited. Trends in Genetics, 29(10), 569-574.
  housekeepers11 <- c(
    "C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2",
    "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29"
  )
  
  ########################################################
  # DEFINE CTRA SET (FROM LITERATURE)
  ########################################################
  
  ctra_original <-
    c(
      "IL1A", "IL1B", "IL6", "IL8", "TNF", "PTGS1", "PTGS2", "FOS", "FOSB", "FOSL1",
      "FOSL2", "JUN", "JUNB", "JUND", "NFKB1", "NFKB2", "REL", "RELA", "RELB",
      "GBP1", "IFI16", "IFI27", "IFI27L1", "IFI27L2", "IFI30", "IFI35", "IFI44",
      "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFIT1L", "IFITM1",
      "IFITM2", "IFITM3", "IFITM4P", "IFITM5", "IFNB1", "IRF2", "IRF7", "IRF8",
      "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IGJ", "IGLL1", "IGLL3"
    )
  
  (ctra_available <- intersect(gene_signature, ctra_original))
  
  #######################################################
  # DEFINE CTRA GENE SET
  # (FROM MY PREVIOUS PROJECT, WHICH UPDATED SOME OUTDATED GENE NAMES)
  ########################################################
  
  # by their HGNC names (HUGO gene nomenclature committee).
  inflammatory <- c(
    "IL1A", "IL1B", "IL6", "IL8", "TNF", "PTGS1", "PTGS2", "FOS",
    "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "NFKB1",
    "NFKB2", "REL", "RELA", "RELB"
  )
  
  interferonTypeI <- c(
    "GBP1", "IFI16", "IFI27", "IFI27L1", "IFI27L2", "IFI30",
    "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1",
    "IFIT2", "IFIT3", "IFIT5", "IFIT1L", "IFITM1", "IFITM2",
    "IFITM3", "IFITM4P", "IFITM5", "IFNB1", "IRF2", "IRF7",
    "IRF8", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL"
  )
  
  # antibody = c("IGJ", "IGLL1", "IGLL3")
  antibody <- c("JCHAIN", "IGLL1", "IGLL3P")
  
  
  AntBIntF <- c(
    "GBP1", "IFI16", "IFI27", "IFI27L1", "IFI27L2", "IFI30",
    "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFIT1L", "IFITM1", "IFITM2", "IFITM3", "IFITM4P", "IFITM5",
    "IFNB1", "IRF2", "IRF7", "IRF8", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "OASL", "JCHAIN", "IGLL1", "IGLL3P"
  )
  
  
  
  
  ctra0 <- c(inflammatory, interferonTypeI, antibody)
  
  ctraCore0 <- c(
    "IRF7", "JUN", "IGJ", "IL8", "IL1B", "FOSB", "FOSL2", "IFIT3",
    "IFI35", "IFI44L", "MX1", "OAS2"
  )
  
  ctraCore <- c(
    "IRF7", "JUN", "JCHAIN", "CXCL8", "IL1B", "FOSB", "FOSL2", "IFIT3",
    "IFI35", "IFI44L", "MX1", "OAS2"
  )
  
  ########################################################
  # INTERSECTION OF ABOVE SETS WITH THE GENES IN THIS SAMPLE
  ########################################################
  ctra_available <- intersect(gene_signature, ctra0)
  inflammatory <- inflammatory %>% intersect(ctra_available)
  interferonTypeI <- interferonTypeI %>% intersect(ctra_available)
  antibody <- antibody %>% intersect(ctra_available) #
  AntBIntF <- AntBIntF %>% intersect(ctra_available)
  ########################################################
  # CELL TYPE MARKERS
  ########################################################
  # From Steve's email: the set of genes used as covariates to crudely control for leukocyte subset prevalence is
  
  cell_types <-
    list(
      T_0 = c("CD3D", "CD3E"), # for T cells (only occasionally the "CD3E")
      T_CD4 = "CD4", # for the CD4+ subset of T cells
      T_CD8 = "CD8A", # for the CD8+ subset of T cells
      B = "CD19", # for B cells
      NK = c("FCGR3A", "NCAM1"), # for NK cells
      M = "CD14"
    ) # for monocytes
  
  cell_types_one <- map_chr(cell_types, 1) # hack: for now just take one gene for each type
  
  ########################################################
  # CECILIA
  ########################################################
  
  # by their HGNC names (HUGO gene nomenclature committee).
  infarction <- c(
    "ABCA1", "ACE", "ADD1", "AGT", "AGTR1", "ALOX5AP", "APOA1", "APOA5", "APOE", "CCL11", "CCR2", "CCR5", "CD14", "CETP", "COMT", "CX3CR1",
    "CYP11B2", "CYP2C9", "ENPP1", "ESR1", "F12", "F13A1", "F2", "F5", "FGB", "FTO", "GJA4", "GP1BA", "NR3C1", "GSTT1", "HFE", "HSL", "HNF1",
    "HTR2A", "ICAM1", "IL1B", "IL6", "IL18", "ITAGA2", "ITGB3", "KIF6", "LDL R", "LIPC", "LPL", "LRP1", "MGP", "MMP3", "MTHFR", "MTP", "MTR",
    "OLR1", "p22-PHOX", "PAI1", "PECAM", "PON1+2", "PPARG", "PTGS2", "RECQL2", "SELE", "SELP", "TFPI", "THBD", "TLR4", "TNF", "TNFRSF1A", "UCP-2"
  )
  
  # depression <- c("MSRA", "FDFT1", "C8orl12", "c8orl13", "MTMR9", "BLK", "MFHAS1")
  loneliness <- c(
    "TCF4", "PHF2", "NMUR2", "EPB41L2", "LOC100499466", "OSTF1", "ARFGEF2", "GBE1", "ERBB4", "NMUR2", "OR1S1", "HIVEP1", "NDUFS3",
    "RAB9BP1", "MBD5"
  )
  social <- c(
    "BARHL2", "DPYD", "SLC4A10", "NUP35", "GNRHR", "ADH1B", "MIR2113", "DPP6", "ZNF462", "PAX2", "C17orf112", "PMAIP1", "TNRC6B", "NFIA",
    "BARHL2", "LRRN2", "KLHL29", "OTX1", "LONRF2", "CNTNAP5", "MST1", "IFT57", "PPA2", "LINC00461", "ODZ2", "FBXL4", "TAC1", "LHFPL3", "TNKS",
    "PCDH17", "TCF4", "KANK4", "THSD7B", "CAMKV", "CADM2", "MIR1275", "MIR147A"
  )
  diabetes <- c(
    "CD101", "CEP68", "EHHADH", "RP11-10L12.4", "ANKH", "POC5", "RREB1", "MICB", "HLA-DQB1", "CENPW", "ARG1", "MED23", "TP53INP1", "RPL8",
    "CAMK1D", "CWF19L1", "SNORA12", "PLEKHA1", "SSSCA1", "ARAP1", "P2RX4", "CAMKK2", "C15orf38", "RCCD1", "ANKFY1", "ATP5G1", "UBE2Z", "PABPC4", "PABPC4",
    "LTA", "CUTA", "ARG1", "TP53INP1", "NUDT5", "CAMK1D", "LOC283070", "CWF19L1", "PLEKHA1", "KLHDC5", "P2RX4", "CAMKK2", "CAMKK2", "C15orf38",
    "ATP5G1", "ATP5G1"
  )
  alzeheimer <- c(
    "CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7", "APOE", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3", "INPP5D",
    "MEF2C", "NME8", "ZCWPW1", "CELF1", "FERMT2", "CASS4"
  )
  
  # gene list saved on the server from outside sources
  
  
  genelist_diabetes <- read_excel(str_c(path, "/share/ext/data/T2Diabetes Xue 2018.xlsx"),
                                  skip = 1, col_types = c("guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "numeric", "text")
  )
  
  genelist_cvd <- read_excel(str_c(path, "/share/ext/data/CVD and MI Nickpay et al 2015.xlsx"),
                             skip = 1, col_types = c("guess", "numeric", "numeric")
  )
  
  prior_genes <- readRDS(str_c(path, "/share/ext/data/gene_list_inflam.rds")) %>%
    as_vector() %>%
    as.character()
  
  
  
  ######## 1k inflame gene
  prior_genes <- intersect(gene_signature, prior_genes)
  ########
  
  ###### diabetes clean gene sets generate
  # construct the dataset getting rid of the duplicates
  genediab1 <- genelist_diabetes %>%
    group_by(`Mapped Gene`) %>%
    filter(n() <= 1) %>%
    ungroup()
  genediab2 <- genelist_diabetes %>%
    group_by(`Mapped Gene`) %>%
    filter(n() > 1) %>%
    arrange(`Mapped Gene`, P) %>%
    ungroup() %>%
    group_by(`Mapped Gene`) %>%
    top_n(-1, P)
  
  genediab <- bind_rows(genediab1, genediab2)
  
  # give different weights to the genes in different quartiles
  genediab <- genediab %>%
    mutate(P_bin = cut(percent_rank(P), c(0, .25, .5, .75, 1.01), right = FALSE, ordered_result = TRUE)) %>%
    mutate(p_quartile = as.integer(P_bin)) %>%
    mutate(stringTF = str_detect(`Mapped Gene`, ",{1}")) %>%
    mutate(weigths = ifelse(p_quartile == 1, 0.4, ifelse(p_quartile == 2, 0.3, ifelse(p_quartile == 3, 0.2, 0.1)))) %>%
    as_tibble()
  
  # duplicate rows with more than one gene
  duplicates1a <- genediab %>% filter(stringTF == TRUE)
  duplicates1a$`Mapped Gene` <- str_extract(duplicates1a$`Mapped Gene`, ".{1,}(?=,)")
  duplicates1a <- duplicates1a %>%
    mutate(stringTF2 = str_detect(.$`Mapped Gene`, ",{1}")) %>%
    as_tibble()
  duplicates1left <- duplicates1a %>%
    filter(stringTF2 == FALSE) %>%
    as_tibble()
  
  duplicates2a <- genediab %>% filter(stringTF == TRUE)
  duplicates2a$`Mapped Gene` <- str_extract(duplicates2a$`Mapped Gene`, "(?<=,).{1,}")
  duplicates2a <- duplicates2a %>%
    mutate(stringTF2 = str_detect(.$`Mapped Gene`, ",{1}")) %>%
    as_tibble()
  duplicates2left <- duplicates2a %>% filter(stringTF2 == FALSE)
  duplicates2b <- duplicates2a %>%
    filter(stringTF2 == TRUE) %>%
    as_tibble()
  duplicates2b$`Mapped Gene` <- str_extract(duplicates2b$`Mapped Gene`, "(?<=,).{1,}")
  
  # aggregating everything
  genediab_def <- genediab %>%
    filter(stringTF == FALSE | is.na(stringTF)) %>%
    bind_rows(duplicates1left, duplicates2left) %>%
    filter(n() > 1) %>%
    arrange(`Mapped Gene`, P) %>%
    ungroup() %>%
    group_by(`Mapped Gene`) %>%
    top_n(-1, P)
  
  # check that the genes are in the list of genes we have in our expression matrix
  (diabetes_available2 <- intersect(gene_signature, genediab_def$`Mapped Gene`))
  diabetes_available2 <- diabetes_available2 %>%
    enframe(value = "genes") %>%
    arrange(genes)
  diabetes2 <- genediab_def$`Mapped Gene` %>% intersect(diabetes_available2$genes)
  
  ######
  
  ##### cvd clean gene sets generate
  genecvd1 <- genelist_cvd %>%
    group_by(`Locus Name`) %>%
    filter(n() <= 1) %>%
    ungroup()
  genecvd2 <- genelist_cvd %>%
    group_by(`Locus Name`) %>%
    filter(n() > 1) %>%
    arrange(`Locus Name`, P) %>%
    ungroup() %>%
    group_by(`Locus Name`) %>%
    top_n(-1, P)
  
  genecvd <- bind_rows(genecvd1, genecvd2)
  
  # give different weights to the genes in different quartiles
  genecvd <- genecvd %>%
    mutate(P_bin = cut(percent_rank(P), c(0, .25, .5, .75, 1.01), right = FALSE, ordered_result = TRUE)) %>%
    mutate(p_quartile = as.integer(P_bin)) %>%
    mutate(stringTF = str_detect(`Locus Name`, ",{1}")) %>%
    mutate(weigths = ifelse(p_quartile == 1, 0.4, ifelse(p_quartile == 2, 0.3, ifelse(p_quartile == 3, 0.2, 0.1)))) %>%
    as_tibble()
  
  # splitting in a very inefficient way the rows with more genes
  duplicates1a <- genecvd %>% filter(stringTF == TRUE)
  duplicates1a$`Locus Name` <- str_extract(duplicates1a$`Locus Name`, ".{1,}(?=,)")
  duplicates1a <- duplicates1a %>%
    mutate(stringTF2 = str_detect(.$`Locus Name`, ",{1}")) %>%
    as_tibble()
  duplicates1left <- duplicates1a %>%
    filter(stringTF2 == FALSE) %>%
    as_tibble()
  
  duplicates1b <- duplicates1a %>%
    filter(stringTF2 == TRUE) %>%
    as_tibble()
  duplicates1b$`Locus Name` <- str_extract(duplicates1b$`Locus Name`, ".{1,}(?=,)")
  
  duplicates2a <- genecvd %>% filter(stringTF == TRUE)
  duplicates2a$`Locus Name` <- str_extract(duplicates2a$`Locus Name`, "(?<=,).{1,}")
  duplicates2a <- duplicates2a %>%
    mutate(stringTF2 = str_detect(.$`Locus Name`, ",{1}")) %>%
    as_tibble()
  duplicates2left <- duplicates2a %>% filter(stringTF2 == FALSE)
  duplicates2b <- duplicates2a %>%
    filter(stringTF2 == TRUE) %>%
    as_tibble()
  duplicates2b$`Locus Name` <- str_extract(duplicates2b$`Locus Name`, "(?<=,).{1,}")
  
  duplicates2c <- duplicates2a %>%
    filter(stringTF2 == TRUE) %>%
    as_tibble()
  duplicates2c$`Locus Name` <- str_extract(duplicates2c$`Locus Name`, ".{1,}(?=,)")
  # aggregating everything
  genecvd_def <- genecvd %>%
    filter(stringTF == FALSE) %>%
    bind_rows(duplicates1left, duplicates1b, duplicates2left, duplicates2b, duplicates2c) %>%
    filter(n() > 1) %>%
    arrange(`Locus Name`, P) %>%
    ungroup() %>%
    group_by(`Locus Name`) %>%
    top_n(-1, P)
  
  
  
  # check that the genes are in the list of genes we have in our expression matrix
  cvd2 <- intersect(gene_signature, genecvd_def$`Locus Name`)
  
  
  
  ######
  
  
  # not all GENES IN THE SET OF INTEREST are here (name change?)
  (heart_available <- intersect(gene_signature, infarction))
  # (depression_available <- intersect(gene_signature, depression))
  (loneliness_available <- intersect(gene_signature, loneliness))
  (social_available <- intersect(gene_signature, social))
  (diabetes_available <- intersect(gene_signature, diabetes))
  (alzhei_available <- intersect(gene_signature, alzeheimer))
  
  # disease signatures from quantile
  
  NSCLC <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "NSCLC")$`Gene Symbol` %>% intersect(gene_signature)
  NSCLC_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "NSCLC_up")$`Gene Symbol` %>% intersect(gene_signature)
  NSCLC_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "NSCLC_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  breast_cancer <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Breast Cancer")$`Gene Symbol` %>% intersect(gene_signature)
  breast_cancer_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Breast Cancer_up")$`Gene Symbol` %>% intersect(gene_signature)
  breast_cancer_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Breast Cancer_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  Lupus <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Lupus")$`Gene Symbol` %>% intersect(gene_signature)
  # Lupus_withdup=read_excel(str_c(path,"/share/projects/DE/data/Chronic disease gene sets.xlsx", sheet = "Lupus_rename")$`Gene Symbol`%>% intersect(gene_signature)
  
  Prostate <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Prostate")$`Gene Symbol` %>% intersect(gene_signature)
  
  Colorectal <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Colorectal")$`Gene Symbol` %>% intersect(gene_signature)
  Colorectal_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Colorectal_up")$`Gene Symbol` %>% intersect(gene_signature)
  Colorectal_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Colorectal_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  Melanoma <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Melanoma")$`Gene Symbol` %>% intersect(gene_signature)
  
  Rheumatoid_Arthritis <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Rheumatoid Arthritis")$`Gene Symbol` %>% intersect(gene_signature)
  Rheumatoid_Arthritis_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Rheumatoid Arthritis_up")$`Gene Symbol` %>% intersect(gene_signature)
  Rheumatoid_Arthritis_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Rheumatoid Arthritis_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  RA_immune_up <- c("CSTF2", "CSF3R", "TGFBR2") %>% intersect(gene_signature)
  RA_metabolism_up <- c("CYP3A4", "HSD11B2") %>% intersect(gene_signature)
  RA_neuromuscular_up <- c("SNTA1", "TNNI2", "TNNT2") %>% intersect(gene_signature)
  RA_transcription_up <- c("ZNF74") %>% intersect(gene_signature)
  
  
  RA_immune_down <- c("CCR1", "PTGES", "LMAN1", "NSEP1", "FKPB1A", "IFI30", "HLA-DP1A", "B2M", "FYB", "HLA-DR4A") %>% intersect(gene_signature)
  RA_cancer_down <- c("SAT", "RAB7", "SSX3", "LAMR1", "M17S2", "PTPRA", "SAS", "LCP1", "LRP8", "S100A10") %>% intersect(gene_signature)
  RA_transcription_down <- c("CAP", "ARPC5", "ARHGDIB", "ARPC3", "HIF1A", "CHES1") %>% intersect(gene_signature)
  RA_growthfactors_down <- c("SNX2", "ZFP364", "LTBP1") %>% intersect(gene_signature)
  RA_metabolism_down <- c("OAZ1", "CYP24", "POR", "METTL1") %>% intersect(gene_signature)
  RA_bone_down <- c("CHI3L1", "BMP4") %>% intersect(gene_signature)
  
  
  Alzheimers <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Alzheimers")$`Gene Symbol` %>% intersect(gene_signature)
  Alzheimers_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Alzheimers_up")$`Gene Symbol` %>% intersect(gene_signature)
  Alzheimers_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Alzheimers_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  Asthma <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Asthma")$`Gene Symbol` %>% intersect(gene_signature)
  Asthma_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Asthma_up")$`Gene Symbol` %>% intersect(gene_signature)
  Asthma_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Asthma_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  Hypertension <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Hypertension")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  Aortic_Aneurysm <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm")$`Gene Symbol` %>% intersect(gene_signature)
  Aortic_Aneurysm_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm_up")$`Gene Symbol` %>% intersect(gene_signature)
  Aortic_Aneurysm_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  Aortic_Aneurysm_DE <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm_DE")$`Gene Symbol` %>% intersect(gene_signature)
  Aortic_Aneurysm_DE_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm_DE_up")$`Gene Symbol` %>% intersect(gene_signature)
  Aortic_Aneurysm_DE_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Aortic Aneurysm_DE_down")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  
  
  kidney_transplant_tolerance <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "kidney transplant tolerance")$`Gene Symbol` %>% intersect(gene_signature)
  
  COPD <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "COPD")$`Gene Symbol` %>% intersect(gene_signature)
  
  
  # aging genes from science paper
  aging_genes <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes.xlsx"), sheet = "aging"))$Gene %>% intersect(gene_signature)
  aging_genes_up <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes.xlsx"), sheet = "aging_up"))$GENE %>% intersect(gene_signature)
  aging_genes_down <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes.xlsx"), sheet = "aging_down"))$GENE %>% intersect(gene_signature)
  
  aging_genes_down_cl1 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl1"))$Gene %>% intersect(gene_signature)
  aging_genes_down_cl1a <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl1a"))$Gene %>% intersect(gene_signature)
  aging_genes_down_cl1b <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl1b"))$Gene %>% intersect(gene_signature)
  aging_genes_down_cl1c <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl1c"))$Gene %>% intersect(gene_signature)
  aging_genes_down_cl2 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl2"))$Gene %>% intersect(gene_signature)
  aging_genes_down_cl3 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "down_cl3"))$Gene %>% intersect(gene_signature)
  
  aging_genes_up_cl1 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "upcl1"))$Gene %>% intersect(gene_signature)
  aging_genes_up_cl2 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "upcl2"))$Gene %>% intersect(gene_signature)
  aging_genes_up_cl3 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "upcl3"))$Gene %>% intersect(gene_signature)
  aging_genes_up_cl4 <- (read_excel(str_c(path, "/share/projects/DE/data/aging_genes_pathways.xlsx"), sheet = "upcl4"))$Gene %>% intersect(gene_signature)
  
  cluster_union = aging_genes_down_cl1 %>%
    union(aging_genes_down_cl2) %>% 
    union(aging_genes_down_cl3) %>% 
    union(aging_genes_up_cl1) %>%
    union(aging_genes_up_cl2) %>%
    union(aging_genes_up_cl3) %>%
    union(aging_genes_up_cl4) 
    
  aging_cluster_complement = setdiff(aging_genes, cluster_union)
  # HPA genes
  hpa_genes <- c(
    "HSD3B2", "CYP11B1", "CYP17A1", "CYP11A1", "PRKAR1A", "MC2R", "LEPR", "TBX19", "POMC", "POU1F1", "NR3C2", "NR3C1", "PROP1", "FKBP5", "CRHR2", "LEP", "CRH", "AVPR1A", "CRHR1", "MC4R", "AVP", "HSD11B1",
    "H6PD", "SRD5A2", "DHRS9", "SRD5A1", "CYP3A4", "AKR1D1", "SERPINA6", "ACE"
  ) %>% intersect(gene_signature)
  
  gluc_genes <- c("NFKB1", "NFKB2", "NR3C1", "MMP9", "TIMP1") %>% intersect(gene_signature)
  
  reproductive_genes <- (read_excel(str_c(path, "/share/projects/DE/data/ageatmenarche.xlsx"), sheet = 1))$GENE %>% intersect(gene_signature)
  
  ### aging immune
  
  IMMAGE <- (read_excel(str_c(path, "/share/projects/DE/data/IMMAGE.xlsx"), sheet = "sTable14"))$Gene %>% intersect(gene_signature)
  IMMAGE_up <- (read_excel(str_c(path, "/share/projects/DE/data/IMMAGE.xlsx"), sheet = "up"))$Gene %>% intersect(gene_signature)
  IMMAGE_down <- (read_excel(str_c(path, "/share/projects/DE/data/IMMAGE.xlsx"), sheet = "down"))$Gene %>% intersect(gene_signature)
  
  # Depression 
  
  
  Depression <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "Depression") %>% 
    dplyr::rename(Symbol = `Symbol `) %>%
    dplyr::filter(!is.na(Symbol)) %>% 
    dplyr::mutate(Symbol = Symbol %>% str_trim())
  
  Depression = Depression$Symbol %>% intersect(gene_signature)
  
  # CKD：Chronic kidney disease
  
  CKD <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "CKD")
  CKD = CKD$Symbol %>% intersect(gene_signature)
  
  
  CKD_up <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "CKD_up")
  CKD_up = CKD_up$Symbol %>% intersect(gene_signature)
  
  CKD_down <- read_excel(str_c(path, "/share/projects/DE/data/Chronic disease gene sets.xlsx"), sheet = "CKD_down")
  CKD_down = CKD_down$Symbol %>% intersect(gene_signature)
  
  
  ############################################################
  # TIDY UP AND CALCULATE MEANS
  ############################################################
  
  
  ########################################################
  # DEFINE DEPENDENT VARIABLE SETS AND LOAD INTO THE GLOBAL ENV
  ########################################################
  outcome_variables <- c(
    "ctra", "inflame", "interferon", "antibody", "AntBIntF",
    "CVD", "diabetes", "inflam1k",
    "NSCLC", "NSCLC_up", "NSCLC_down", "NSCLC_ud",
    "breast_cancer", "breast_cancer_up", "breast_cancer_down", "breast_cancer_ud",
    "Lupus",
    "Prostate",
    "Colorectal", "Colorectal_up", "Colorectal_down", "Colorectal_ud",
    "Melanoma",
    "Rheumatoid_Arthritis", "Rheumatoid_Arthritis_up", "Rheumatoid_Arthritis_down", "Rheumatoid_Arthritis_ud",
    "RA_immune_up",
    # "RA_metabolism_up",# empty after removing genes with lower counts
    "RA_neuromuscular_up",
    # "RA_transcription_up",# empty after removing genes with lower counts
    "RA_immune_down", "RA_cancer_down", "RA_transcription_down",
    "RA_growthfactors_down", "RA_metabolism_down", "RA_bone_down",
    "Alzheimers", "Alzheimers_up", "Alzheimers_down", "Alzheimers_ud",
    "Asthma", "Asthma_up", "Asthma_down", "Asthma_ud",
    "Hypertension",
    "Aortic_Aneurysm", "Aortic_Aneurysm_up", "Aortic_Aneurysm_down", "Aortic_Aneurysm_ud", # classify TAA patient and normal
    "Aortic_Aneurysm_DE", "Aortic_Aneurysm_DE_up", "Aortic_Aneurysm_DE_down", "Aortic_Aneurysm_DE_ud", # DE genes distinguishing TAA vs. control
    "kidney_transplant_tolerance",
    "COPD",
    "aging", "aging_up", "aging_down", "aging_ud", "aging_cluster_complement",
    "IMMAGE", "IMMAGE_up", "IMMAGE_down", "IMMAGE_ud",
    "aging_down_cl1", "aging_down_cl1a", "aging_down_cl1b", "aging_down_cl1c",
    "aging_down_cl2", "aging_down_cl3",
    "aging_up_cl1", "aging_up_cl2", "aging_up_cl3", "aging_up_cl4",
    "Depression",
    "CKD", "CKD_up", "CKD_down"
  )
  
  outcome_set <-
    list(
      # cell_type   = cell_types_one, #only for compositional
      ctra = ctra_available,
      inflame = inflammatory,
      interferon = interferonTypeI,
      # housekeepers11,
      antibody = antibody, # only 1 element doesnot fit compositional method
      AntBIntF = AntBIntF,
      CVD = cvd2,
      # # # CECILIAS
      # heart = heart_available,
      # depression = depression_available,
      # RANK TOO LOW FOR THE FOLLOWING TWO APARENTLY
      # loneliness_available = loneliness_available,
      # social_available  = social_available,
      diabetes = diabetes2,
      inflam1k = prior_genes,
      
      NSCLC = NSCLC,
      NSCLC_up = NSCLC_up,
      NSCLC_down = NSCLC_down,
      
      breast_cancer = breast_cancer,
      breast_cancer_up = breast_cancer_up,
      breast_cancer_down = breast_cancer_down,
      
      Lupus = Lupus,
      # Lupus_withdup=Lupus_withdup,
      
      Prostate = Prostate,
      
      Colorectal = Colorectal,
      Colorectal_up = Colorectal_up,
      Colorectal_down = Colorectal_down,
      
      
      Melanoma = Melanoma,
      
      Rheumatoid_Arthritis = Rheumatoid_Arthritis,
      Rheumatoid_Arthritis_up = Rheumatoid_Arthritis_up,
      Rheumatoid_Arthritis_down = Rheumatoid_Arthritis_down,
      
      RA_immune_up = RA_immune_up,
      # RA_metabolism_up=RA_metabolism_up,
      RA_neuromuscular_up = RA_neuromuscular_up,
      # RA_transcription_up=RA_transcription_up,
      
      
      RA_immune_down = RA_immune_down,
      RA_cancer_down = RA_cancer_down,
      RA_transcription_down = RA_transcription_down,
      RA_growthfactors_down = RA_growthfactors_down,
      RA_metabolism_down = RA_metabolism_down,
      RA_bone_down = RA_bone_down,
      
      
      
      
      Alzheimers = Alzheimers,
      Alzheimers_up = Alzheimers_up,
      Alzheimers_down = Alzheimers_down,
      
      Asthma = Asthma,
      Asthma_up = Asthma_up,
      Asthma_down = Asthma_down,
      
      Hypertension = Hypertension,
      
      
      Aortic_Aneurysm = Aortic_Aneurysm,
      Aortic_Aneurysm_up = Aortic_Aneurysm_up,
      Aortic_Aneurysm_down = Aortic_Aneurysm_down,
      
      Aortic_Aneurysm_DE = Aortic_Aneurysm_DE,
      Aortic_Aneurysm_DE_up = Aortic_Aneurysm_DE_up,
      Aortic_Aneurysm_DE_down = Aortic_Aneurysm_DE_down,
      
      kidney_transplant_tolerance = kidney_transplant_tolerance,
      
      COPD = COPD,
      
      aging = aging_genes,
      aging_up = aging_genes_up,
      aging_down = aging_genes_down,
      aging_cluster_complement = aging_cluster_complement,
      aging_down_cl1 = aging_genes_down_cl1,
      aging_down_cl1a = aging_genes_down_cl1a,
      aging_down_cl1b = aging_genes_down_cl1b,
      aging_down_cl1c = aging_genes_down_cl1c,
      aging_down_cl2 = aging_genes_down_cl2,
      aging_down_cl3 = aging_genes_down_cl3,
      aging_up_cl1 = aging_genes_up_cl1,
      aging_up_cl2 = aging_genes_up_cl2,
      aging_up_cl3 = aging_genes_up_cl3,
      aging_up_cl4 = aging_genes_up_cl4,
      
      IMMAGE = IMMAGE,
      IMMAGE_up = IMMAGE_up,
      IMMAGE_down = IMMAGE_down,
      
      Depression = Depression,
      CKD = CKD,
      CKD_up = CKD_up,
      CKD_down = CKD_down
    )
  #### simple take the mean over the genes in each diesease signature for each AID
  
  ##### add mRNA to each signature to avoid misunderstanding
  outcomes_newname <- outcome_set %>%
    names() %>%
    map(~ .x %>% str_c("_mRNA"))
  outcome_set <- outcome_set %>% setNames(outcomes_newname)
  outcome_variables <- outcome_variables %>% map(~ .x %>% str_c("_mRNA"))
  
  if(0){ 
    # NOW TAKE THE MEAN OF EACH GENE SET AND STORE IT IN A DATAFRAME
    if (method == "steve") {
      if (only2019 <- FALSE) {
        outcome_signature <- outcome_set %>%
          map(
            ~ dat[.x] %>%
              exprs() %>%
              apply(2, mean, na.rm = TRUE)
          ) %>%
          as.data.frame() %>%
          rownames_to_column(var = "VialID") %>%
          left_join(dat@phenoData@data %>% dplyr::select(AID, VialID), by = "VialID") %>%
          dplyr::select(-matches("VialID"))
        num_genes <- outcome_set %>%
          map(
            ~ length(.x)
          ) %>%
          as.data.frame() %>%
          gather(outcome_set, n_genes) %>%
          as_tibble()
      } else {
        outcome_signature <- outcome_set %>%
          map(
            ~ dat[.x] %>%
              exprs() %>%
              apply(2, mean, na.rm = TRUE)
          ) %>%
          as.data.frame() %>%
          rownames_to_column(var = "AID") %>%
          mutate(AID = substring(AID, 2, 100) %>% as.character())
        
        num_genes <- outcome_set %>%
          map(
            ~ length(.x)
          ) %>%
          as.data.frame() %>%
          gather(outcome_set, n_genes) %>%
          as_tibble()
      }
    }
    
    if (method == "quantile") {
      outcome_signature <- outcome_set %>%
        map(
          ~ dat[.x] %>%
            exprs() %>%
            apply(2, mean, na.rm = TRUE)
        ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "AID") %>%
        mutate(AID = substring(AID, 2, 100) %>% as.character())
      
      num_genes <- outcome_set %>%
        map(
          ~ length(.x)
        ) %>%
        as.data.frame() %>%
        gather(outcome_set, n_genes) %>%
        as_tibble()
    }
    
    
    if (method == "raw") {
      outcome_signature <- outcome_set %>%
        map(
          ~ dat[.x] %>%
            exprs() %>%
            apply(2, sum, na.rm = TRUE)
        ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "AID") %>%
        mutate(AID = substring(AID, 2, 100) %>% as.character())
      
      num_genes <- outcome_set %>%
        map(
          ~ length(.x)
        ) %>%
        as.data.frame() %>%
        gather(outcome_set, n_genes) %>%
        as_tibble()
    }
    
    
    if(jc_suggestion <- FALSE){ 
      if (method == "steve") {
        if (only2019 <- FALSE) {
          outcome_signature <- outcome_set %>%
            map(
              ~ dat_ref[.x] %>%
                exprs() %>%
                apply(2, mean, na.rm = TRUE)
            ) %>%
            as.data.frame() %>%
            rownames_to_column(var = "VialID") %>%
            left_join(dat_ref@phenoData@data %>% dplyr::select(AID, VialID), by = "VialID") %>%
            dplyr::select(-matches("VialID"))
          num_genes <- outcome_set %>%
            map(
              ~ length(.x)
            ) %>%
            as.data.frame() %>%
            gather(outcome_set, n_genes) %>%
            as_tibble()
        } else {
          outcome_signature <- outcome_set %>%
            map(
              ~ dat_ref[.x] %>%
                exprs() %>%
                apply(2, mean, na.rm = TRUE)
            ) %>%
            as.data.frame() %>%
            rownames_to_column(var = "AID") %>%
            mutate(AID = substring(AID, 2, 100) %>% as.character())
          
          num_genes <- outcome_set %>%
            map(
              ~ length(.x)
            ) %>%
            as.data.frame() %>%
            gather(outcome_set, n_genes) %>%
            as_tibble()
        }
      }
      
      if (method == "quantile") {
        outcome_signature <- outcome_set %>%
          map(
            ~ dat[.x] %>%
              exprs() %>%
              apply(2, mean, na.rm = TRUE)
          ) %>%
          as.data.frame() %>%
          rownames_to_column(var = "AID") %>%
          mutate(AID = substring(AID, 2, 100) %>% as.character())
        
        num_genes <- outcome_set %>%
          map(
            ~ length(.x)
          ) %>%
          as.data.frame() %>%
          gather(outcome_set, n_genes) %>%
          as_tibble()
      }
      
      
    }
  }
  
  # return(output = list(outcome_signature = outcome_signature,
  #                      num_genes = num_genes, 
  #                      outcome_set = outcome_set))
  return(output = list(outcome_set = outcome_set))
}