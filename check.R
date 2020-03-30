
# Sanity checks which currently fail for me:

# 1) log cpm and quantile normalization really can't yield the same values.
# 2) Steve apparently uses a different normalization between his original batch 1 release and the latest batch 1 release (which is now perhaps log cpm). (This is established by the code below. The code does require access to our server to execute, but the filename should be familiar.)

# The only thing that makes sense for me right now is to work from the raw count data, having applied your filters for dodgy  gene and subjects. 
# Regarding the latter, I have provisionally taken the good genes and subject ids directly from the rows and collumns of your dataset, which I currently believe is called "gene.expression.batch1.batch2.030320" (despite the wording of Kathie's email). 

# For @Wenjia and @Cecilia, this is in /home/share/preprocessed_two_batches/dt_two.rds.
# There are 8139 features, 2492 samples.

##############
# Revisit Steve's normalization(s)
##############

# Steve's original batch1 release
steve_batch1 <- read.table("../../share/data_input/AddHealthYr1 Plates1-12 - Gene CPM ReferenceGeneNormalized Log2.txt", stringsAsFactors = FALSE, header = TRUE)
# Steve's new batch1 release 
steve_batch1_new <- read.table("../../share/data_input/gene.expression.cole.batch1.030320", stringsAsFactors = FALSE, header = TRUE)

# different dimensions, gene nomenclature and subject id (row and colnames)
dim(steve_batch1)
dim(steve_batch1_new)
steve_batch1[1:10,1:10]
steve_batch1_new[1:10,1:10]

# different values
who1 = steve_batch1 %>% colnames()
who1_new = steve_batch1_new %>% colnames()
who1 %>% glimpse()
who1_new %>% glimpse()

# BECAUSE THE ORIGINAL BATCH 1 RELEASE USED VIAL ID NOT AID FOR COLNAMES, I USE dt_steve.rds, AN EXPRESSIONSET OBJECT I MADE THAT MAKES THE AID-VIALID CORRESPONDENCE EASY AND SECURE.
# Different values for the normalized data, between Steve's original normalization and 
# For example, below we one gene SCYL3 aka ENSG00000000457 for the same 5 subjects across the two versions of batch 1: the values of the normalized data are different.
dt_steve <- readRDS("/home/share/preprocessed/dt_steve.rds")
vial_id =  colnames(exprs(dt_steve["SCYL3",1:5]))
aid     = str_c("X", dt_steve$AID[dt_steve$VialID %in% sampleNames(dt_steve[,1:5])]) # corresponding aids
rbind( aid, exprs(dt_steve["SCYL3",1:5]))
steve_batch1_new[rownames(steve_batch1_new) == "ENSG00000000457", 1:5]

###############
# OTHER DETAILS
###############

# The current batch1 release contains 8139 Hugo gene IDs. 
# The intersection between the 8139 genes of the current batch1 release and the original batch1 release is only 7905 genes
length(intersect(featureNames(dt_steve), rownames(dt_steve_batch1)))

