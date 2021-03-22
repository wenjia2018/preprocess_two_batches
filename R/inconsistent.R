library(tidyverse)

############################################################
# 2) STEVES NORMALIZATION HAS CHANGED
############################################################

############################################################
# 2019 batch 1
############################################################

# load 2019 data
batch1.1 <- "/home/share/data_input/AddHealthYr1 Plates1-12 - Gene CPM ReferenceGeneNormalized Log2.txt" 
batch1.1 <- read.table(batch1.1, stringsAsFactors = FALSE, header = TRUE)

# remove duplicate genes for simplicity of exposition
dups <- unique(batch1.1[duplicated(batch1.1[, 1]), 1])  
batch1.1 <- batch1.1[! batch1.1$Gene %in% dups, ] %>% as_tibble %>% column_to_rownames("Gene")

# convert colnames from vialID to AID
ids <- "/home/share/data_input/AID_VialID.csv"
ids <- 
  read.csv(ids, stringsAsFactors = FALSE, header = TRUE)  %>%
  mutate_all(~ str_c( "X", .))
id <- `names<-`(ids$AID, ids$VialID)
colnames(batch1.1) <- id[colnames(batch1.1)]

############################################################
# 2020 batch 1
############################################################

# load 2020 data
batch1.2 <- "/home/share/data_input/two_batches/gene.expression.genenormalized.040820/gene.expression.genenormalized.cole.batch1.040820"
batch1.2 <- read.table(batch1.2, stringsAsFactors = FALSE, header = TRUE)

# consider genes/subjects common to the two datasets
joint_cols <- intersect(colnames(batch1.1), colnames(batch1.2))
joint_rows <- intersect(rownames(batch1.1), rownames(batch1.2))

# behold: any random sample of genes/subjects have DIFFERENT values between the 2019 and 2020 datasets
i = sample(joint_rows, 8)
j = sample(joint_cols, 8)
batch1.1[i, j]
batch1.2[i, j]

# not all equal
all.equal(batch1.1[joint_rows, joint_cols], 
          batch1.2[joint_rows, joint_cols])

############################################################
# 1) THE RAW COUNTS HAVE NOT CHANGED
###########################################################

# raw 2019 data
batch1.1_raw <- "/home/share/data_input/AddHealthYr1 Plates1-12 - ReadsPerGene ENSG Raw.txt"
batch1.1_raw <- 
  read.table(batch1.1_raw, stringsAsFactors = FALSE, header = TRUE) %>% 
  as_tibble(ex) %>% 
  column_to_rownames("Gene")
colnames(batch1.1_raw) <- id[colnames(batch1.1_raw)]

# raw 2020 data
batch1.2_raw <- "/home/share/data_input/two_batches/gene.expression.batch1.batch2.notnormalizedorfiltered.030320"
batch1.2_raw <- read.table(batch1.2_raw, stringsAsFactors = FALSE, header = TRUE)

# consider genes/subjects common to the two datasets
joint_cols <- intersect(colnames(batch1.1_raw), colnames(batch1.2_raw))
joint_rows <- intersect(rownames(batch1.1_raw), rownames(batch1.2_raw))

# behold: any random sample of genes/subjects have DIFFERENT values between the 2019 and 2020 datasets
i = sample(joint_rows, 8)
j = sample(joint_cols, 8)
batch1.1_raw[i, j]
batch1.2_raw[i, j]

# all equal?
a = all.equal(batch1.1_raw[joint_rows, joint_cols], 
              batch1.2_raw[joint_rows, joint_cols])

# so... these 4 subjects have different raw counts
# eyeball their counts (col 1 is from batch 1 2019, col 2 is from batch 1 2020)
map(c("X98578343", "X57258951", "X57258951", "X91712842"), 
    function(x) {
      cbind(batch1.1_raw[joint_rows, joint_cols][x],
            batch1.2_raw[joint_rows, joint_cols][x]) %>%
        head %>% 
        print
    }) 
