library(haven)
library(tidyverse)

#load the data
ancestralPCA <- read.delim("/Volumes/Data/addhealth/addhealthdata/waves_1_4/DNA_PCA_ancestory/ancestralPCA.eigenStratPCs.allparticipants.082719")
biorel <- read.delim("/Volumes/Data/addhealth/addhealthdata/waves_1_4/DNA_PCA_ancestory/relatedness.allparticipants.082719")

#change aid label for the coherence
ancestralPCA <- ancestralPCA  %>% rename("AID" ="aid")

ancestralPCA <-ancestralPCA %>% rename_at(vars( starts_with("X")), list( ~(str_replace(., "X", "AncestryPC") ) ))

#change aid label for the coherence and create a dummy for biological relatedness
biorel <- biorel %>% mutate(biorels = case_when(sample=="AddHealthUnrel" ~ 0,
                                      sample=="RelsAddHealth" ~ 1)) %>% 
                     rename("AID" ="aid")

ancestry_biorel <- ancestralPCA %>% left_join(biorel, by="AID") 
ancestry_biorel$AID<-as.character(ancestry_biorel$AID)

ancestry_biorel %>% saveRDS(file = "/Volumes/share/preprocess/data/ancestry_biorel.rds")

