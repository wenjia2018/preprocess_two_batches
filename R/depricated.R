# INSTEAD OF SOURCING init_data.R, just read dat.rds from disk

# DEPRICATED!! PUT THE FOLLOWING LINES IN YOUR ANALYSIS SCRIPT

########################################################
# EXTRACT PHENOTYPE DATA FROM EXPRESSION SET
########################################################

phen <- dat %>% pData()
phen$male <- phen$sex_interv_m # BEWARE NAME CHANGE!!!!!!
phen$raceth <- phen$race_interv # BEWARE NAME CHANGE!!!!!!
