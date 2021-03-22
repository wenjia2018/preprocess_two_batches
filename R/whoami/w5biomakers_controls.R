library(haven)
library(tidyverse)
library(lubridate)

w5biocovars <- read_sas("/Volumes/Data/addhealth/addhealthdata/wave5/biomarkers/w5biocovars.sas7bdat", NULL)
waves  <- readRDS("/Volumes/share/preprocess/data/waves.rds")  # path relative to server root

w5biocovars_wave <- w5biocovars %>% left_join(waves, by="AID")

w5biocovars_wave <- w5biocovars_wave %>% mutate(pregnant_biow5     = case_when (sex_interv_m==1 ~ 0, #male by default not pregnant, don't know are assigned as missing
                                                                               sex_interv_m==0 & Q012==2  ~ 0, #female and not pregnant
                                                                               sex_interv_m==0 & Q012==1  ~ 1), #female and pregnant 
                                                illness_4wks_biow5 = case_when (Q013a==1 | Q013b==1 | Q013c==1 | Q013d==1 | Q013e==1 | Q013f==1  ~ 1, #any illness in the past 4 weeks
                                                                               Q013a>=2 | Q013b>=2 | Q013c>=2 | Q013d>=2 | Q013e>=2 | Q013f>=2  ~ 0),
                                                illness_2wks_biow5 = case_when (Q014a==1 | Q014b==1 | Q014c==1 | Q014d==1 | Q014e==1 | Q014f==1  ~ 1, #any illness in the past 4 weeks
                                                                               Q014a>=2 | Q014b>=2 | Q014c>=2 | Q014d>=2 | Q014e>=2 | Q014f>=2  ~ 0),
                                                smoking_biow5     = case_when (Q015==1  ~ 1, #smoking at the time of interview
                                                                               Q015==2  ~ 0),
                                                time_biow5        =  EXAMDATE, #if we want to include as linear predictor, otherwise we should create categories
                                                kit_biow5         = case_when (KITCOND==0 ~ 0, # kit ok
                                                                               KITCOND>0 ~ 1), #some problems with the kit
                                                tube_biow5        = case_when (TUBECOND==0 ~ 0, # tube ok
                                                                               TUBECOND>0 ~ 1), #some problems with the tube
                                                more48h_biow5     = case_when (AlqQuality=="> 48 hours" ~ 1, #something that went wrong with timing (ask Brandt, no info in the documentation)
                                                                                AlqQuality=="Normal" ~ 0),
                                                months_biow5      = as.factor (month(EXAMDATE)),
                                                hour_biow5        = as.factor (hour(EXAMTIME)),
                                                travel_biow5      = case_when (Q016==1  ~ 1,
                                                                               Q016>1   ~ 0)
                                                )
w5biocovars_wave$AID<-as.character(w5biocovars_wave$AID)

w5biocovars_constr <- w5biocovars_wave %>% select(pregnant_biow5, illness_4wks_biow5, illness_2wks_biow5, smoking_biow5, kit_biow5, 
                                                  tube_biow5, more48h_biow5, FastHrs, months_biow5, hour_biow5, travel_biow5, AID)

w5biocovars_constr %>% saveRDS(file = "/volumes/share/preprocess/data/w5biovars.rds")