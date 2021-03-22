init_pheno <- function(output_fname) {
  
  ########################################################
  # AUTHORS, JUSTIN CHUMBLEY, CECILIA POTENTE
  ########################################################

  ########################################################
  # LOAD DATA
  ########################################################
  # MULTIPLE AUTHORS OF THIS CODE: NEEDS MORE STYLE.

  # load 5 waves
  # load 5 waves
  wave1 <- read.xport(str_c(data_input, "allwave1.xpt")) %>% as_tibble
  wave2 <- read.xport(str_c(data_input, "wave2.xpt")) %>% as_tibble
  wave3 <- read.xport(str_c(data_input, "wave3.xpt")) %>% as_tibble
  wave4 <- read.xport(str_c(data_input, "wave4.xpt")) %>% as_tibble
  wave5 <- read.xport(str_c(data_input, "wave5.xpt")) %>% as_tibble
  
  # other
  w1_context     = read.xport(str_c(data_input, "Context1.xpt")) 
  w4_constructed = read.xport(str_c(data_input, "w4vars.xpt"))
  wave1_friends  = read.xport(str_c(data_input, "hfriend1.xpt")) %>% as_tibble
  wave2_friends  = read.xport(str_c(data_input, "hfriend2.xpt")) %>% as_tibble
  glucose        = read.xport(str_c(data_input, "glu_a1c.xpt")) %>% as_tibble
  sibling_w3 = read.xport(str_c(data_input, "sibling3.xpt")) %>% as_tibble
  PGS       = read.xport(str_c(data_input, "PGS_AH1.xpt")) %>% as_tibble %>% mutate_at(.vars=vars("AID"),.funs=(as.character))
  obesityclass   = readRDS(str_c(data_input, "obesityclass.rds")) 
  shawn <- haven::read_dta(str_c(data_input, "soc2000-mean-median-educ-wages.dta"))
  SEI_occ10<- read_excel(str_c(data_input,"PRESTG10SEI10_supplement-1.xls")) %>% as.tibble
  famstr_w1=read.xport(str_c(data_input, "famst.xpt")) %>% as.tibble() %>% mutate_at(.vars=vars(matches("FAMST")), .funs = list(~ .x %>% as.factor))
  w4meds=read.xport(str_c(data_input, "w4meds.xpt")) %>% as_tibble
  cvd_med <- read_excel(str_c(data_input, "cvd_med.xlsx"))
  w5biocovars <- read_sas(str_c(data_input, "w5biocovars.sas7bdat", NULL))
  ancestralPCA <- read.delim(str_c(data_input, "ancestralPCA.eigenStratPCs.allparticipants.082719"))
  biorel <- read.delim(str_c(data_input, "relatedness.allparticipants.082719"))
  ########################################################region information at wave3 and wave4 wave5
  # WX: IDEALLY HAVE ALL LOADS IN ONE PLACE
  wave1_region = read.xport(str_c(data_input, "w1homewc.xpt")) %>% as_tibble %>% transmute(AID=AID,W1REGION=REGION)
  wave2_region = read.xport(str_c(data_input, "w2homewc.xpt")) %>% as_tibble %>% transmute(AID=AID,W2REGION=REGION)
  wave3_region = read.xport(str_c(data_input, "w3region.xpt")) %>% as_tibble
  wave4_region = read.xport(str_c(data_input, "w4region.xpt")) %>% as_tibble
  
  
  # wave 5 region is included in the wave5 file, 
  # wave5_region = read.xport(str_c(data_input, "wgtsw5s1.xpt")) %>% as_tibble 
  # wave5_region = wave5_region %>% transmute(AID=AID, W5REGION=REGION %>% factor %>% fct_recode("NE" = "1", # northeast
  #                                                                                          "MW" = "2",   #midwest
  #                                                                                              "S" = "3",    #south
  #                                                                                              "W"  = "4") ) #west
  
  # blood 
  w5AID <- wave5$AID
  
  ########### quality measures (get Plate and AvgCorrelogram100)
  quality = readRDS(str_c(data_input,"quality.rds"))
  
  ###########
  anthropometric <- read.xport(str_c(data_input, "/new_phenotype/banthro5.xpt")) %>%
    as_tibble %>%
    mutate_at(.vars = c("H5BMI"), .funs = ~ifelse(H5BMI>999,NA,H5BMI))
  ########################################################region information at all waves
  
  
  waves_region=wave1_region %>% left_join(wave2_region, by = "AID") %>% left_join(wave3_region, by = "AID") %>% left_join(wave4_region, by = "AID")%>% 
    # left_join(wave5_region, by = "AID") %>% 
    mutate_at(.vars = c("W1REGION","W2REGION"), .funs =(. %>% factor%>%fct_recode("W" = "1", # west
                                                                                  "MW" = "2", #midwest
                                                                                  "S" = "3", #south
                                                                                  "NE"  = "4") )) %>% #northeast
    mutate_at(.vars = c("W3REGION","W4REGION"), .funs =(. %>% factor%>%fct_recode("NE" = "1", # northeast
                                                                                  "MW" = "2", #midwest
                                                                                  "S" = "3", #south
                                                                                  "W"  = "4") ))    #west
  ################
  
  waves = 
    wave1 %>%
    left_join(wave2, by = "AID") %>%
    left_join(wave3, by = "AID") %>%
    left_join(wave4, by = "AID") %>%
    left_join(wave5, by = "AID") %>%
    left_join(w1_context, by = "AID") %>% 
    left_join(w4_constructed, by = "AID") %>% 
    left_join(glucose, by="AID") %>% 
    left_join(quality, by="AID") %>% 
    left_join(waves_region) %>% 
    left_join(anthropometric, by ="AID") %>% 
    left_join(sibling_w3, by ="AID")
  
  waves_full = waves # full copy, waves will be trimmed below
  
  # call reformat_rna_counts.R here
  # i.e. add laurens file, make an expression set object dat
  # finally overwrite in memo.Rmd : a <- dat %>% pData
  
  rm(wave1, wave2, wave3, wave4, wave5, w4_constructed, w1_context) # duplicate
  
  
  
  
  ########################################################
  # CONSTRUCT CLEAN METRICS
  # ***NOTE*** THE FINAL DATAFRAME OF INTEREST "waves" IS DERIVED 
  # IN 4 CODE CHUNKS: THEY HAVE 
  ########################################################
  ########################################################
  # CHUNK 1: MICHAELS INITIAL FAVORITE VARIABLES
  ########################################################
  
  waves = 
    waves_full %>%
    transmute(wave5 = AID %in% w5AID, 
              # ID
              AID = AID,
              # SSS
              sss_4 = H4EC19 %>% factor(levels = 1:10),
              sss_5 = H5EC9 %>% factor(levels = 1:10), 
              sss_change=case_when(H5EC9 - H4EC19==0 ~ "nochange",
                                   H5EC9 - H4EC19> 0 ~ "poschange",
                                   H5EC9 - H4EC19< 0 ~ "negchange"),
              sss_mobility=case_when(sss_4%in%c("1","2","3") & sss_5%in%c("1","2","3") ~ "stable_low",
                                     sss_4%in%c("8","9","10") & sss_5%in%c("8","9","10") ~ "stable_high",
                                     !is.na(sss_4) & !is.na(sss_5) &  H5EC9-H4EC19 > 2 ~ "increase",
                                     !is.na(sss_4) & !is.na(sss_5) &  H4EC19-H5EC9 > 2 ~ "descrease",
                                     sss_4%in%c("4","5","6","7") & sss_5%in%c("4","5","6","7") ~ "stable_middle"
                                     
                                     
              ),
              
              # Income personal and household
              # income_4 = H4EC2, 
              income_4_house = H4EC1,
              income_5 = H5EC1,
              income_5_house = H5EC2,
              # resident mother education wave 1 and wave 2, etc
              edu_rm1 = H1RM1, 
              edu_rf1 = H1RF1,
              edu_rm2 = H2RM1, 
              edu_rf2 = H2RF1,
              poor_1 = CST90598,
              # Wave 4-5 Education
              edu_4 = H4ED2, 
              edu_5 = H5OD11, 
              # Race, interviewer rated
              race_interv = H1GI9, 
              # Sex, interviewer rated
              sex_interv = BIO_SEX,
              # wave 1 parental education
              parent_self_report = PA12,
              parent_partner_report = PB8,
              parent_income = PA55,
              # wave 1 IQ
              verbal_iq = AH_PVT,
              # wave 5 debt
              debt = H5EC6,
              # wave 4 job variables
              job_autonomy = H4LM23, 
              job_nonrepetitive = H4LM24, 
              job_supervisory = H4LM25, 
              job_satisfaction = H4LM26,
              # wave 5
              job_satisfaction_w5 = H5LM26,
              
              # wave 4 variables constructed by add health 
              # wealth
              earnings_4 = C4VAR041,# "mid point personal earnings": from: H4EC2, H4LM6, H4LM11
              income_4 = C4VAR040, 
              assets_4 = C4VAR042, 
              house_4 = C4VAR043,
              # wave 5 wealth
              new_assets_5 = H5EC4, # large gifts and inheritances
              # psychological
              depression_ah = C4VAR044, # depression 
              stress_ah = C4VAR001, # Cohen Perceived Stress Scale
              depression2_ah = C4VAR002, # CESD Depression Scale
              # cardiovascular health
              hypertension1 = C4VAR045, #hypertension stage 1
              hypertension2 = C4VAR046, #hypertension stage 2
              cholesterol = C4VAR047 #High blood cholesterol or triglycerides or lipids
    )   %>%  
    replace_with_na_all(~.x %in% c(94:99,9996, 999996, 9999997, 999998, 9999998)) %>%   # missing
    mutate(parent_income_categories = parent_income %>% ntile(5)) %>% # make quantile groups for comparability with childrens income
    left_join(dplyr::select(waves_full,AID, H4LM18), by = "AID") %>%  # wave 4 occupatianal - income data: shawn
    mutate(SOC2000 = H4LM18) %>%
    mutate( edu_rm = ifelse(!is.na(edu_rm1),edu_rm1,edu_rm2) %>% factor %>% 
              fct_collapse(
                high =  c("10",
                          "1", 
                          "2",
                          "3",
                          "4",
                          "5"),
                votec  = c("6", 
                           "7"),
                college = c("8"),
                post = c("9"),
                NULL=c("11","12")), 
            edu_rf = ifelse(!is.na(edu_rf1),edu_rf1,edu_rf2) %>%  factor %>% 
              fct_collapse(
                high =  c("10",
                          "1", 
                          "2",
                          "3",
                          "4",
                          "5"),
                votec  = c("6", 
                           "7"),
                college = c("8"),
                post = c("9"),
                NULL=c("11","12")),
            race_interv =
              race_interv %>%
              factor %>% 
              fct_recode("w" = "1", # white
                         "b" = "2", #black
                         "n" = "3", # native
                         "a"  = "4",  #asian
                         "o"  = "5",  #other
                         NULL = "6",
                         NULL = "8", 
                         NULL = "9"), 
            parent_self_report =
              parent_self_report %>% 
              factor %>%
              fct_recode("8th grade or less" ="1" ,
                         "more than 8th grade, but did not graduate from high school" = "2",
                         "went to a business, trade, or vocational school instead of high school" = "3",
                         "high school graduate" = "4",
                         "completed a GED" = "5",
                         "went to a business, trade or vocational school after high school" = "6",
                         "went to college, but did not graduate" = "7",
                         "graduated from a college or university" = "8",
                         "professional training beyond a 4-year college or university" = "9",
                         "refused" = "96", 
                         "never went to school" = "10"),
            parent_self_report_col = 
              parent_self_report %>% 
              fct_collapse(
                high =  c("never went to school",
                          "8th grade or less", 
                          "more than 8th grade, but did not graduate from high school",
                          "went to a business, trade, or vocational school instead of high school",
                          "high school graduate",
                          "completed a GED"),
                votec  = c("went to a business, trade or vocational school after high school", 
                           "went to college, but did not graduate"),
                college = c("graduated from a college or university"),
                post = c("professional training beyond a 4-year college or university")), 
            parent_partner_report = 
              parent_partner_report %>% 
              factor %>%
              fct_recode("8th grade or less" ="1" ,
                         "more than 8th grade, but did not graduate from high school" = "2",
                         "went to a business, trade, or vocational school instead of high school" = "3",
                         "high school graduate" = "4",
                         "completed a GED" = "5",
                         "went to a business, trade or vocational school after high school" = "6",
                         "went to college, but did not graduate" = "7",
                         "graduated from a college or university" = "8",
                         "professional training beyond a 4-year college or university" = "9",
                         "refused" = "96", 
                         "never went to school" = "10",
                         "school don't know" = "11",
                         "don't know" = "12",
                         "skip" = "97"),
            parent_partner_report_col = 
              parent_partner_report %>% 
              fct_collapse(
                high =  c("never went to school",
                          "8th grade or less", 
                          "more than 8th grade, but did not graduate from high school",
                          "went to a business, trade, or vocational school instead of high school",
                          "high school graduate",
                          "completed a GED",
                          "never went to school",
                          "school don't know",
                          "don't know"),
                votec  = c("went to a business, trade or vocational school after high school", 
                           "went to college, but did not graduate"),
                college = c("graduated from a college or university"),
                post = c("professional training beyond a 4-year college or university")), 
            sex_interv  =
              sex_interv  %>%
              factor %>% 
              fct_recode("m" = "1",
                         "f" = "2",
                         NULL = "6",
                         NULL = "8"),
            edu_4 = 
              edu_4 %>% 
              factor %>% 
              fct_recode("8th grade or less"  = "1"	,
                         "some high school" = "2",
                         "high school graduate"	= "3",
                         "some vocational/technical training" = "4",
                         "completed vocational/technical training" = "5",
                         "some college"	= "6",
                         "completed college (bachelor's degree)"	= "7",
                         "some graduate school" =	"8",
                         "completed a master's degree" = "9",
                         "some graduate training beyond a master's degree" = "10",
                         "completed a doctoral degree"	= "11", 
                         "some post baccalaureate professional" = "12"	,
                         "completed post baccalaureate professional"	= "13"),
            edu_5 = 
              edu_5 %>% 
              factor %>% 
              fct_recode("8th grade or less"  = "1"	,
                         "some high school" = "2",
                         "high school graduate"	= "3",
                         "GED" = "4",
                         "some vocational/technical training"  = "5",
                         "some community college" = "6",
                         "completed vocational/technical training"	= "7",
                         "associate or junior college degree" =	"8",
                         "some college" = "9",
                         "completed college (bachelor's degree)" = "10",
                         "some graduate school"	= "11", 
                         "completed a master's degree" = "12"	,
                         "some graduate training beyond a master's degree"	= "13",
                         "completed a doctoral degree" = "14",
                         "some post baccalaureate professional" = "15",
                         "completed a post baccalaureate professional" = "16"),
            edu_rm1 = 
              edu_rm1 %>% 
              factor %>% 
              fct_recode("eighth grade or less." = "1",  
                         ">eighth grade,  high school unfinished."  =  "2" ,  
                         " vocational school not high school." ="3",    
                         "high school graduate."  ="4"  ,  
                         "completed a GED." ="5"   , 
                         " vocational school after high school."   =  "6"   , 
                         " college, not graduated."=  "7"    ,
                         " graduated from a college "   =   "8"    ,
                         "professional training beyond college"    =    "9"    ,
                         "He never went to school."    = "10"   ,
                         "school, but unknown level."  =  "11"   ,
                         "schooling unknown."   ="12" ),
            edu_rf1 = 
              edu_rf1 %>% 
              factor %>% 
              fct_recode("eighth grade or less." = "1",  
                         ">eighth grade,  high school unfinished."  =  "2" ,  
                         " vocational school not high school." ="3",    
                         "high school graduate."  ="4"  ,  
                         "completed a GED." ="5"   , 
                         " vocational school after high school."   =  "6"   , 
                         " college, not graduated."=  "7"    ,
                         " graduated from a college "   =   "8"    ,
                         "professional training beyond college"    =    "9"    ,
                         "He never went to school."    = "10"   ,
                         "school, but unknown level."  =  "11"   ,
                         "schooling unknown."   ="12" ),
            income_5 = 
              income_5 %>% 
              factor %>%
              fct_recode("2500" = "1", 
                         "7500" = "2",
                         "12500" = "3", 
                         "17500" = "4",
                         "22500" = "5", 
                         "27500" = "6" ,
                         "35000" = "7",
                         "45000" = "8",
                         "62500" = "9",
                         "87500" = "10" ,
                         "125000" = "11",
                         "150000" = "12", 
                         "150000" = "13" ),
            black_white =
              race_interv %>%
              fct_collapse(o = c("n", "a", "o"))) %>% # coarser race categories 
    mutate_at(c("edu_4", "edu_5"), 
              function(x) x = factor(x, levels = levels(fct_c(.$edu_4, .$edu_5))))  %>%  
    # or use fct_c: union of levels across 2 factors
    mutate_at(.vars = c("edu_4", "edu_5"), # collapse to coarser educational categorization
              .funs = funs(col = fct_collapse(.,
                                              high = c("8th grade or less", 
                                                       "some high school", 
                                                       "high school graduate", 
                                                       "GED"),
                                              voctec  = c("some vocational/technical training", 
                                                          "completed vocational/technical training", 
                                                          "some college",
                                                          "some community college", 
                                                          "associate or junior college degree"),
                                              college = c("completed college (bachelor's degree)"),
                                              post = c("some post baccalaureate professional", 
                                                       "completed a post baccalaureate professional",
                                                       "completed post baccalaureate professional",
                                                       "completed a master's degree",
                                                       "completed a doctoral degree", 
                                                       "some graduate training beyond a master's degree",  
                                                       "some graduate school" )))) %>% 
    mutate_at(.vars = c("edu_rm1", "edu_rf1"), # collapse parental education categories
              .funs = funs(col = fct_collapse(., high = c("eighth grade or less.", 
                                                          ">eighth grade,  high school unfinished.",
                                                          " vocational school not high school.",
                                                          "high school graduate.",
                                                          "completed a GED.",
                                                          "He never went to school."),
                                              votec = c(" vocational school after high school.",
                                                        " college, not graduated."),
                                              college = c(" graduated from a college "),
                                              post = c("professional training beyond college")))) %>% 
    replace_with_na(list(edu_rm1_col = "school, but unknown level.",
                         edu_rf1_col = "school, but unknown level.")) %>%
    replace_with_na(list(edu_rm1_col = "schooling unknown.",
                         edu_rf1_col = "schooling unknown.")) %>%
    mutate_at(c("edu_rm1_col", "edu_rf1_col"), fct_drop) %>% 
    mutate(edu_par_max =  # maximum of ordinal vars. Strategy: convert to int, then do max.
             dplyr::select(., edu_rm1_col,  edu_rf1_col) %>% 
             mutate_all(as.integer) %>% 
             apply(1,max, na.rm = 1)) %>%  
    mutate(edu_par_max_parental_report_col =  # maximum of ordinal vars. Strategy: convert to int, then do max.
             dplyr::select(., parent_partner_report_col,  parent_self_report_col) %>% 
             mutate_all(as.integer) %>% 
             apply(1, max, na.rm = 1)) %>%  
    mutate(edu_max =  
             dplyr::select(., edu_4_col, edu_5_col) %>% 
             mutate_all(as.integer) %>% 
             apply(1, max, na.rm = 1 ))  %>%   
    replace_with_na(replace = list(edu_par_max= -Inf))  %>% # the max of two NAs in Inf 
    replace_with_na(replace = list(edu_par_max_parental_report_col = -Inf)) %>% 
    replace_with_na(replace = list(edu_max= -Inf)) %>% 
    mutate_at(c("edu_max", "edu_par_max", "edu_par_max_parental_report_col"), # convert back from int to fct
              ~  factor(.)  %>% 
                fct_recode("high" = "1",
                           "votec" = "2", 
                           "college" = "3", 
                           "post" = "4")) %>% 
    # https://stackoverflow.com/questions/13593703/replace-in-one-factor-by-another-factor-inside-dataframe 
    mutate(edu_p = 
             factor(if_else(is.na(edu_par_max_parental_report_col), 
                            as.character(edu_par_max), 
                            as.character(edu_par_max_parental_report_col)),
                    levels = c("high" ,
                               "votec" , 
                               "college", 
                               "post"))) %>% 
    mutate(debt = 
             debt %>% 
             factor %>% 
             fct_recode("net surplus" = "1",
                        "breaking even"  = "2",
                        "net debt" = "3")) %>% 
    mutate_at(c("income_4", "income_5"), as.character) %>%
    mutate_at(c("income_4", "income_5"), as.numeric) %>%
    mutate_at(c("income_4", "income_5"), factor, ordered = T)
  
  ########################################################
  # CHUNK 2: A COUPLE OF NEW VARIABLES (MICHAEL'S AFTERTHOUGHT)
  ########################################################
  
  # Job satisfaction h5lm26 
  # Personal earnings (use midpoints) h5ec1, (redundant, already defined above I call this income 5)
  # Household assets (h5ec2+4)-(h5ec5a+b+c)
  # Depressed mood scale based on H5SS0A-E.
  
  new_variables_for_michael = 
    waves_full %>%
    dplyr::select(AID, 
                  matches("H5SS0[A-E]"),
                  matches("H5EC[2-4]"), 
                  matches("H5EC5[A-C]"),
                  matches("H5EC[1-2]"))  %>% 
    mutate(H5EC4=ifelse(H5EC4==999999997,0,H5EC4)) %>% 
    replace_with_na_all(~.x %in% c(997,998)) %>% 
    # Depressed mood scale based on H5SS0A-E. (usually people just average these items)
    mutate(depressed_w5  = 
             waves_full %>% 
             dplyr::select(matches("H5SS0[A-E]")) %>% 
             apply(1,mean, na.rm = 1)) %>% 
    # Personal earnings (use midpoints) h5ec1, and assets h5ec2
    # mutate_at(c("H5EC1", "H5EC2"), 
    mutate_at(vars(matches("H5EC[1-2]")), 
              funs( . %>% 
                      factor %>% 
                      fct_recode("2500" = "1", 
                                 "7500" = "2", 
                                 "12500" = "3", 
                                 "17500"  = "4",  
                                 "22500"  = "5",  
                                 "27500" = "6",
                                 "35000" = "7", 
                                 "45000" = "8",
                                 "62500" = "9",
                                 "87500" = "10",
                                 "125000" = "11",
                                 "175000" = "12",
                                 "200000" = "13",
                                 NULL = "97",
                                 NULL = "98"
                      ) %>% 
                      as.character %>% 
                      as.numeric)) %>% 
    # more assets
    mutate_at(vars(matches("H5EC5[A-C]")), 
              funs( . %>% 
                      factor %>% 
                      fct_recode("0" = "1", 
                                 "5000" = "2", 
                                 "17500" = "3", 
                                 "37500"  = "4",  
                                 "75000"  = "5", 
                                 "175000" = "6",
                                 "375000" = "7", 
                                 "750000" = "8",
                                 "1000000" = "9"
                      ) %>% 
                      as.character %>% 
                      as.numeric)) %>% 
    mutate(#wave 5 household income, replace H5EC2 with H5EC1 when H5EC2 is missing(real missing and missing because of only 1 person in hh)
      income_hh_ff5=log(ifelse(is.na(H5EC2),H5EC1,H5EC2)), 
      # Household assets (h5ec2+c3?+4)-(h5ec5a+b+c) # JC indeed h5ec3 is not to be included
      assets_w5_pos           = dplyr::select(., matches("H5EC[2,4]")) %>% apply(1, sum, na.rm = 1),
      assets_w5_neg           = dplyr::select(., matches("H5EC5[A-C]")) %>% apply(1, sum, na.rm = 1),
      assets_household_net_w5 = assets_w5_pos - assets_w5_neg) %>% 
    # create household assets mike in the ses model and taking log
    mutate(asset_hh_ff5=case_when(!is.na(income_hh_ff5)&!is.na(H5EC4) ~ (income_hh_ff5+H5EC4),
                                  is.na(income_hh_ff5) ~ H5EC4,
                                  is.na(H5EC4) ~ income_hh_ff5), 
           asset_hh_ff5_log=ifelse(asset_hh_ff5>0, log(asset_hh_ff5), ifelse(asset_hh_ff5==0,0,NA)))
  waves = 
    waves %>% 
    left_join(new_variables_for_michael, by = "AID")
  
  ########################################################
  # CHUNK 4: ADD SHAWN'S INCOME-OCCUPATION DATA
  ########################################################
  
  
  waves = 
    waves %>% 
    mutate(SOC2000 = as.character(SOC2000)) %>% 
    left_join(shawn, by = "SOC2000")
  
  
  ########################################################
  # CHUNK 5: CECILIA
  ########################################################
  # wave5$H5LIFE5
  # NULL
  # Warning message:
  #   Unknown or uninitialised column: 'H5LIFE5'.
  ph <- 
    waves_full %>%
    dplyr::select(AID, 
                  # wave 5
                  H5EL5, H5PG5, H5EL6A,  H5EL6B, H5EL6C, H5EL6D, H5EL6E, H5EL6F, H5EL6G, H5EL6H,
                  H5EL6I, H5EL6J, H5EL6K, H5EL6L, H5EL6M, H5EL6N, H5EL6O, H5EL6P, H5EL6Q,
                  H5ID6AA, H5ID6BA, H5ID6CA, H5ID6DA, H5ID6EA, H5ID6FA, H5ID6GA, H5ID6HA, H5ID6IA, H5ID6JA,
                  H5ID6KA, H5ID6LA, H5ID6MA, H5ID6NA, H5ID6OA, H5ID6PA, H5ID6QA, H5ID6RA, H5ID6SA, H5ID6TA,
                  H5OD11, H5EC7, H5EC8, H5HR1, H5TR1, H5TR2, H5TR3, H5LM9, H5LM10, H5LM11, H5LM12, H5ID1,
                  H5TO2, H5TO4, H5TO10, H5TO13, H5TO14, H5TO21, H5TO26A, H5TO26B, H5TO26C, H5TO26D, H5TO27A,
                  H5ID6B, H5ID6C, H5ID6BM, H5ID6CM, H5ID6F, H5ID6FM, H5ID6Q,H5ID6QM, H5ID6A, H5ID6AM,
                  # wave 1
                  H5TO27B, H5TO27C, H5TO27D, AID, PA56 , PB8, PC31, PA10 , PA12, PA40, PA31, PA32, PA33, PA34, PA63, PA62,
                  PA22, PC37, PC38, PC53, PC54, PC49A_1, PC49B_1, PC49C_1, PC49D_1, PC49E_1, PC49F_1, H1GI9, H1RM4, H1RF4) %>% 
    #handling missing data 
    replace_with_na_all(condition = ~.x == 96) %>% 
    replace_with_na_at(.vars = c("PA63","PA56", "PA10", "PA31", "PA32", "PA33", "PA34"),# remove"PA62" from this line as PA62==6 means five or more times having more than 5 drinks
                       condition = ~ .x == 6) %>% 
    replace_with_na_at(.vars = c("H1RM4", "H1RF4"),
                       condition = ~ .x >= 98) %>% 
    #CHILDHOOD
    mutate(
      #education of the parent (HS/GED or less =1; 0 otherwise)
      ed_pa1 = as.factor(ifelse(is.na(PA12), NA, ifelse(PA12<= 5, 1, 0))),
      ed_pa2 = as.factor(ifelse(is.na(PB8), NA, ifelse(PB8<= 5, 1, 0))),
      #parental alcohol
      palchohol = as.factor(ifelse(is.na(PA62), NA, ifelse(PA62==1, 0, 1))),
      # parental marital status (1 if married, 0 otherwise)
      pmarital = as.factor(ifelse(is.na(PA10), NA, ifelse(PA10==2, 1, 0))),
      #parental marriage like relationships
      marriagelike = ifelse(is.na(PA40), NA, PA40),
      #paying the bills
      pbills = as.factor(ifelse(is.na(PA56), NA, ifelse(PA56==1,1, 0))),
      #educational aspirations (very disappointed =1 )
      peduasp = as.factor(ifelse(is.na(PC31), NA, ifelse(PC31==1, 1, 0))),
      #parental smoking
      psmoking = as.factor(ifelse(is.na(PA63) , NA, ifelse(PA63==1, 1, 0))),
      #parental occupation resident mother (note 97 is a legitimate skip meaning that the parent might be absent)
      poccrm = factor(H1RM4),
      #parental occupation resident father
      poccrf = factor(H1RF4),
      #neighborhood quality 1 (scale 1-5)
      pneighqual1= factor(PA32, ordered = TRUE),
      #neighborhood quality 2  (scale 1-5)
      pneighqual2= factor(PA33, ordered = TRUE),
      #neighboorhood quality 3 (scale 1-3)
      pneighqual3=factor(PA34, ordered = TRUE),
      #respondent
      #childhood health
      #disability
      cdisab = as.factor(ifelse(H5EL5>= 3 & (PC37==1 | PC38==1 | PC53==1 | PC54==1), 1, 0)),
      #health measure 1
      # JCQ (apparently _ not . on JC server!!)
      # chealth1 = as.factor(ifelse(H5EL5<= 3 & (PC49A.1 ==1 | PC49B.1==1 | PC49C.1==1 | PC49D.1==1 | PC49E.1==1 | PC49F.1==1), 1, 0)),
      chealth1 = as.factor(ifelse(H5EL5>= 3 & (PC49A_1 ==1 | PC49B_1==1 | PC49C_1==1 | PC49D_1==1 | PC49E_1==1 | PC49F_1==1), 1, 0)),
      # health measure 2
      chealth2 = as.factor(ifelse(H5EL5>= 3 & (H5EL6A==1 | H5EL6B==1 | H5EL6C==1 | H5EL6D==1 |
                                                 H5EL6E==1 |  H5EL6F==1 | H5EL6G==1 | H5EL6H==1 |
                                                 H5EL6I==1 | H5EL6L==1 | H5EL6M==1 | H5EL6N==1 | 
                                                 H5EL6O==1 | H5EL6P==1 | H5EL6Q==1), 1, 0)),
      #health measure 3
      #I excluded H5ID6QA and H5ID6RA since there were no observations. Moreover, different age threshold were applied
      #to some of the scores because the minimum age available was higher than 15.
      # In particular, age 20 for H5ID6EA, H5ID6KA, H5ID6LA, H5ID6MA, H5ID6SA
      #age 25 for H5ID6NA, H5ID6OA, H5ID6PA, and age 30 for H5ID6JA
      chealth3 = as.factor(ifelse(H5EL5>= 3 & (H5ID6AA <=15 | H5ID6BA <=15 | H5ID6CA <=15 | 
                                                 H5ID6DA <=15 | H5ID6EA <=20 | H5ID6FA <=15 | 
                                                 H5ID6GA <=15 | H5ID6HA <=15 | H5ID6IA <=15 | 
                                                 H5ID6JA <=30 | H5ID6KA <=20 | H5ID6LA <=20 | 
                                                 H5ID6MA <=20 | H5ID6NA <=25 | H5ID6OA <=25 | 
                                                 H5ID6PA <=25  | H5ID6SA <=20 | H5ID6TA <=15), 1, 0))) %>% 
    # ADULTHOOD
    mutate(
      ed = as.factor(ifelse(!is.na(H5OD11) , ifelse(H5OD11<= 4, 1, 0), NA)),
      bills = as.factor(ifelse(!is.na(H5EC7) , ifelse(H5EC7== 1, 1, 0), NA)),
      evict = as.factor(ifelse(!is.na(H5EC8) , ifelse(H5EC8== 1, 1, 0), NA)),
      married = as.factor(ifelse(!is.na(H5HR1) , ifelse(H5HR1== 1, 1, 0), NA)),
      #health problem between age 15 and 30 (for some variable only able to determine less than 30)
      ahealth2 = as.factor(ifelse(((H5ID6AA >15 & H5ID6AA <=30) |    # cancer
                                     (H5ID6BA>15 & H5ID6BA<=30) |    # high blood cholesteral
                                     (H5ID6CA >15 & H5ID6CA<=30) |   # hypertension
                                     (H5ID6DA>15 & H5ID6DA <=30) |   # high blood sugar or diabetes
                                     H5ID6EA <=30 |                  # heart attack
                                     (H5ID6FA >15 &  H5ID6FA <=30) | # asma
                                     (H5ID6GA >15 & H5ID6GA <=30) |  # depression
                                     (H5ID6HA >15 & H5ID6HA <=30) |  # ptsd
                                     (H5ID6IA >15 & H5ID6IA <=30) |  # anxiety
                                     H5ID6JA <=30 |                  # HIV
                                     H5ID6KA <=30 |                  # hepatitis  
                                     H5ID6LA <=30 |                  # chronic kidney disease
                                     H5ID6MA <=30 |                  # blood clot in lung
                                     H5ID6NA <=30 |                  # stroke
                                     H5ID6OA <=30 |                  # heart failure
                                     H5ID6PA <=30 |                  # atrial fibrilation
                                     H5ID6SA <=30 |                  # sleep apnea
                                     H5ID6TA <=30),                  # eating disorder
                                  1, 0)),
      #health problem age 29 + (for some variables greater than 30)
      # "97" is a legitimate skip, i.e. not had this diagnosis to date
      ahealth3 = as.factor(ifelse (((H5ID6AA >30 & H5ID6AA<96) |   # cancer
                                      (H5ID6BA>30 & H5ID6BA<96) |  # high blood cholesteral
                                      (H5ID6CA >30 & H5ID6CA <96) | # hypertension
                                      (H5ID6DA>30 & H5ID6DA<96) |  # high blood sugar or diabetes
                                      (H5ID6EA >30 & H5ID6EA<96) | # heart attack
                                      (H5ID6FA >30 & H5ID6FA<96) | # asthma
                                      (H5ID6GA >30 & H5ID6GA<96) | # depression
                                      (H5ID6HA >30 & H5ID6HA<96) | # ptsd
                                      (H5ID6IA >30 & H5ID6IA<96) | # anxiety
                                      (H5ID6JA >30 & H5ID6JA<96) | # HIV
                                      (H5ID6KA >30 & H5ID6KA<96) | # hepatitis
                                      (H5ID6LA >30 & H5ID6LA<96) | # chronic kidney disease
                                      (H5ID6MA >30 & H5ID6MA<96) | # blood clot in lung
                                      (H5ID6NA >30 & H5ID6NA<96) | # stroke
                                      (H5ID6OA >30 & H5ID6OA<96) | # heart failure
                                      (H5ID6PA >30 & H5ID6PA<96) | # atrial fibrilation
                                      (H5ID6SA >30 & H5ID6SA<96) | # sleep apnea
                                      (H5ID6TA >30 & H5ID6TA<96)), # eating disorder
                                   1, 0)),
      #own occupation at wave V (97 stands for no job, code 1-9)
      occw5=H5LM9,
      #self-reported health
      srh = as.factor(ifelse(H5ID1>=4, 1, 0)), # 0:excellent very good
      #smoker in the past 30 days more than 4 cigarettes
      tobacco = as.factor(ifelse((H5TO2>4 | H5TO4>4 | H5TO10>4), 1,0)),
      #alcohol (dummy =1 if number of days drank in the last 30 days >3 and average number of drinks per occasion>2)
      alcohol = as.factor(ifelse((H5TO13>2 & H5TO14>2), 1,0)), ##H5TO13 level changed 0 is none, 1 is one day
      #drug_use (dummy =1 if number of days smoked marijuana in the past 30 days or any misuse of Rx's medications or any misuse of hard drugs)
      drug = as.factor(ifelse((H5TO21>1 | H5TO26A==1 |  H5TO26B==1 | H5TO26C==1 |  H5TO26D==1 |
                                 H5TO27A==1 | H5TO27B==1 | H5TO27C==1  | H5TO27D==1), 1, 0)))
  
  waves = 
    waves %>% 
    left_join(ph, by = "AID")
  
  ########################################################
  # EXAMPLE SUMMARIES
  ########################################################
  
  #waves %>% 
  #  filter(AID %in% dat$AID) %>% 
  #  dplyr::select(ed_pa1, ed_pa2, palchohol, pmarital) %>% 
  #  summary
  
  ########################################################
  # CECILIAS CONTROL VARIABLES
  ########################################################
  
  waves_cecilia  =
    waves_full %>%
    dplyr::select(AID, H1GH59A, H1GH59B, H1GH60, H2WS16W, H2WS16HF,
                  H2WS16HI, H3WGT, H3HGT_F, H3HGT_I, H4BMI, H5ID3, H5ID2F, H5ID2I, H5ID6D, H5ID6E, H5ID6DM, H5ID6EM, PC19A_P, PC19B_O,
                  C_MED, C_HBA1C, C_JOINT, C_NFGLU, C_FGLU, H5TR4, Plate, AvgCorrelogram100, H5BMI, H5BMICLS)  %>% #biomarkers wave 5 and #quality control varialbes
    
    #handling missing data
    replace_with_na_at(.vars = c("H1GH59A","H1GH59B", "H3HGT_F","H3HGT_I", "H2WS16HF", "H2WS16HI", "PC19A_P", "PC19B_O"),
                       ~.x %in% c(96:99)) %>%
    replace_with_na_at(.vars = c("H5ID2I","H1GH60", "H3WGT", "H4BMI",  "H2WS16W"),
                       ~.x %in% c(996:999)) %>%
    replace_with_na_at(.vars = c("H5ID2F"),
                       ~.x == 98) %>%
    replace_with_na_at(.vars = c("H4BMI"),
                       ~.x %in% c(888:889)) %>%
    mutate(
      #bmi at different waves    
      #impute H3WGT to 330 if the code is 888 (which means 330 or more)
      H3WGT = ifelse(H3WGT==888, 330, H3WGT)  ,
      w1bmi = ifelse((is.na(H1GH60) | is.na(H1GH59A)  | is.na(H1GH59B)), NA, (H1GH60 / ((H1GH59A*12+ H1GH59B)*(H1GH59A*12+ H1GH59B)))*703) ,
      w2bmi = ifelse((is.na(H2WS16W) | is.na(H2WS16HF)  | is.na(H2WS16HI)), NA, (H2WS16W / ((H2WS16HF*12+ H2WS16HI)*(H2WS16HF*12+ H2WS16HI)))*703),
      w3bmi = ifelse((is.na(H3WGT) | is.na(H3HGT_F)  | is.na(H3HGT_I)), NA, (H3WGT / ((H3HGT_F*12+ H3HGT_I)*(H3HGT_F*12+ H3HGT_I)))*703) ,
      w4bmi = H4BMI ,
      w5bmi_selfreport = ifelse((is.na(H5ID3) | is.na(H5ID2F)  | is.na(H5ID2I)), NA, (H5ID3 / ((H5ID2F*12+ H5ID2I)*(H5ID2F*12+ H5ID2I)))*703),
      w5bmi = ifelse(!is.na(H5BMI),H5BMI,w5bmi_selfreport), 
      #birthweight
      birthweight = (PC19A_P*16 + PC19B_O),
      #obesity at different waves
      obesew5 = as.numeric(w5bmi>=30),
      obesew4 = as.numeric(w4bmi>=30),
      obesew3 = as.numeric(w3bmi>=30),
      obesew1or2 = ifelse(!is.na(w1bmi), as.numeric(w1bmi>=30), as.numeric(w2bmi>=30)),
      w1or2bmi  = ifelse(!is.na(w1bmi), w1bmi, w2bmi),
      #only overweight
      overweightw5 = as.numeric(w5bmi>=25 & w5bmi<=30),
      overweightw4 = as.numeric(w4bmi>=25 & w4bmi<=30),
      overweightw3 = as.numeric(w3bmi>=25 & w3bmi<=30),
      overweightw1or2 = ifelse(!is.na(w1bmi), as.numeric(w1bmi>=25 & w1bmi<=30), as.numeric(w2bmi>=25 & w2bmi<=30)),
      
      #overweight and obese
      overw5 = as.numeric(w5bmi>=25),
      overw4 = as.numeric(w4bmi>=25),
      overw3 = as.numeric(w3bmi>=25),
      overw1or2 = ifelse(!is.na(w1bmi), as.numeric(w1bmi>=25), as.numeric(w2bmi>=25)),
      
      #overweight at different waves
      overweight_ff5 = as.numeric(w5bmi>=25),
      
      # having or had diabetes and heart attack  
      diabetes = ifelse((H5ID6D==1 | H5ID6DM==1), 1, 0),
      heartatk = ifelse((H5ID6E==1 | H5ID6EM==1), 1, 0)
    ) 
  
  
  # adding to the dataset
  waves = waves %>% 
    left_join(waves_cecilia, by = "AID")
  
  
  ########################################################
  # Creating new VARIABLES Wenjia
  ########################################################
  
  #calc_age is a function to calculate the age based on the date of birth and the given reference date
  calc_age <- function(birthDate, refDate) {
    period <- as.period(interval(birthDate, refDate), unit="years")
    period$year
  }
  
  waves_wenjia =
    waves_full%>%
    dplyr::select(AID, H1GI1Y, H1GI1M, H2GI1Y, H2GI1M, H3OD1Y, H3OD1M, H4OD1Y, H4OD1M, H5OD1Y, H5OD1M,
                  IYEAR, IMONTH, IYEAR2, IMONTH2, IYEAR3, IMONTH3, IYEAR4, IMONTH4, IYEAR5, IMONTH5,
                  H1GI4, H3OD2, H5OD4C, H1GI6A, H3OD4A, H5OD4A, H1GI6B, H3OD4B, H5OD4B, H1GI6D,H3OD4D,
                  H5OD4D, H1GI9, H3IR4,H4IR4, W1REGION,W2REGION,W3REGION,W4REGION,REGION)%>% ###add region information for create variable for people in south
    dplyr::rename(W5REGION=REGION) %>% mutate_at(.var = "W5REGION", .funs = . %>%
                                            factor %>% fct_recode("W" = "1", #west
                                                                  "MW" = "2",  #midwest 
                                                                  
                                                                  "NE" = "3", # northeast
                                                                  
                                                                  "S" = "4", #south
                                                                  
                                                                  NULL = "98", #international
                                                                  NULL = "99")) %>%  # no location                                         
    
    #set birth month and year to na when the response refused to provide
    replace_with_na_at(c("H1GI1Y", "H1GI1M","H2GI1Y","H2GI1M"), ~.x%in% c("96", "98"))%>%
    replace_with_na_at(.vars = c("H1GI4","H3OD2","H5OD4C","H1GI6A", "H3OD4A","H5OD4A", "H1GI6B", "H3OD4B","H5OD4B",
                                 "H1GI6D","H3OD4D", "H5OD4D","H1GI9"), condition = ~.x %in% c(5, 6, 8, 9, 12, 13))%>%
    #replace na value in wave 1 with information in w2 or w3 w4 or w5
    mutate(BirthY=coalesce(as.integer(H1GI1Y),as.integer(H2GI1Y),as.integer(H3OD1Y),as.integer(H4OD1Y),as.integer(H5OD1Y)) %>% as.character,
           BirthY=ifelse(nchar(BirthY)==2, str_c("19",BirthY,""), BirthY),
           BirthM=coalesce(as.integer(H1GI1M),as.integer(H2GI1M),as.integer(H3OD1M),as.integer(H4OD1M),as.integer(H5OD1M)),
           birth_winter=ifelse(BirthM %in%c(1,2,3,12),1,0) %>% as.factor, #variable to indicate birth in winter
           data_winter_w1=ifelse(IMONTH %in%c(1,2,3,12),1,0) %>% as.factor, #variable to indicate data collected in winter
           data_winter_w2=ifelse(IMONTH2 %in%c(1,2,3,12),1,0) %>% as.factor,
           data_winter_w3=ifelse(IMONTH3 %in%c(1,2,3,12),1,0) %>% as.factor,
           data_winter_w4=ifelse(IMONTH4 %in%c(1,2,3,12),1,0) %>% as.factor,
           data_winter_w5=ifelse(IMONTH5 %in%c(1,2,3,12),1,0) %>% as.factor,
           
           birth_winter_exclsouth=ifelse(BirthM %in%c(1,2,3,12) & W1REGION!="S",1,0) %>% as.factor, #variable to indicate birth in winter people in south
           data_winter_w1_exclsouth=ifelse(IMONTH %in%c(1,2,3,12) & W1REGION!="S",1,0) %>% as.factor, #variable to indicate data collected in winter
           data_winter_w2_exclsouth=ifelse(IMONTH2 %in%c(1,2,3,12) & W2REGION!="S",1,0) %>% as.factor,
           data_winter_w3_exclsouth=ifelse(IMONTH3 %in%c(1,2,3,12) & W3REGION!="S",1,0) %>% as.factor,
           data_winter_w4_exclsouth=ifelse(IMONTH4 %in%c(1,2,3,12) & W4REGION!="S",1,0) %>% as.factor,
           data_winter_w5_exclsouth=ifelse(IMONTH5 %in%c(1,2,3,12) & W5REGION!="S",1,0) %>% as.factor,
           # reformate the time for conducting each wave xxxx-xx-xx(y-m-d)
           date_w1=paste(IYEAR, IMONTH,"1", sep = "-"),
           date_w2=paste(IYEAR2, IMONTH2,"1", sep = "-"),
           date_w3=paste(IYEAR3, IMONTH3,"1", sep = "-"),
           date_w4=paste(IYEAR4, IMONTH4,"1", sep = "-"),
           date_w5=paste(IYEAR5, IMONTH5,"1", sep = "-"),
           dob=paste (BirthY, BirthM, "1", sep = "-"),
           age_w1=calc_age(dob,date_w1),
           age_w2=calc_age(dob,date_w2),
           age_w3=calc_age(dob,date_w3),
           age_w4=calc_age(dob,date_w4),
           age_w5=calc_age(dob,date_w5),
           nh=coalesce(H1GI4,H3OD2, H5OD4C),#0:nonhispanic 1:hispanic
           rw=coalesce(H1GI6A, H3OD4A, H5OD4A, H1GI9, H3IR4, H4IR4), # white 1
           rb=coalesce(H1GI6B, H3OD4B, H5OD4B, H1GI9, H3IR4, H4IR4), # black: 1 or 2
           ra=coalesce(H1GI6D, H3OD4D, H5OD4D, H1GI9, H3IR4, H4IR4), # asian: 1 or 4
           #create variable re to reprsent the 5 combinations of race and enthnicity
           re=case_when(nh==0 & rw==1~1, #white nonhispanic
                        nh==0 & (rb==1|rb==2)  ~2, # black nonhispanic
                        nh==0 & (ra==1|ra==4) ~3, #asian nonhispanic
                        nh==0 & rw!=1 & rb!=1 & rb!=2 & ra!=4 & ra!=1 ~4, #other nonhispanic
                        nh==1 ~5) %>% as.factor)%>% #hispanic
    dummy_cols(select_columns =c( "re"))
  
  
  # #combine friends information from wave 1 and wave 2
  # friends=left_join(wave1_friends,wave2_friends,by="AID")%>%
  #   mutate_at(vars(-AID),funs(.%>%as.character%>%as.numeric))%>%
  #   replace_with_na_all(~.x %in% c(99999999, 88888888, 77777777, 55555555))%>%
  #   dplyr::select(-DAT)
  # 
  # #calculate the number of friends for each reponse(attention "" is counted in when data type is character)
  # 
  # 
  # friends$friends_num=apply(dplyr::select(friends,-AID), 1, function (x) n_distinct(x, na.rm = TRUE))
  # #replace duplicated nominated friends ids with NA
  # 
  # friends[]=t(apply(friends, 1, function(x) replace(x, duplicated(x), NA)))
  # friends_bmi=friends
  # #replace the friends AID with their bmi in wave 1 or wave2
  # friends_bmi[2:21] <- waves$w1or2bmi[match(as.character(unlist(friends_bmi[2:21])),
  #                                           as.character(waves$AID))]
  # 
  # friends_bmi=friends_bmi%>%mutate(avgbmi=apply(friends_bmi[2:21], 1, function (x) mean(x, na.rm = TRUE)),
  #                                  avgbmi_ff=apply(dplyr::select(.,starts_with("FF")), 1, function (x) mean(x, na.rm = TRUE)),
  #                                  avgbmi_mf=apply(dplyr::select(.,starts_with("MF")), 1, function (x) mean(x, na.rm = TRUE)))
  # 
  ########################################################birth weight and low birthweight
  
  
  waves_weight=waves_full%>%dplyr::select(AID, BIO_SEX, PC19B_O, PC19A_P, H5EL1P, H5EL1O, H5EL2, H4WAIST)%>%
    replace_with_na_all(~.x %in% c(98, 99, 998))%>%
    replace_with_na_at(.vars = "H5EL2", ~.x== 97)%>%
    replace_with_na_at(.vars = c("H4WAIST"), ~.x %in% c(996:999)) %>% 
    mutate(birthweight_P=coalesce(PC19A_P,H5EL1P),
           birthweight_O=coalesce(PC19B_O,H5EL1O), # suppliment birth weight info using self report info.
           birthweight_new=case_when(
             !is.na(birthweight_O) & !is.na(birthweight_P) ~ birthweight_P+birthweight_O/16,
             is.na(birthweight_O) & !is.na(birthweight_P) ~ birthweight_P),
           highbirthweight=case_when(
             birthweight_new < 8.8 ~ 0,
             birthweight_new >= 8.8 ~ 1, # birthweight high
             birthweight_P == 9 ~ 1 # for 9 pound birthweight, if birth ounce information is not avaiable 
           ),
           lowbirthweight=case_when(
             birthweight_new > 5.5 ~ 0, # birthweight more than 5.5 definitely not lowbirthweight
             birthweight_P!=5 & birthweight_new <= 5.5 ~ 1, # birthweight 0-4 pounds definitely lowbirthweigh
             !is.na(birthweight_O) & birthweight_P==5 ~ 1, # for 5 pound birthweight, if birth ounce information is avaiable, then it is less than 5.5, so it is lowbirthweight 
             is.na(birthweight_O) & birthweight_P==5 ~ H5EL2,# for 5 pound birthweight, if birth ounce information is na, use selfreport info
             is.na(birthweight_new) & (H5EL2==1|H5EL2==0) ~ H5EL2 #for birthweigh na case, use seflreport lowbirthweight information
           ),
           high_lowbirth= case_when( #either low or high birth weight
             (lowbirthweight==1 | highbirthweight==1) ~ 1,
             (lowbirthweight==0 | highbirthweight==0) ~ 0))
  
  
  
  
  ########################################################wave5 occupation
  # SEI_occ10=readRDS("/Volumes/Share/SEI_occ10.rds")
  
  
  w5occupation=transmute(waves_full,AID,H5LM12,H5LM27,H5LM5, H5LM22)%>%###H5LM22: occupation code for people who used to work in the past
    mutate_at(.vars = c("H5LM12","H5LM22"), as.character)%>%
    mutate(occ10=ifelse(H5LM12!=99997,case_when(
      H5LM12 %in% c("0010-3540","0010-0950","0010-0430")~".",
      H5LM12 %in% c("3600-3655","3600-4650")~ ".",
      H5LM12 %in% c("4700-4965","4700-5940")~ ".",
      H5LM12=="6005-7630"~ ".",
      H5LM12 %in% c("7700-8965","7700-9750")~".",
      nchar(H5LM12)==4~H5LM12
      
    ),case_when(
      H5LM22 %in% c("0010-3540","0010-0950","0010-0430")~".",
      H5LM22 %in% c("3600-3655","3600-4650")~ ".",
      H5LM22 %in% c("4700-4965","4700-5940")~ ".",
      H5LM22=="6005-7630"~ ".",
      H5LM22 %in% c("7700-8965","7700-9750")~".",
      H5LM22=="99997"~ ".",
      nchar(H5LM22)==4~H5LM22)))%>%
    mutate_at(.vars = "occ10", as.integer)%>%
    mutate(w5occupgroup=case_when(
      occ10 %in% c(10:3540)~"manbussciart", #management,business,science and art
      occ10 %in% c(3600:4650)~"service", #service
      occ10 %in% c(4700:5940)~"saleoffice", #sales and office
      occ10 %in% c(6005:7630)~"construction", #natural resources, construction and maintenance
      occ10 %in% c(7700:9750)~"production", #production, transportation, material moving
      occ10 %in% c(9800:9830)~"military"#military specific
    ))%>%
    mutate(w5occupgroup=if_else(is.na(w5occupgroup),case_when(
      H5LM27 %in% c(997,2,3) ~"undefined",# job category not defined, sick or maternity leave
      H5LM27 %in% c(1,4,5,6,7,9,10) ~"unemployed", #unemployed, disabled, student, retired, other
      H5LM27==8 ~"keepinghouse"#keeping house
    ),w5occupgroup)%>%as.factor) %>% 
    left_join(SEI_occ10 %>% dplyr::select(occ10,SEI10), by="occ10") %>% 
    mutate(SEI_ff5=SEI10)
  
  
  ######laura's substance variables and biological controls if pregnant
  
  #eversmoke 1, non eversmoke 0, current smoke 1, current not smoke 0
  #eversmoke is defined as smoke at least 1 cigaretter every day for 30 days
  #currentsmoke is defined in the past 30 days smoke at least one day
  #drinkeveryday is defined drink more than 3-5 days per week in the past 30 days
  #bingedrink is defined drink days in a row in the past 12 months.
  waves_laura=waves_full%>%dplyr::select(AID, BIO_SEX, H5TO1, H5TO2, H5TO13, H5TO15, H5PG5)%>%
    mutate(eversmoke=as.factor(H5TO1),
           currentsmoke=case_when(
             H5TO1 == 0 | H5TO2==0 ~ 0,
             H5TO2 >0 & !is.na(H5TO2) ~ 1
           ) %>% as.factor,
           
           ###drinkeveryday 1, not drink everyday 0
           drinkeveryday=case_when(
             H5TO13 >= 5 & H5TO13!=97 & !is.na(H5TO13) ~ 1,
             H5TO13 == 97 |  H5TO13 < 4 ~ 0
           )%>% as.factor,
           
           ###bingedrink 1, non bingerdrink 0
           
           bingedrink= case_when(
             H5TO15 ==0 | H5TO15 ==97 ~ 0,
             H5TO15 >0 & H5TO15 < 97 ~ 1
           )%>% as.factor,
           
           
           ## if is pregnant
           w5pregnant=case_when(
             H5PG5==0 ~0,
             H5PG5==1&BIO_SEX=="2" ~1,
             H5PG5==1&BIO_SEX=="1" ~0
           )%>% as.factor)
  #####wenjia's substance variables
  ##bingedrinker
  w5substanceuse=waves_full%>%dplyr::select(AID, BIO_SEX,H5TO1, H5TO13, H5TO14, H5TO15, H5TO20, H5TO21)%>%
    mutate(
      bingedrink_year=case_when(#use drinking information in the past 12 months
        H5TO15==0 | H5TO15 ==97 ~ "nonbinge",#non binge
        H5TO15==1 | H5TO15==2 ~ "occasionbinge", #occasional binge
        H5TO15==3  ~ "weeklybinge", #approaching weekly  binge
        H5TO15==4 |  H5TO15 ==5| H5TO15 ==6 ~ "frequentbinge" # frequent binge
        
      )%>%as.factor,
      bingedrink_month=case_when( #use drinking information in the past 30 days
        H5TO14==997 | H5TO13==97 ~ "nonbinge",# non binge
        BIO_SEX=="2" & H5TO14<4  ~ "nonbinge",# non binge
        BIO_SEX=="1" & H5TO14<5  ~ "nonbinge",# non binge
        BIO_SEX=="2" & H5TO14>=4 & H5TO13==1 ~ "occasionbinge",
        BIO_SEX=="1" & H5TO14>=5 & H5TO13==1 ~ "occasionbinge",
        
        BIO_SEX=="2" & H5TO14>=4 & H5TO13>1 ~ "regularbinge",#  binge
        BIO_SEX=="1" & H5TO14>=5 & H5TO13>1 ~ "regularbinge"#  binge
        
      )%>%as.factor,
      calchohol=case_when(# dichotomous variable corresponding to palchohol
        H5TO14==997 | H5TO13==97 ~ 0,
        H5TO14<5 ~ 0,
        H5TO14>=5 & H5TO13>0 ~ 1)
      %>%as.factor
    )
  
  ############
  # bingedrink_month_frequency=waves_full %>% filter( (BIO_SEX=="2" & H5TO14>=4)|(BIO_SEX=="1" & H5TO14>=5)) %>% right_join(aD[1], by="AID") %>% 
  # dplyr::select(BIO_SEX, H5TO13) %>% mutate(H5TO13=H5TO13 %>% as.factor %>% fct_recode("non"="1","one day"="2","2 or 3 days"="3","1 day a week"="4", "2 days a
  # week"="5", "3 to 5 days a week"="6", "every day or almost every day"="7", "refused"="97"))
  ############pubertal development for boys H1MP1-4 H2MP1-4, for girls H1FP1-4, H2FP1-4, H1FP6, H2FP9
  #############H1IR5, H2IR5 are judgement from the interviewer
  waves_pubertal=waves_full%>%dplyr::select(AID, BIO_SEX, matches("H1MP[1-4]"), matches("H2MP[1-4]"),
                                            H1FP1, H1FP2, H1FP3, H1FP4, H1FP6,
                                            H2FP1, H2FP2, H2FP3, H2FP4, H2FP9,
                                            H1IR5, H2IR5)%>%
    replace_with_na_at(.vars= c("H1MP1", "H1MP2","H1MP3","H1MP4",
                                "H2MP1", "H2MP2","H2MP3","H2MP4",
                                "H1FP1", "H1FP2","H1FP3","H1FP6",
                                "H2FP1", "H2FP2","H2FP3","H2FP9",
                                "H1IR5", "H2IR5"), ~.x%in% c(6,7,8,9))%>%
    replace_with_na_at(.vars = c("H1FP4", "H2FP4"), ~.x%in% c(96, 97, 98))%>%
    mutate_at(.vars = c("H1FP4", "H2FP4"), .funs = (. %>%
                                                      factor %>%
                                                      fct_collapse("1" = c("15","16","17"),
                                                                   "2" = c("13","14"),
                                                                   "3" = c("12"),
                                                                   "4" = c("11","10"),
                                                                   "5" = c("0","1","2","3","4","5","6","7","8","9")
                                                      )%>%
                                                      as.character%>%
                                                      as.numeric))%>%
    mutate(H2FP3=if_else(is.na(H2FP3) & H1FP3==1, H1FP3, H2FP3),
           H2FP4=if_else(is.na(H2FP4) & !is.na(H1FP4), H1FP4, H2FP4))
  
  waves_pubertal=waves_pubertal%>%
    mutate(w1pd=dplyr::select(., starts_with("H1"))%>%apply(1, mean, na.rm = TRUE),
           w2pd=dplyr::select(., starts_with("H2"))%>%apply(1, mean, na.rm = TRUE),
           w1pd=if_else(is.na(w1pd), H1IR5, w1pd),
           w2pd=if_else(is.na(w2pd), H2IR5, w2pd)
    )
  
  ############parents social economic status at wave 1 and 2
  
  waves_pses=waves_full %>% dplyr::select(AID, PA21, matches("PA57[A-F]"), PA59, PC21_6, PC22, PA20, PB17, PB18, H5EL3,
                                        H1RF9,H2RF9,H1RM9, H2RM9,
                                        H1RM4,H2RM4,H1RF4,H2RF4,
                                        H1RM14,H1RF14,H2RM14,H2RF14, PA63
  )%>%
    replace_with_na_at(.vars=c("H1RM4","H2RM4","H1RF4","H2RF4"), ~.x%in%c(96, 98, 99))%>%
    mutate(
      #occupation of parents
      poccrm_w12=coalesce(H1RM4,H2RM4)%>%
        factor %>%
        fct_collapse("manbussciart" = c("1","2","3","4"),
                     "service" = c("7"),
                     "saleoffice" = c("5","6"),
                     "construction" = c("9","10","14"),
                     "production" = c("8","11","12"),
                     "military" =c("13"),
                     "undefined"=c("15"),
                     "unemployed"=c("16"),
                     "norm"=c("97")),
      poccrf_w12=coalesce(H1RF4,H2RF4)%>%
        factor %>%
        fct_collapse("manbussciart" = c("1","2","3","4"),
                     "service" = c("7"),
                     "saleoffice" = c("5","6"),
                     "construction" = c("9","10","14"),
                     "production" = c("8","11","12"),
                     "military" =c("13"),
                     "undefined"=c("15"),
                     "unemployed"=c("16"),
                     "norf"=c("97")),
      #receive food stamp in one month before the interview in 1955
      pfoodstamp=case_when(PA57D==0 ~0,
                           PA57D==1 ~1)%>%as.factor,
      #other types of public assitance
      ppublicassit=case_when(
        PA21==1 | PA57B==1 | PA57C==1 | PA57D==1  | PA57F==1 ~ 1,
        PA57B==0 & PA57C==0 & PA57D==0  & PA57F==0 ~0,
        H1RF9==1 | H2RF9==1 | H1RM9==1 | H2RM9==1 ~1,
        H1RF9==0 & H2RF9==0 & H1RM9==0 & H2RM9==0 ~0)%>%as.factor,
      #children has health insurance or not
      chealthinsurance=case_when(
        PC21_6==0 ~1,
        PC21_6==1  ~ 0
      )%>%as.factor,
      #children has ever lost health insurance
      chealthinsurance_everlost=case_when(
        PC22==1 |PC22==7 ~ 1,
        PC22==0 ~0
      )%>%as.factor,
      #parents feel hard to get medical care for family
      phard_healthinsurance=case_when(
        PA59==3|PA59==4 ~ 1,
        PA59==1|PA59==2 ~ 0
      )%>%as.factor,
      #parent ever smoke
      peversmoke=case_when(
        PA63==1 ~1,
        H1RM14==1|H1RF14==1|H2RM14==1|H2RF14==1 ~1,
        (H1RM14==0 |H1RM14==7) & (H1RF14==0|H1RF14==7) & (H2RM14==0|H2RM14==7) & (H2RF14==0|H2RF14==7) ~0,
        PA63==0 ~0
      )%>%as.factor,
      #preterm delivery
      cpreterm=as.factor(H5EL3)
    )
  ############
  ####################Mike's added variable for long arm 
  waves_longarm_M=waves_full%>% 
    transmute(AID=AID,
              smoke_cf_p1=case_when(PC40==1 ~1,
                                    PC40==0 ~0) %>% as.factor,
              health_cf_p1=case_when(PC18==1 ~"excellent",
                                     PC18==2 ~"very good",
                                     PC18==3 ~"good",
                                     PC18==4 ~"fair",
                                     PC18==5 ~"poor")%>% as.factor,
              #occupation SEI of parents using w1 and w2 folks response
              temp1=coalesce(H1RM4,H2RM4)%>%
                as.factor,
              temp2=coalesce(H1RF4,H2RF4)%>%
                as.factor,
              
              
              SEI_rm_w12=case_when(temp1== "3"~ 40.22,
                                   temp1 %in% c("1","2") ~ 60.92,
                                   temp1== "4"~ 47.05,
                                   temp1== "6"~ 36.24,
                                   temp1== "5"~ 32.24,
                                   temp1== "7"~ mean(16.87,21.96),
                                   temp1== "13"~ 39.29,
                                   temp1== "14"~ 23.34,
                                   temp1 %in% c("8","9","10") ~ 31.51,
                                   temp1== "11"~ 22.58,
                                   temp1== "12"~ 26.50,
                                   temp1== "15"~ 36.81
              ) %>% as.numeric, #%>% replace_na(0),
              SEI_rf_w12=case_when(temp2== "3"~ 40.22,
                                   temp2 %in% c("1","2") ~ 60.92,
                                   temp2== "4"~ 47.05,
                                   temp2== "6"~ 36.24,
                                   temp2== "5"~ 32.24,
                                   temp2== "7"~ mean(16.87,21.96),
                                   temp2== "13"~ 39.29,
                                   temp2== "14"~ 23.34,
                                   temp2 %in% c("8","9","10") ~ 31.51,
                                   temp2== "11"~ 22.58,
                                   temp2== "12"~ 26.50,
                                   temp2== "15"~ 36.81
              ) %>% as.numeric,# %>% replace_na(0),
              # SEI_rf_w12=coalesce(H1RF4,H2RF4)%>%
              # as.factor %>% 
              #   fct_recode("40.22" = "3",
              #              "60.92" = c("1","2"),
              #              "47.05" = "4",
              #              "36.24" = "6",
              #              "32.24" = "5",
              #              min(16.87,21.96) ="7",
              #              "39.29"="13",
              #              "23.34"="14",
              #              "31.51"=c("8","9","10"),
              #              "22.58"="11",
              #              "26.50"="12",
              #              "36.81"="15"
              #   ) %>% as.numeric,
              # 
              SEI_max_p_w12=case_when(
                
                is.na(SEI_rm_w12) ~SEI_rf_w12,
                is.na(SEI_rf_w12) ~SEI_rm_w12,
                !is.na(SEI_rm_w12) & !is.na(SEI_rf_w12) ~ if_else(SEI_rm_w12>SEI_rf_w12, SEI_rm_w12,SEI_rf_w12)
                # is.na(SEI_rm_w12)|is.na(SEI_rf_w12) ~ if_else(is.na(SEI_rm_w12), SEI_rf_w12, SEI_rm_w12)
              ),
              SEI_max_p_w12_log=log(SEI_max_p_w12),
              income_pp1_log=case_when(PA55==0 ~0,
                                       PA55<9996 ~ log(PA55*1000)) ###unit is thousand dollars
              
    ) %>% dplyr::select(-starts_with("t"))
  
  ############blue collar and white collar
  waves_workcollar=waves_full%>% mutate_at(.vars = c("H5LM12","H5LM22"), as.character)%>%
    mutate(occ10=ifelse(H5LM12!=99997,case_when(
      H5LM12 %in% c("0010-3540","0010-0950","0010-0430")~".",
      H5LM12 %in% c("3600-3655","3600-4650")~ ".",
      H5LM12 %in% c("4700-4965","4700-5940")~ ".",
      H5LM12=="6005-7630"~ ".",
      H5LM12 %in% c("7700-8965","7700-9750")~".",
      nchar(H5LM12)==4~H5LM12
      
    ),case_when(
      H5LM22 %in% c("0010-3540","0010-0950","0010-0430")~".",
      H5LM22 %in% c("3600-3655","3600-4650")~ ".",
      H5LM22 %in% c("4700-4965","4700-5940")~ ".",
      H5LM22=="6005-7630"~ ".",
      H5LM22 %in% c("7700-8965","7700-9750")~".",
      H5LM22=="99997"~ ".",
      nchar(H5LM22)==4~H5LM22)))%>%
    mutate_at(.vars = "occ10", as.integer)%>%
    mutate(work_collar_ff5=case_when(
      occ10 %in% c(10:3540)~"white_collar", #management,business,science and art
      occ10 %in% c(3600:4650)~"service", #service
      occ10 %in% c(4700:5940)~"white_collar", #sales and office
      occ10 %in% c(6000:6130)~"farming", #natural resources, construction and maintenance
      occ10 %in% c(6200:9750)~"blue_collar", #production, transportation, material moving
      occ10 %in% c(9800:9830)~"military"#military specific
    ),
    #occupation SEI of parents using w1 and w2 folks response
    temp1=coalesce(H1RM4,H2RM4)%>%
      as.factor,
    temp2=coalesce(H1RF4,H2RF4)%>%
      as.factor,
    
    work_collar_rm_f12=case_when(temp1 %in% c("1","2","3","4","5","6") ~ "white_collar",
                                 temp1 %in% c("8","9","10","11","12")  ~ "blue_collar",
                                 temp1 %in% c("7","13") ~ "service",
                                 temp1 %in% c("14") ~ "farming"
    ) %>% as.factor, 
    work_collar_rf_f12=case_when(temp2 %in% c("1","2","3","4","5","6") ~ "white_collar",
                                 temp2 %in% c("8","9","10","11","12")  ~ "blue_collar",
                                 temp2 %in% c("7","13") ~ "service",
                                 temp2 %in% c("14") ~ "farming"
    ) %>% as.factor)%>% 
    dplyr::select(AID,work_collar_ff5,work_collar_rm_f12,work_collar_rf_f12)
  
  
  
  
  
  
  
  
  
  
  
  ############
  
  
  
  
  
  
  
  
  
  
  #################### birth order, number of big siblings, breast feeding time
  
  waves_birthorder =  
    waves_full %>% 
    dplyr::select(AID, H1HR14, H1HR15, matches("H1HR3[A-T]"), matches("H1HR5[A-T]"),matches("H1HR7[A-T]")) %>% 
    left_join(waves_wenjia %>% dplyr::select(AID,age_w1)) %>% 
    replace_with_na_all(~.x %in% c(996,997,998,999)) %>% 
    mutate(bio_sibling_ff1=ifelse(H1HR14<96,H1HR14,NA),
           birth_order_f1=case_when(H1HR14==1 ~ "single", #how many children biological parents have
                                    H1HR15==1 ~ "1st",
                                    H1HR15==2 ~ "2nd",
                                    H1HR15==3 ~ "3rd",
                                    H1HR15==4 ~ "4th",
                                    H1HR15==5 ~ "5th",
                                    H1HR15==6 ~ "6th",
                                    H1HR15==7 ~ "7th",
                                    H1HR15==8 ~ "8th",
                                    H1HR15==9 ~ "9th",
                                    H1HR15==10 ~ "10th",
                                    H1HR15==11 ~ "11th",
                                    H1HR15==12 ~ "12th",
                                    H1HR15==13 ~ "13th",
                                    H1HR15==14 ~ "14th",
                                    H1HR15==15 ~ "15th") %>% as.factor,
           birth_order_aggregate_f1=fct_collapse(birth_order_f1, "single"="single", "1st"="1st",
                                                 "2nd"="2nd", "3rd"="3rd", "4th"="4th",
                                                 "5th and more"=c("5th","6th","7th","8th","9th","10th","11th","12th","13th","14th","15th")))
  
  waves_birthorder_temp =  
    waves_birthorder %>%  
    gather(H1HR3A:H1HR7T, key= key, value = value) %>% 
    separate(key, sep = -1, into = c("var", "num")) %>% 
    split(.$var) %>%
    map(~mutate(., !!.$var[1] := paste0(var, num), !!paste0(.$var[1], "_values") := value)) %>%
    map(~dplyr::select(., -var, -value)) %>%
    Reduce(f = merge, x = .) %>%
    dplyr::select(-num) %>% 
    mutate(big_brothers_f1=case_when(H1HR14==1 | H1HR15==1 ~ 0, #single child; first child
                                     H1HR3_values==5 & H1HR5_values==1 & H1HR7_values>age_w1 ~ 1, # +1 big brother
                                     H1HR3_values==5 & H1HR5_values==1 & H1HR7_values<=age_w1 ~ 0, # no big brother
                                     H1HR3_values==5 & H1HR5_values!=1 ~ 0, # not full brother
                                     H1HR3_values<96 & H1HR3_values!=5 ~ 0
                                     
    ),
    big_sisters_f1=case_when(H1HR14==1 | H1HR15==1 ~ 0, #single child; first child
                             H1HR3_values==8 & H1HR5_values==7 & H1HR7_values>age_w1 ~ 1, # +1 big sis
                             H1HR3_values==8 & H1HR5_values==7 & H1HR7_values<=age_w1 ~ 0, # no big sis
                             H1HR3_values==8 & H1HR5_values!=7 ~ 0, # not full sis
                             H1HR3_values<96 & H1HR3_values!=7 ~ 0
                             
                             
                             
    )
    ) %>% 
    group_by(AID)%>%summarise_at(.vars=c("big_brothers_f1","big_sisters_f1"), .funs=funs(ifelse(all(is.na(.)), NA, sum(., na.rm=TRUE))))
  
  
  waves_birthorder=waves_birthorder %>% left_join(waves_birthorder_temp) %>%  
    mutate(big_sib_f1=big_brothers_f1+big_sisters_f1,
           big_brothers_aggregate_f1=case_when(big_brothers_f1 %in% c(0,1) ~ big_brothers_f1,
                                               big_brothers_f1>1 ~ 2),
           big_sisters_aggregate_f1=case_when(big_sisters_f1 %in% c(0,1) ~ big_sisters_f1,
                                              big_sisters_f1>1 ~ 2)    
    ) %>% 
    left_join(waves_full %>% transmute(AID=AID,
                                       breastfeed_length_fp1=PC20 %>% as.factor %>% fct_recode("less than 3 months"="1",
                                                                                               "3 to 6 months"="2",
                                                                                               "6 to 9 months"="3",
                                                                                               "9 to 12 months"="4",
                                                                                               "12 to 24 months"="5",
                                                                                               "more than 24 months"="6",
                                                                                               "not breastfed"="7",
                                                                                               NULL="96",
                                                                                               NULL="98"),
                                       breastfeed_length_aggregate_fp1=fct_collapse(breastfeed_length_fp1, "less than 3 months"="less than 3 months",
                                                                                    "3 to 6 months"="3 to 6 months",
                                                                                    "6 to 9 months"="6 to 9 months",
                                                                                    "9 to 12 months"="9 to 12 months",
                                                                                    "more than 12 months"=c("12 to 24 months","more than 24 months"),
                                                                                    "not breastfed"="not breastfed")   
                                       
    )
    )  %>% dplyr::select(AID, birth_order_f1, birth_order_aggregate_f1, big_brothers_f1, big_brothers_aggregate_f1,
                         big_sisters_f1, big_sisters_aggregate_f1, big_sib_f1,breastfeed_length_fp1, breastfeed_length_aggregate_fp1,bio_sibling_ff1)
  
  
  #######family pattern: two parents or single parent, number of siblings in the household (which is a little different from H1HR14
  ####### number of children biological parents have together)
  #family structure
  #1. number of children in the household
  #2. two parents (including two biological parents, one biological parent and his or her partner or spouse, adoptive parents)
  #   or single parent
  # H1HR3: relation; H1HR6: father or mother; H1HR9: whether or not always lived with this person (in the future in case needed to be refining)
  # waves_famstr =  
  #   waves_full %>% 
  #   dplyr::select(AID, H1HR14, H1HR15, matches("H1HR3[A-T]"), matches("H1HR6[A-T]"),matches("H1HR9[A-T]"))  %>% 
  #   gather(H1HR3A:H1HR9T, key= key, value = value) %>% 
  #   separate(key, sep = -1, into = c("var", "num")) %>% 
  #   split(.$var) %>%
  #   map(~mutate(., !!.$var[1] := paste0(var, num), !!paste0(.$var[1], "_values") := value)) %>%
  #   map(~dplyr::select(., -var, -value)) %>%
  #   Reduce(f = merge, x = .) %>%
  #   dplyr::select(-num) %>% 
  #   mutate(parents_num_f1=case_when(
  #     H1HR3_values %in% c(11:16) ~ 1 # +1 parents figure
  #   )   ) %>% 
  #   group_by(AID)%>%
  #   summarise_at(.vars=c("parents_num_f1"), .funs=funs(ifelse(all(is.na(.)), NA, sum(., na.rm=TRUE)))) %>% 
  #   mutate(famstr_ff1=case_when(parents_num_f1>1 ~ "two_parents",# any type
  #                               parents_num_f1==1 ~ "single_parent",
  #                               is.na(parents_num_f1) ~ "other"))#no parents figure in the household
  ###use add health calculated family structure data instead, they have more detailed calculation for family structure, 
  
  
  
  ####################
  #children maltreatment
  
  waves_maltreat=waves_full %>% dplyr::select(AID,matches("H.MA."), "H4WP3", "H4WP12", "H4CJ17", "H4RD19", "H4SE32", "H4DS14", "H4DS16", "PA10") %>% dplyr::select(-18) %>% 
    replace_with_na_all(~.x %in% c(95:99)) %>% 
    mutate_at(.vars=c("H3MA1","H3MA2","H3MA3","H3MA4","H4MA1", "H4MA3", "H4MA5"), ~ ifelse(.x==6,0,.x)) %>%
    #replace never happened to 0 instead of 6
    
    mutate(supervision_neglect_f3=H3MA1,
           physical_neglect_f3=H3MA2,
           physical_abuse_f3=H3MA3,
           sexual_abuse_f3=H3MA4,
           social_service_investigate_f3=H3MA5,
           social_service_intervene_f3=H3MA6,
           emotional_abuse_f4=H4MA1,
           emotional_abuse_age_f4=H4MA2,
           physical_abuse_f4=H4MA3,
           physical_abuse_age_f4=H4MA4,
           sexual_abuse_f4=H4MA5,
           sexual_abuse_age_f4=H4MA6,
           mother_jail_f4=H4WP3,
           father_jail_f4=H4WP12,
           jail_f4=H4CJ17,
           partner_abuse_f4=H4RD19,
           sexual_abuse_ever_f4=H4SE32,
           community_violence_f4=H4DS14,
           community_violence_f4=H4DS16,
           parental_marit_p1=PA10
    ) %>% 
    mutate_at(.vars = c("supervision_neglect_f3",
                        "physical_neglect_f3",
                        "physical_abuse_f3",
                        "sexual_abuse_f3",
                        "social_service_investigate_f3",
                        "social_service_intervene_f3",
                        "emotional_abuse_f4",
                        "physical_abuse_f4",
                        "sexual_abuse_f4"),function(.x) ifelse(!is.na(.x),ifelse(.x==0,0,1),.x)) %>%
    #coded in 0 (never happened) and 1(happened) 
    mutate(physical_abuse_f34=ifelse(is.na(physical_abuse_f3),physical_abuse_f4,physical_abuse_f3),
           sexual_abuse_f34=ifelse(is.na(sexual_abuse_f3),sexual_abuse_f4,sexual_abuse_f3)) %>% 
    # physical and sex abuse are both measured in wave 3 and 4, use data in wave 4 to suppliment wave 3
    mutate(
      maltreat_cum_index=dplyr::select(., .vars = c("supervision_neglect_f3",
                                                    "physical_neglect_f3",
                                                    "physical_abuse_f34",
                                                    "sexual_abuse_f34",
                                                    "social_service_investigate_f3",
                                                    "social_service_intervene_f3",
                                                    "emotional_abuse_f4")) %>%
        apply(1, function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))), #create an accumulative index of 7 index
      maltreat_cum_index_aggregate=case_when(maltreat_cum_index %in% c(0,1) ~ maltreat_cum_index, 
                                             # recode the index to 0(none of the abuse ever happened, 1 (one of those happened),
                                             #2 (2 and more than 2 of those happened))
                                             maltreat_cum_index>1 ~ 2)) %>% 
    dplyr::select(-starts_with("H"))
  ######add folk family situation at wave 5
  #the number of children add health folk has (including adpoted children)
  #####add folk relationship mark
  waves_offsprings=waves_full %>% dplyr::select(AID, H5PG2,H5PG5, H5PG7, H5PG30,H5PG31,H5PG32,H5PG33,H5PG34A) %>% mutate(
    num_child_ff5=case_when(
      H5PG7 %in%c(0:12) ~H5PG7,
      H5PG7==997 ~ 0,
      H5PG2==3 ~0,
      H5PG34A %in% c(1:5) ~1 # if they answer this question means at least they have 1 child
      
    )# number of children add health folk has at wave 5
  )
  
  waves_famstr_ff5=waves_full %>% dplyr::select(AID,H5HR1:H5HR6O)%>% 
    gather(H5HR4A:H5HR6O, key= key, value = value) %>% 
    separate(key, sep = 5, into = c("var", "num")) %>% 
    split(.$var) %>%
    map(~mutate(., !!.$var[1] := paste0(var, num), !!paste0(.$var[1], "_values") := value)) %>%
    map(~dplyr::select(., -var, -value)) %>%
    Reduce(f = merge, x = .) %>%
    dplyr::select(-num) %>% 
    mutate(partner=case_when(
      H5HR6_values %in% c(1:2) ~ 1, # +1 partner or spouse
      H5HR6_values >2 ~0 
    )   ) %>% group_by(AID)%>%
    summarise_at(.vars=c("partner"), .funs=funs(ifelse(all(is.na(.)), NA, sum(., na.rm=TRUE)))) %>% 
    mutate(famstr_ff5=case_when(partner >=1 ~ "partnership",# any type
                                partner==0 ~ "single",#no partner living together
                                is.na(partner) ~ "other"))#other                               
  
  waves_physact=waves_full %>% dplyr::select(AID, H5ID25, H5ID26, H5ID27,H5ID29)%>%
    ##another way to define physical activity
    # mutate_at(vars(matches("H5ID2")), .funs = funs(case_when(
    #   . ==0 ~ 0,
    #   . %in% c(1,2) ~ 1,
    #   . >2 ~ 2))) %>%
    # mutate(
    #   temp= dplyr::select(.,matches("H5ID2")) %>% apply(1,max, na.rm = T)) %>%
    # mutate(
    #   phys_activ=case_when(temp ==0 ~ "none",
    #                        temp ==1 ~ "moderate",
    #                        temp ==2 ~ "intens"))
  
  mutate(temp=dplyr::select(.,matches("H5ID2")) %>% apply(1,  function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))),
         temp2=temp/4) %>% 
    mutate(
      phys_activ_ff5=case_when(temp2 ==0 ~ "none",
                               temp2<3 ~ "moderate",
                               temp2 >=3 ~ "intens"))
  
  
  #############social integration w5 and cvd medication taken at w4
  
  waves_t1=waves_full %>% dplyr::select(AID, H5SS1,H5SS2,H5ID9)%>%
    
    mutate(#H5SS1_rescale=H5SS1 %>% scales::rescale(to=c(1,5)),
      H5SS1_sd=H5SS1 %>% scale(),
      H5SS2_sd=H5SS2 %>% scale(),
      integration_ff5=case_when(!is.na(H5SS1) & !is.na(H5SS2) ~ as.numeric(H5SS1_sd) + as.numeric(H5SS2_sd),
                                is.na(H5SS1) ~  as.numeric(H5SS2_sd),
                                is.na(H5SS2) ~ as.numeric(H5SS1_sd)),
      
      ins_ff5=case_when(H5ID9 %in% c(1:7,11:13) ~ "private_ins",
                        H5ID9 %in%c(8:10) ~ "aid_ins",
                        H5ID9 ==15 ~ "no_ins"
                        
      ),
      integration_quantile_ff5=integration_ff5 %>% percent_rank() %>%
        cut(c(0, .25, .5, .75, 1.01), right = FALSE, labels=c("high_iso","mid_iso","mid_int","high_int")),
      
      isolation_ff5=case_when(H5SS1 %in% c(1:3) & H5SS2 ==1 ~ "total_isolation",
                              H5SS1 %in% c(1:3) & H5SS2 %in% c(2:5) ~ "friends_iso",
                              H5SS1 %in% c(4:7) & H5SS2 ==1 ~ "neigbor_iso",
                              H5SS1 %in% c(4:7) & H5SS2%in% c(2:5) ==1 ~ "no_isolation")
      
      
    )
  
  w4meds=w4meds %>% mutate(
    cvd_med_taken_ff4=ifelse(SET1%in%cvd_med$Var1,1,0)) %>%
    group_by(AID) %>% summarise_at(.vars =vars(cvd_med_taken_ff4), .funs = funs(sum))
  
  
  ######################################################## create waves_wenjia by combine all the newly created variables
  waves_wenjia=waves_wenjia%>%
    # left_join(friends_bmi[c("AID","friends_num","avgbmi", "avgbmi_ff", "avgbmi_mf")], by= "AID")%>%
    left_join(obesityclass[c("AID", "obesityclass", "Prob_C1", "Prob_C2", "Prob_C3")], by= "AID")%>%
    left_join(waves_weight[c("AID", "birthweight_new", "lowbirthweight", "H4WAIST", "highbirthweight","high_lowbirth")], by="AID")%>%
    left_join(dplyr::select(w5occupation, AID, w5occupgroup, SEI_ff5), by="AID")%>%
    left_join(dplyr::select(waves_pubertal,AID, w1pd, w2pd), by="AID")%>%
    left_join(dplyr::select(w5substanceuse,AID, bingedrink_year, bingedrink_month, calchohol), by="AID")%>%
    left_join(dplyr::select(waves_pses, AID, poccrm_w12,poccrf_w12,pfoodstamp,ppublicassit,
                            chealthinsurance,chealthinsurance_everlost, phard_healthinsurance,peversmoke,cpreterm), by="AID") %>% 
    left_join(waves_longarm_M,by="AID") %>% 
    left_join(waves_birthorder,by="AID") %>% 
    left_join(waves_maltreat,by="AID") %>% 
    left_join(famstr_w1, by="AID") %>% 
    left_join(waves_offsprings %>% dplyr::select(AID, num_child_ff5), by="AID") %>% 
    left_join(waves_famstr_ff5 %>% dplyr::select(AID, famstr_ff5), by="AID") %>% 
    left_join(waves_physact %>% dplyr::select(AID, phys_activ_ff5), by="AID") %>% 
    left_join(waves_workcollar, by="AID") %>% 
    left_join(waves_t1, by="AID") %>%
    left_join(w4meds, by="AID") 
   
  
  waves=waves%>%left_join(waves_wenjia, by="AID")%>%
    left_join(dplyr::select(waves_laura, -starts_with(match = "H5")),by="AID")%>%
    dummy_cols(select_columns =c( "edu_max", "sex_interv", "obesityclass", "w5occupgroup","W5REGION","bingedrink_year","bingedrink_month","Plate")) %>% 
    left_join(PGS,by="AID")
  
  
  ######################################################## put new variables in to waves and make some of them dummy
  ########## derived variables (when using variables from multiple chuncks)
  waves=waves%>% mutate(
    age_w1orw2= ifelse(!is.na(w1bmi), age_w1, age_w2))
  
  
  ######################## create SES composite, interaction
  
  
  waves_ses_composite=waves %>% mutate_at(.vars = c("edu_max","edu_p"),
                                          .funs = list(~ .x %>% factor %>%
                                                         fct_recode("1" = "high", 
                                                                    "2" = "votec",
                                                                    "3" = "college",
                                                                    "4" = "post") %>% 
                                                         as.numeric)) %>% 
    mutate_at(.vars = c("sss_5","edu_max","income_hh_ff5","SEI_ff5","SEI_max_p_w12","edu_p","income_pp1_log"), 
              .funs = list(~ .x %>% scale(center = TRUE, scale = TRUE)%>% as.numeric)) %>% 
    mutate(ses_sss_composite = (sss_5 + SEI_ff5 + edu_max + income_hh_ff5),
           ses_composite_ff5 = (SEI_ff5 + edu_max + income_hh_ff5),
           ses_composite_pp1 = (SEI_max_p_w12 + edu_p + income_pp1_log),
           ses_interaction_1 = ses_sss_composite * ses_composite_pp1,
           ses_interaction_2 = ses_composite_ff5 * ses_composite_pp1) %>% 
    dplyr::select(AID, ses_sss_composite, ses_composite_ff5, ses_composite_pp1, ses_interaction_1, ses_interaction_2)
  
  
  
  
  ### funs(name = f(.)) ---->  list(name = ~f(.))
  #############wave 5 bio marker controls and ancestry
  ######################
  
  
  
  
  w5biocovars_wave <- w5biocovars %>% left_join(waves %>% dplyr::select(AID, sex_interv_m), by="AID")
  
  w5biocovars_wave <- w5biocovars_wave %>% mutate(pregnant_biow5     = case_when (sex_interv_m==1 ~ 0, #male by default not pregnant, don't know are assigned as missing
                                                                                  sex_interv_m==0 & Q012==2  ~ 0, #female and not pregnant
                                                                                  sex_interv_m==0 & Q012==1  ~ 1) %>% as.factor, #female and pregnant 
                                                  illness_4wks_biow5 = case_when (Q013a==1 | Q013b==1 | Q013c==1 | Q013d==1 | Q013e==1 | Q013f==1  ~ 1, #any illness in the past 4 weeks
                                                                                  Q013a>=2 | Q013b>=2 | Q013c>=2 | Q013d>=2 | Q013e>=2 | Q013f>=2  ~ 0)%>% as.factor,
                                                  illness_2wks_biow5 = case_when (Q014a==1 | Q014b==1 | Q014c==1 | Q014d==1 | Q014e==1 | Q014f==1  ~ 1, #any illness in the past 4 weeks
                                                                                  Q014a>=2 | Q014b>=2 | Q014c>=2 | Q014d>=2 | Q014e>=2 | Q014f>=2  ~ 0)%>% as.factor,
                                                  smoking_biow5     = case_when (Q015==1  ~ 1, #smoking at the time of interview
                                                                                 Q015==2  ~ 0)%>% as.factor,
                                                  time_biow5        =  EXAMDATE, #if we want to include as linear predictor, otherwise we should create categories
                                                  kit_biow5         = case_when (KITCOND==0 ~ 0, # kit ok
                                                                                 KITCOND>0 ~ 1)%>% as.factor, #some problems with the kit
                                                  tube_biow5        = case_when (TUBECOND==0 ~ 0, # tube ok
                                                                                 TUBECOND>0 ~ 1)%>% as.factor, #some problems with the tube
                                                  more48h_biow5     = case_when (AlqQuality=="> 48 hours" ~ 1, #something that went wrong with timing (ask Brandt, no info in the documentation)
                                                                                 AlqQuality=="Normal" ~ 0)%>% as.factor,
                                                  months_biow5      = as.factor (month(EXAMDATE)),
                                                  hour_biow5        = as.factor (hour(EXAMTIME)),
                                                  travel_biow5      = case_when (Q016==1  ~ 1,
                                                                                 Q016>1   ~ 0)%>% as.factor
  )
  
  w5biocovars_wave$AID<-as.character(w5biocovars_wave$AID)
  
  w5biocovars_constr <- w5biocovars_wave %>%dplyr::select(pregnant_biow5, illness_4wks_biow5, illness_2wks_biow5, smoking_biow5, kit_biow5, 
                                                          tube_biow5, more48h_biow5, FastHrs, months_biow5, hour_biow5, travel_biow5, AID)
  
  
  
  #change aid label for the coherence
  ancestralPCA <- ancestralPCA  %>% dplyr::rename("AID" ="aid")
  
  ancestralPCA <-ancestralPCA %>% dplyr::rename_at(vars( starts_with("X")), list( ~(str_replace(., "X", "AncestryPC") ) ))
  
  #change aid label for the coherence and create a dummy for biological relatedness
  biorel <- biorel %>% mutate(biorels = case_when(sample=="AddHealthUnrel" ~ 0,
                                                  sample=="RelsAddHealth" ~ 1)) %>% 
    dplyr:: rename("AID" ="aid")
  
  ancestry_biorel <- ancestralPCA %>% left_join(biorel, by="AID") 
  ancestry_biorel$AID<-as.character(ancestry_biorel$AID)
  
  
  
  waves=waves %>% left_join(w5biocovars_constr,by="AID") %>% left_join(ancestry_biorel,by="AID") %>% left_join(waves_ses_composite, by="AID")
  
  
  ########################
  ########################
  ##create wave 1 neighbourhood varialbes using wave 1 contextual data
  
  ####doing pca for these variable could be done as
  ###
  # dt=w1neighbour %>% dplyr::select(-1)
  # pca=FactoMineR::PCA(dt%>% na.omit())
  # barplot(pca$eig[,1],main="Eigenvalues",names.arg=1:nrow(pca$eig))
  # #get rotated/transformed data
  # ind=get_pca_ind(pca)
  # new_dt=ind$coord
  
  
  
  ###
  w1neighbour = waves_full %>% transmute(
    AID=AID, #variable initial means: S = State, C = County, T = Tract, B = Block Group
    poverty1989_w1 = BST90598, #Proportion Persons with Income in 1989 Below Poverty Level
    unemployed_w1=BST90754,    #Unemployment Rate
    welfare_receive_w1 = SHC95A31, #Proportion Adult Medicaid Population Receiving AFDC - State
    female_headed_hh_w1 = BST90455, #Proportion Households that are Female Householder, No Husband Present, Households with Own Children Under 18 Years
    nohighschool_25plus_w1 = BST90680, #Proportion Aged 25+ with No High School Diploma or Equivalency
    college_25plus_w1 = BST90686, #Proportion Aged 25+ with College Degree or More
    managerial_prof_employ_w1 = BST90795 #Proportion Employed in Managerial and Professional Specialty Occupations
  ) %>% mutate(
    neighbourhood_composite_w1 = ###as all the columns is percentage, i.e. on the same scale, i did not scale them before pca.
      (prcomp(~ poverty1989_w1 + unemployed_w1 + female_headed_hh_w1 + nohighschool_25plus_w1 + college_25plus_w1 + managerial_prof_employ_w1,
              data = ., na.action = na.exclude)$x)[,1]) #use summary of the pca result could get the proportion of variance: 1st pc 70%, first 2 86% 
  
  
  waves=waves %>% left_join(w1neighbour)
  
  ########################
  ########################
  
  # depression scales for Justin's paper possilbe example.
  # the method to create these variables from
  #Chen, P.; & Harris, K. M. (2019). Association of Positive Family Relationships With Mental Health Trajectories From Adolescence To Midlife.
  # We summed the responses on the 3 CES-D items (range, 0-9, with 0 indicating never having any depressive symptoms and
  #9 indicating having symptoms most or all of the time) at each wave, with a higher score representing greater depressive symptoms
  #(1) could not shake off the blues,evenwith help fromfamily and friends; H1FS3, H2FS3, H3SP6, H4MH19, H5SS0A
  #(2) felt depressed during the past 7 days;H1FS6, H2FS6, H3SP9, H4MH22, H5SS0B
  #and (3) felt sad, H1FS16, H2FS16, H3SP12, H4MH26, H5SS0D
  depression_CESD= waves_full %>% dplyr::select(AID, H1FS3, H2FS3, H3SP6, H4MH19, H5SS0A,
                                                H1FS6, H2FS6, H3SP9, H4MH22, H5SS0B,
                                                H1FS16, H2FS16, H3SP12, H4MH26, H5SS0D) %>%
    replace_with_na_at(.vars = vars(-matches("AID")), ~.x %in% c("6", "8" ,"9")) %>% 
    mutate_at(.vars= vars(matches("H5")), .funs = list(~ (as.numeric(.x) - 1))) %>% 
    mutate_at(.vars = vars(-matches("AID")),.funs = list(~ (as.numeric(.x)))) %>% 
    mutate(CESD_w1 = (as.numeric(H1FS3) + as.numeric(H1FS6) + as.numeric(H1FS16)),
           CESD_w2 = (as.numeric(H2FS3) + as.numeric(H2FS6) + as.numeric(H2FS16)),
           CESD_w3 = (as.numeric(H3SP6) + as.numeric(H3SP9) + as.numeric(H3SP12)),
           CESD_w4 = (as.numeric(H4MH19) + as.numeric(H4MH22) + as.numeric(H4MH26)),
           CESD_w5 = (as.numeric(H5SS0A) + as.numeric(H5SS0B) + as.numeric(H5SS0D)),
           CESD_composite_w1 = CESD_w1 %>% scale()%>% as.numeric,
           CESD_composite_w2 = CESD_w2 %>% scale()%>% as.numeric,
           CESD_composite_w3 = CESD_w3 %>% scale()%>% as.numeric,
           CESD_composite_w4 = CESD_w4 %>% scale()%>% as.numeric,
           CESD_composite_w5 = CESD_w5 %>% scale()%>% as.numeric) %>% dplyr::select(- starts_with("H"))
  
  
  physical_attractive = waves_full %>% dplyr::select(AID, "H1IR1","H2IR1","H3IR1","H4IR1") %>% 
    replace_with_na_at(.vars = c("H1IR1","H2IR1","H3IR1","H4IR1"),
                       condition = ~ .x %in% c(6,7,8,9)) %>% 
    mutate(attractive_w1 = H1IR1,
           attractive_w2 = H2IR1,
           attractive_w3 = H3IR1,
           attractive_w4 = H4IR1) %>% dplyr::select(- starts_with("H"))
  
  
  
  
  waves=waves %>% left_join(depression_CESD, by = "AID")%>% left_join(physical_attractive, by = "AID")
  
  # cibersort
  waves = waves %>%  
    mutate(time_biow5 = case_when( (hour_biow5==6 | hour_biow5==7) ~ "6_7",
                                   (hour_biow5==8 | hour_biow5==9) ~ "8_9",
                                   (hour_biow5==10 | hour_biow5==11) ~ "10_11",
                                   (hour_biow5==12 | hour_biow5==13) ~ "12_13",
                                   (hour_biow5==14 | hour_biow5==15) ~ "14_15",
                                   (hour_biow5==16 | hour_biow5==17) ~ "16_17",
                                   (hour_biow5==18 | hour_biow5==19 | hour_biow5==20) ~ "18_19_20"))
  
  cibersort               = read.csv(str_c(data_input,"two_batches_cibersort/CIBERSORT.Output_Job17.csv"), header = TRUE, sep = ",")
  cibersort_crosswalk     = read.csv(str_c(data_input,"two_batches/gene.expression.cpm.log2.040120/addhealth.cibersort.crosswalk.032720.txt"), header = TRUE, sep = "\t")
  cibersort_crosswalk$aid = as.character(cibersort_crosswalk$aid)
  cibersort_crosswalk     = cibersort_crosswalk %>% dplyr::rename(AID=aid)
  cibersort_final         = cibersort %>%  dplyr::rename(sample=Input.Sample) %>% left_join(cibersort_crosswalk, by="sample")
  waves                   = waves %>% left_join(cibersort_final, by = "AID") 
  
  # potential mediators for SES paper
  waves_mediators = waves_full %>%
    dplyr::select(AID, H5MN1, H5MN2, H5MN3, H5MN4, H5ID10A, H5ID10B, H5ID10C) %>% 
    dplyr::mutate_at(.vars = c("H5MN2", "H5MN3"),
              .funs = list(~ 6 - .x)) %>% 
    dplyr::mutate(stress_perceived = dplyr::select(., H5MN1, H5MN2, H5MN3, H5MN4,) %>% apply(1, mean, na.rm = TRUE),
           insurance_lack = ifelse(H5ID10A ==1| H5ID10B ==1| H5ID10C ==1, 1, 0) %>% as.factor)
  waves                   = waves %>% left_join(waves_mediators, by = "AID") 
  ########################################################
  # add CP's ACE
  ########################################################
  waves_ace=waves_full %>% dplyr:: transmute(AID = AID,
                                             ace_physical_neglet_br=case_when(H3MA2==6~0,
                                                                              H3MA2==1~0,
                                                                              H3MA2>1 & H3MA2<6 ~1),
                                             ace_foster_br=case_when(H3OD31==0~0,
                                                                     H3OD31==1~1),
                                             ace_sexual_br=case_when(H3MA4>0 & H3MA4<6 & H4MA5>0 & H4MA5<6 ~1,
                                                                     H4MA5>0 & H4MA5<6 ~1,
                                                                     H3MA4==6 & H4MA5==6~0,
                                                                     H3MA4>6 & H4MA5==6~0,
                                                                     H3MA4==6 & H4MA5>6~0),
                                             ace_physical_abuse_br=case_when(H4MA3>1 & H4MA3<6 ~1,
                                                                             H4MA3==6~0,
                                                                             H4MA3==1~0),
                                             ace_alcoholism_br=case_when(PC49E_2==1 ~1,
                                                                         PC49E_3==1 ~1,
                                                                         PC49E_2==0 ~0,
                                                                         PC49E_3==0 ~0),
                                             BirthY1=coalesce(as.integer(H1GI1Y),as.integer(H2GI1Y),as.integer(H3OD1Y),as.integer(H4OD1Y),as.integer(H5OD1Y)),
                                             BirthM1=coalesce(as.integer(H1GI1M),as.integer(H2GI1M),as.integer(H3OD1M),as.integer(H4OD1M),as.integer(H5OD1M)),
                                             dates_w1=paste(IYEAR, IMONTH,"1", sep = "-"),
                                             dob1=paste (BirthY1, BirthM1, "1", sep = "-"),
                                             ages_w1=calc_age(dob1,dates_w1),
                                             ace_jail_br=case_when(H4WP9 ==1 & H4WP11<=ages_w1~1,
                                                                   H4WP3 ==1 & H4WP5<= ages_w1~1,
                                                                   H4WP16==1 & H4WP18<=ages_w1~1,
                                                                   H4WP30==1 & H4WP32<=ages_w1~1,
                                                                   H4WP9 ==0 | H4WP11>ages_w1~0,
                                                                   H4WP3 ==0 | H4WP5>ages_w1~0,
                                                                   H4WP16==0 | H4WP18>ages_w1~0,
                                                                   H4WP30==0 | H4WP32>ages_w1~0),
                                             parental_warmth1=case_when(H1PF1<6 ~H1PF1),
                                             parental_warmth2=case_when(H1PF4<6 ~H1PF4),
                                             parental_warmth3=case_when(H1PF5<6 ~H1PF5),
                                             parental_warmth4=case_when(H1PF23<6 ~H1PF23),
                                             parental_warmth5=case_when(H1PF24<6~H1PF24),
                                             parental_warmth6=case_when(H1PF25<6~H1PF25),
                                             community_violence1=case_when(H1FV1==0~ 0,
                                                                           H1FV1<6 ~ 1),
                                             community_violence2=case_when(H1FV2==0~ 0,
                                                                           H1FV2<6 ~ 1),
                                             community_violence3=case_when(H1FV3==0~ 0,
                                                                           H1FV3<6 ~ 1),
                                             community_violence4=case_when(H1FV4==0~ 0,
                                                                           H1FV4<6 ~ 1),
                                             ace_community_violence_br=case_when(community_violence1==1  ~ 1,
                                                                                 community_violence2==1  ~ 1,
                                                                                 community_violence3==1  ~ 1,
                                                                                 community_violence4==1  ~ 1,
                                                                                 community_violence1==0  ~ 0,
                                                                                 community_violence2==0  ~ 0,
                                                                                 community_violence3==0  ~ 0,
                                                                                 community_violence4==0  ~ 0),
                                             underres_school1=case_when(PC29A<6 ~ PC29A),
                                             underres_school2=case_when(PC29B<6 ~ PC29B),
                                             underres_school3=case_when(PC29C<6 ~ PC29C),
                                             underres_school_cum=underres_school1+underres_school2+underres_school3,
                                             underres_school_quantile= underres_school_cum %>% ntile(5),
                                             ace_underres_school_br=case_when(underres_school_quantile==1~1,
                                                                              underres_school_quantile>1~0),
                                             neighborhood1=case_when(PA33<6 ~ PA33),
                                             neighborhood2=case_when(PA34<6 ~ PA34),
                                             neighborhood_cum=neighborhood1+neighborhood2,
                                             neighborhood_quantile= neighborhood_cum %>% ntile(3),
                                             ace_neighborhood_br=case_when(neighborhood_quantile<3~0,
                                                                           neighborhood_quantile==3~1),
                                             hhincome=case_when(PA55<9996~PA55),
                                             hhincome1000=hhincome*1000,
                                             size=S27)
  
  waves_ace <- waves_ace %>% mutate(ace_parental_warmth=dplyr::select(., starts_with("parental_warmth"))%>%apply(1, mean, na.rm = TRUE),
                                    ace_parental_warmth_quantile=ace_parental_warmth %>% ntile(5),
                                    ace_parental_warmth_quantile_br=case_when(ace_parental_warmth_quantile==1~1,
                                                                              ace_parental_warmth_quantile>1~0),
                                    ace_poverty_br=case_when(hhincome1000<7360+2480*(size-1) ~ 1,
                                                             hhincome1000>7360+2480*(size-1) ~ 0))
  waves_ace <- waves_ace %>% mutate(ace_br_any=case_when(ace_physical_neglet_br==1 |
                                                           ace_foster_br==1 |
                                                           ace_sexual_br==1 |
                                                           ace_physical_abuse_br==1 |
                                                           ace_alcoholism_br==1 |
                                                           ace_parental_warmth_quantile_br==1 |
                                                           ace_community_violence_br==1 |
                                                           ace_underres_school_br==1 |
                                                           ace_jail_br==1 |
                                                           ace_poverty_br==1 |
                                                           ace_neighborhood_br==1 ~ 1,
                                                         ace_jail_br!=1 |
                                                           ace_poverty_br!=1 |
                                                           ace_physical_neglet_br!=1 |
                                                           ace_foster_br!=1 |
                                                           ace_sexual_br!=1 |
                                                           ace_physical_abuse_br!=1 |
                                                           ace_alcoholism_br!=1 |
                                                           ace_parental_warmth_quantile_br!=1 |
                                                           ace_community_violence_br!=1 |
                                                           ace_underres_school_br!=1 |
                                                           ace_neighborhood_br!=1 ~ 0),
                                    ace_br_cum=dplyr::select(., .vars = c("ace_physical_neglet_br",
                                                                          "ace_physical_abuse_br",
                                                                          "ace_jail_br",
                                                                          "ace_foster_br",
                                                                          "ace_sexual_br",
                                                                          "ace_alcoholism_br",
                                                                          "ace_community_violence_br",
                                                                          "ace_parental_warmth_quantile_br",
                                                                          "ace_underres_school_br",
                                                                          "ace_neighborhood_br",
                                                                          "ace_poverty_br")) %>%
                                      apply(1, function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))))
  
  waves_ace <- waves_ace %>% mutate(ace_br_cum2=case_when(ace_br_cum==0  ~ 0,
                                                          ace_br_cum==1  ~ 1,
                                                          ace_br_cum==2  ~ 2))
  
  
  waves                   = waves %>% left_join(waves_ace, by = "AID") 
  
  ########################################################
  # race discrimnation variables
  ########################################################
  
  waves_race_discrim = waves_full %>% 
    dplyr::select(AID, matches("H5MN5[A-E]"), matches("H5MN6[A-L]"), H5MN7) %>% 
    dplyr::mutate(discrim1 = (H5MN5A - 1) + (H5MN5B -1) + (H5MN5C -1) + (H5MN5D -1) + (H5MN5E -1),
                  discrim2 = ifelse(H5MN6A==1 | H5MN6D==1, discrim1, 0)) %>% 
    dplyr::mutate_at(vars(matches("H5MN5[A-E]")), ~ case_when(. %in% c(3,4) ~ 1,
                                                              . %in% c(1,2) ~ 0)) %>% 
    dplyr::mutate(totdiscrim1 = H5MN5A + H5MN5B + H5MN5C + H5MN5D + H5MN5E + H5MN7,
                  totdiscrim2 = ifelse(H5MN6A==1 | H5MN6D==1, totdiscrim1, 0),
                  countdiscrimwhy = case_when(totdiscrim1 == 0 ~ "no_discrim",
                                              totdiscrim1 == 1 & H5MN7 == 1 ~ "0",
                                              totdiscrim1 > 0 ~ select(., matches("H5MN6[A-L]")) %>%
                                                rowSums(na.rm = TRUE) %>% as.character())) %>% 
    dplyr::select(-starts_with("H")) %>% 
    dplyr::left_join(waves_full %>% 
                       dplyr::select(AID, matches("H5MN5[A-E]"), matches("H5MN6[A-L]"), H5MN7, H3IR17, H4DA17)) %>% 
    mutate_at(vars(H4DA17), ~ ifelse(. %in% c(995, 996, 997, 998), NA, H4DA17))
  
  waves                   = waves %>% left_join(waves_race_discrim, by = "AID") 
  
  ########################################################
  # full siblings family ID
  ########################################################
  
  
  # keep only the AID who has full siblings, convert the rest to NA, AID is the same as SIB_AID1 for some records, survey mistake 
  fullsibling = sibling_w3 %>% 
    mutate_at(vars(c(SIB_AID1, SIB_AID2, SIB_AID3, SIB_AID4)), list(~ .x %>% as.character())) %>%
    mutate(SIB_AID1 = ifelse(SIB_REL1 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID1 != AID, SIB_AID1,NA), NA))%>%
    mutate(SIB_AID2 = ifelse(SIB_REL2 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID2 != AID, SIB_AID2,NA), NA))%>%
    mutate(SIB_AID3 = ifelse(SIB_REL3 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID3 != AID, SIB_AID3,NA), NA))%>%
    mutate(SIB_AID4 = ifelse(SIB_REL4 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID4 != AID, SIB_AID4,NA), NA)) %>% 
    filter(!((is.na(SIB_AID1) & is.na(SIB_AID2) & is.na(SIB_AID3) & is.na(SIB_AID4)))) %>% 
    mutate(famid_fullsib = str_c("fam_",c(1:dim(.)[1])))
  
  # reconstruct a dataframe for all the AID who has full siblings as one column, and their famid as another column
  # remove some duplicated rows
  temp = list(t0 = fullsibling %>% select(AID, famid_fullsib),
              t1 = fullsibling %>% select(SIB_AID1, famid_fullsib) %>% dplyr::rename(AID = SIB_AID1),
              t2 = fullsibling %>% select(SIB_AID2, famid_fullsib) %>% dplyr::rename(AID = SIB_AID2),
              t3 = fullsibling %>% select(SIB_AID3, famid_fullsib) %>% dplyr::rename(AID = SIB_AID3),
              t4 = fullsibling %>% select(SIB_AID4, famid_fullsib) %>% dplyr::rename(AID = SIB_AID4)) %>% 
    purrr::reduce(bind_rows) %>% 
    filter(!is.na(AID)) %>% 
    nest(AIDs = AID) %>%
    mutate(AIDs = AIDs %>% unlist(recursive = F))  %>% 
    mutate(AIDs = AIDs %>% map(str_replace_na) %>% map(as.numeric) %>% map(sort)) %>%
    column_to_rownames(var = "famid_fullsib") %>%
    unique() %>%
    tidyr::separate(col = AIDs, sep = ",", into = c("AID1", "AID2", "AID3")) %>%
    rownames_to_column("famid_fullsib") %>% 
    mutate_at(vars(c("AID1", "AID2", "AID3")), ~ gsub("[^0-9.-]", "",.))
  
  # full siblings AIDs and their family id
  waves_fullsibling = purrr::reduce(list(t1 = temp %>% select(AID1, famid_fullsib) %>% dplyr::rename(AID = AID1),
                                         t2 = temp %>% select(AID2, famid_fullsib) %>% dplyr::rename(AID = AID2),
                                         t3 = temp %>% select(AID3, famid_fullsib) %>% dplyr::rename(AID = AID3)),
                                    bind_rows) %>%
    filter(!is.na(AID)) %>%
    filter(!(famid_fullsib %in% c("fam_585", "fam_658")))
  # test duplicate cases
  # bbb=waves_fullsibling %>% duplicated()
  # bbbbb=bbb %>% as.data.frame()
  # temp = waves_fullsibling %>% cbind(bbbbb)
  
  # remove manuelly fam_658 fam_585	and keep fam_2264
  # 97572230	91712236 91582938 
  # 1	91582938	NA	sister	97572230	twin brother	NA		NA		1	fam_585
  # 2	91712236	NA	brother	97572230	brother	NA		NA		1	fam_658
  # 3	97572230	91582938	sister	91712236	twin brother	NA		NA		1	fam_2264s
  
  waves                   = waves %>% left_join(waves_fullsibling, by = "AID") 
  
  ########################################################
  # SAVE
  ########################################################
  waves %>% saveRDS(file = str_c(data_output, "/", output_fname))
}
