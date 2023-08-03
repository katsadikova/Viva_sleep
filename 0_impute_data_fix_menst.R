#~~ 0_impute_data_20x.R
#~~ Date: 4/27/2023
#----------------------------------------

library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)
library(gtools)
library(tableone)
library(corrplot)
library(reshape2)
library(mice)
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)

#--- Data set sent by Sheryl - contains 2128 rows, 177 variables ---#
d<-read_sas("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/new pull - 3-30-2023/ks_093022.sas7bdat") 
names(d)
hist(d$avgsleeptime_wd/60)
summary(d$sleep_schoolday_cq14)
summary(d$sleep_weekday_child_cq15)
summary(d$sleep_weekday_cq16)
summary(d$teen_sleep_weekday_tq17)


d %>% summarise(total_non_na = sum(!is.na(c_age_days_comp_d_qu_11y) & !is.na(female_d))) #704
d %>% summarise(total_non_na = sum(!is.na(c_age_days_comp_d_qu_12y) & !is.na(female_d))) #1138
d %>% summarise(total_non_na = sum(!is.na(age_days_cq14))) #409
d %>% summarise(total_non_na = sum(!is.na(age_days_cq15))) #467
d %>% summarise(total_non_na = sum(!is.na(age_days_cq16))) #726

d_actig <- d %>% filter(!is.na(avgsleeptime_wd))
d_actig %>% summarise(total_non_na = sum(!is.na(c_age_days_comp_d_qu_12y))) #871
d_actig %>% summarise(total_non_na = sum(!is.na(age_days_cq14))) #320
d_actig %>% summarise(total_non_na = sum(!is.na(age_days_cq15))) #388
d_actig %>% summarise(total_non_na = sum(!is.na(age_days_cq16))) #572


d_pre_imp <- d %>%
  mutate(
    csection = case_when(csection_d=="Y"~1,
                         csection_d=="N"~0,
                         T~NA_real_),
    race_amind = ifelse(race2_mom_epi_epia_d=="amind",1,0),
    race_asian = ifelse(race2_mom_epi_epia_d=="asian",1,0),
    race_black = ifelse(race2_mom_epi_epia_d=="black",1,0),
    race_hisp = ifelse(race2_mom_epi_epia_d=="hispa",1,0),
    race_morethan1 = ifelse(race2_mom_epi_epia_d=="more than 1 race",1,0),
    race_other = ifelse(race2_mom_epi_epia_d=="other",1,0),
    ADULT_QU7Y = case_when(ADULT_QU7Y==-9 ~ NA_real_,
                           ADULT_QU7Y==-2 ~ NA_real_,
                           T~ADULT_QU7Y),
    CHILD_QU7Y = case_when(CHILD_QU7Y==-9 ~ NA_real_,
                           CHILD_QU7Y==-2 ~ NA_real_,
                           T~CHILD_QU7Y),
    HHsize_QU7Y = ADULT_QU7Y+CHILD_QU7Y,
    MARITAL_QU7Y = case_when(MARITAL_QU7Y==-9 ~ NA_real_,
                             T~MARITAL_QU7Y),
    HINCOME_QU7Y = case_when(HINCOME_QU7Y==-9 ~ NA_real_,
                             T~HINCOME_QU7Y),
    WELFARE_QU7Y = case_when(WELFARE_QU7Y==-9 ~ NA_real_,
                             WELFARE_QU7Y==2 ~ 0,
                             T~WELFARE_QU7Y),
    #--Impute tanner at 12 using nearest available values
    tannerstg_12y = case_when(is.na(tannerstg_12y) & !is.na(tannerstg_11y) & !is.na(tannerstg_14y) ~ round((tannerstg_11y+tannerstg_14y)/2,0),
                              T ~ tannerstg_12y),
    #--Impute PDS at 12 using nearest available values
    puberty_12y = case_when(is.na(puberty_12y) & !is.na(puberty_11y) & !is.na(puberty_14y) ~ round((puberty_11y+puberty_14y)/2,0),
                            T ~ puberty_12y),
    sex = case_when(female_d==1~"Female",
                    female_d==0~"Male"),
    BMI_7y = weight_7y/((height_7y/100)*(height_7y/100)),
    C_PERCENT_CA7Y = case_when(C_PERCENT_CA7Y==-2~NA_real_,
                               T~C_PERCENT_CA7Y),
    tv_schoolday_12y = case_when(tv_schoolday_12y==-9~NA_real_,
                                 T~tv_schoolday_12y),
    tv_weekend_12y = case_when(tv_weekend_12y==-9~NA_real_,
                               T~tv_weekend_12y),
    mother_weight_12y_bin = case_when(mother_weight_12y==-9 ~ NA_real_,
                                      mother_weight_12y>1 ~ 1,
                                      mother_weight_12y==1 ~ 0),
    c_age_days_men_date_10y = case_when(!is.na(c_age_days_men_date_9y) ~ c_age_days_men_date_9y,
                                        T ~ c_age_days_men_date_10y),
    c_age_days_men_date_11y = case_when(!is.na(c_age_days_men_date_10y) ~ c_age_days_men_date_10y,
                                        T ~ c_age_days_men_date_11y),
    c_age_days_men_date_12y = case_when(!is.na(c_age_days_men_date_11y) ~ c_age_days_men_date_11y,
                                        T ~ c_age_days_men_date_12y),
    men_date_12y = case_when(is.na(c_age_days_men_date_12y) & !is.na(child_age_days_period_qu14) ~ child_age_days_period_qu14,
                             T ~ c_age_days_men_date_12y)
  ) %>%
  #-- Only keep rows with data at ET visit
  filter(!is.na(c_age_days_comp_d_qu_12y))%>%
  select(familyid,aid,sex,age_mom_enroll_d,
         race_amind,race_asian,race_black,
         race_hisp,race_morethan1,race_other,
         ses_index_new_03,
         BMI_7y,C_PERCENT_CA7Y,
         puberty_7y,puberty_8y,puberty_9y,
         c_age_days_COMP_D_QU7Y,c_age_days_comp_d_qu_9y,c_age_days_comp_d_qu_10y,
         tannerstg_10y,puberty_10y,
         c_age_days_comp_d_qu_11y,
         tannerstg_11y,puberty_11y,feeling_score_11y,
         age_peak_velocity_years,height_peak_height_velocity_cm,
         c_age_days_comp_d_qu_12y,
         #tv_schoolday_12y,tv_weekend_12y,
         mother_weight_12y_bin,
         tannerstg_12y,men_date_12y, puberty_12y,avgsleeptime,avgsleeptime_wd,avgsleeptime_we,avgwaso,avgwaso_wd,avgwaso_we,
         avgefficiency,avgefficiency_wd,avgefficiency_we,agey_actigraph,feeling_score_ET,
         age_days_cq14,
         tannerstg_14y,puberty_14y,sleep_schoolday_cq14,feeling_score_14y,
         age_days_cq15,
         tannerstg_15y,sleep_weekend_child_cq15,sleep_weekday_child_cq15,sleep_weekend_child_cq15,feeling_score_15y,
         age_days_cq16,
         tannerstg_16y,sleep_weekday_cq16,sleep_weekend_cq16,feeling_score_16y,
         teen_sleep_weekday_tq17,teen_sleep_weekend_tq17,C5_depr_tscore)
names(d_pre_imp)
summary(d_pre_imp)
save(d_pre_imp, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_pre_imp.Rdata")


#-------------------------------------------------------------------#
#-- Imputation code source: https://rpubs.com/kaz_yos/mice-exclude

## Proportion missing
prom_miss <- gather(d_pre_imp, key="aid") %>%
  group_by(aid) %>%
  summarize(prop_na = mean(is.na(value))) 

## Configure parallelization
## Detect core count
nCores <- min(parallel::detectCores(), 8)
## Used by parallel::mclapply() as default
options(mc.cores = nCores)
## Used by doParallel as default
options(cores = nCores)
## Register doParallel as the parallel backend with foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Extract all variable names in dataset
allVars <- names(d_pre_imp)
allVars

## names of variables with missingness
missVars <- names(d_pre_imp)[colSums(is.na(d_pre_imp)) > 0]
missVars

predictorMatrix <- matrix(0, ncol = length(allVars), nrow = length(allVars))
rownames(predictorMatrix) <- allVars


colnames(predictorMatrix) <- allVars

names(d_pre_imp)
cat("
###  Specify Variables informing imputation\n")
## These can be either complete variables or variables with missingness.
## Those with missingness must be imputed - need to exclude men_date_12y
## Explicitly specify.
imputerVars <- names(d_pre_imp)[c(3:30,32:61)]
imputerVars 
## Keep variables that actually exist in dataset
imputerVars <- intersect(unique(imputerVars), allVars)
imputerVars
imputerMatrix <- predictorMatrix
imputerMatrix[,imputerVars] <- 1
imputerMatrix


cat("
###  Specify variables with missingness to be imputed \n")
## Could specify additional variables that are imputed,
## but does not inform imputation.
imputedOnlyVars <- c("men_date_12y")
## Imputers that have missingness must be imputed.
imputedVars <- intersect(unique(c(imputedOnlyVars, imputerVars)), missVars)
imputedVars
imputedMatrix <- predictorMatrix
imputedMatrix[imputedVars,] <- 1
imputedMatrix

cat("
###  Construct a full predictor matrix (rows: imputed variables; cols: imputer variables)\n")
## Keep correct imputer-imputed pairs only
predictorMatrix <- imputerMatrix * imputedMatrix
## Diagonals must be zeros (a variable cannot impute itself)
diag(predictorMatrix) <- 0
predictorMatrix


cat("
###  Dry-run mice for imputation methods\n")
dryMice <- mice(data = d_pre_imp, m = 1, predictorMatrix = predictorMatrix, maxit = 0)
## Update predictor matrix
predictorMatrix <- dryMice$predictorMatrix
cat("###   Imputers (non-zero columns of predictorMatrix)\n")
imputerVars <- colnames(predictorMatrix)[colSums(predictorMatrix) > 0]
imputerVars
cat("###   Imputed (non-zero rows of predictorMatrix)\n")
imputedVars <- rownames(predictorMatrix)[rowSums(predictorMatrix) > 0]
imputedVars
cat("###   Imputers that are complete\n")
setdiff(imputerVars, imputedVars)
cat("###   Imputers with missingness\n")
intersect(imputerVars, imputedVars)
cat("###   Imputed-only variables without being imputers\n")
setdiff(imputedVars, imputerVars)
cat("###   Variables with missingness that are not imputed\n")
setdiff(missVars, imputedVars)
cat("###   Relevant part of predictorMatrix\n")
predictorMatrix[rowSums(predictorMatrix) > 0, colSums(predictorMatrix) > 0]


## Empty imputation method to really exclude variables
## http://www.stefvanbuuren.nl/publications/MICE%20in%20R%20-%20Draft.pdf
##
## MICE will automatically skip imputation of variables that are complete.
## One of the problems in previous versions of MICE was that all incomplete
## data needed to be imputed. In MICE 2.0 it is possible to skip imputation
## of selected incomplete variables by specifying the empty method "".
## This works as long as the incomplete variable that is skipped is not being
## used as a predictor for imputing other variables.
## Note: puttting zeros in the predictorMatrix alone is NOT enough!
##
dryMice$method[setdiff(allVars, imputedVars)] <- ""
cat("###   Methods used for imputation\n")
dryMice$method[sapply(dryMice$method, nchar) > 0]


cat("
###  Run mice\n")
M <- 20
cat("### Imputing", M, "times\n")

## Set seed for reproducibility
#set.seed(3561126)
set.seed(317)
## Parallelized execution
miceout <- foreach(i = seq_len(M), .combine = ibind) %dorng% {
  cat("### Started iteration", i, "\n")
  miceout <- mice(data = d_pre_imp, m = 1, print = TRUE,
                  predictorMatrix = predictorMatrix, method = dryMice$method,
                  MaxNWts = 2000)
  cat("### Completed iteration", i, "\n")
  ## Make sure to return the output
  miceout
}


cat("
###  Show mice results\n")
## mice object ifself
miceout
## Variables that no longer have missingness after imputation
cat("###   Variables actually imputed\n")
actuallyImputedVars <-
  setdiff(names(d_pre_imp)[colSums(is.na(d_pre_imp)) > 0],
          names(complete(miceout, action = 1))[colSums(is.na(complete(miceout, action = 1))) > 0])
actuallyImputedVars

## Examine discrepancies
cat("###   Variables that were unexpectedly imputed\n")
setdiff(actuallyImputedVars, imputedVars)
cat("###   Variables that were planned for MI but not imputed\n")
setdiff(imputedVars, actuallyImputedVars)

## Still missing variables
cat("###   Variables still having missing values\n")
names(complete(miceout, action = 1))[colSums(is.na(complete(miceout, action = 1))) > 0]

#-- Save Imputed Data!
d_imp <- miceout
save(d_imp, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_imp_nofactors_fixed_chron.Rdata")

names(complete(d_imp))

hist(d_pre_imp$avgsleeptime_wd/60)
summary(as.factor(d_pre_imp$avgsleeptime_wd/60))
hist(d_pre_imp$sleep_schoolday_cq14)
summary(as.factor(d_pre_imp$sleep_schoolday_cq14))
hist(d_pre_imp$sleep_weekday_child_cq15)
summary(as.factor(d_pre_imp$sleep_weekday_child_cq15))


