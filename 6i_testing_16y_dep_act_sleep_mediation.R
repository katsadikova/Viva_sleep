#--------------------------------------------------------------------------#
#--- 2i_repeated_cross_sectional_models.R
#--- Date: 4/30/2023
#--- Description: Run linear models with repeated cross-sectional outcomes
#--------------------------------------------------------------------------#

library(tidyverse)
library(quadprog)
library(haven)
library(ggplot2)
library(tableone)
library(kableExtra)
library(raters)
library(gmodels)
library(ggpubr)
library(stargazer)
library(gtsummary)
library(tableone)
library(nlme)
library(nlstools)
library(CMAverse)
library(mice)

#-- Load in long-form data (out)
load(file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_imp_long.Rdata")
#-- Load in original mice object (d_imp)
load(file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_imp_nofactors_fixed_chron.Rdata")
#-- Load pre-imputation data
load(file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_pre_imp.Rdata")

summary(d_pre_imp$feeling_score_16y)
nrow(d_pre_imp)


#----------------------------------------------------------------------------------------#
#-- Make edits to the mids object d_imp

#--- convert to long-from data set to add the "max_problems" variable
long <- complete(d_imp, action="long", include=TRUE)
long$ses_index_new_10pct <- ifelse(long$ses_index_new_03<=7,1,0)
long$miss_actigraphy <-ifelse(is.na(d_pre_imp$avgsleeptime),1,0)
long$miss_ET <-ifelse(is.na(d_pre_imp$c_age_days_comp_d_qu_12y),1,0)
long$avgsleeptime_wd_hr <- long$avgsleeptime_wd/60
long$tannerstg_12y <- as.numeric(long$tannerstg_12y)
long$overweight <- ifelse(long$BMI_7y>=25,1,0)
long$obese <- ifelse(long$BMI_7y>=30,1,0)
long$men_date_12y <- long$men_date_12y/365.25
long$act_sl_8hr <- ifelse(long$avgsleeptime_wd_hr>=8,1,0)
long$sr_sl14_8hr <- ifelse(long$sleep_schoolday_cq14>=8,1,0)
long$sr_sl15_8hr <- ifelse(long$sleep_weekday_child_cq15>=8,1,0)
long$sr_sl16_8hr <- ifelse(long$sleep_weekday_cq16>=8,1,0)
long$promis_modsev <- ifelse(long$C5_depr_tscore>=65,1,0)

#--- convert back to mids object
d_imp <- as.mids(long)
#----------------------------------------------------------------------------------------#

#-- Subset the imputed mids object by sex
d_imp_b <- filter(d_imp, sex == "Male")
nrow(complete(d_imp_b)) #585
d_imp_g <- filter(d_imp, sex == "Female")
nrow(complete(d_imp_g)) #553


b_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+as.factor(sex) + ses_index_new_10pct + c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y") 
b_dep16adj_pds

b_sleep16adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+as.factor(sex) + ses_index_new_10pct + c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y") 
b_sleep16adj_pds

b_depsleep16adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_16y ~ avgsleeptime_wd_hr + puberty_12y+BMI_7y+as.factor(sex) + ses_index_new_10pct + c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y") 
b_depsleep16adj_pds
