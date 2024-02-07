#----------------------------------------------------------------------------#
#--- 1i_describe_data.R
#--- Date: 6/20/2023
#--- Dat updated: 2/7/2024
#--- Description: Create descriptive tables looking at distributions of 
#---              pubertal timing measures, covariates, and outcomes
#----------------------------------------------------------------------------#

library(tidyverse)
library(quadprog)
library(haven)
library(ggplot2)
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
library(lme4)
library(broom.mixed)


#-- Load in original mice object (d_imp)
load(file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_imp_nofactors_fixed_chron.Rdata")
#-- Load pre-imputation data
load(file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/d_pre_imp.Rdata")
#-- Load in new data pull from Sheryl (from 7/5/2023) - with age at the mid-teen visit
age17 <- read_sas("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/new pull - 07-05-2023/ks_093022.sas7bdat") %>%
  dplyr::select(aid,age_years_comp_d_tq17)


#-- Observed Ns across time points
d_pre_imp %>% filter(!is.na(feeling_score_ET)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(feeling_score_14y)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(feeling_score_15y)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(feeling_score_16y)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(C5_depr_tscore)) %>% group_by(sex) %>% summarize(n=n())

d_pre_imp %>% filter(!is.na(avgsleeptime_wd)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(sleep_schoolday_cq14)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(sleep_weekday_child_cq15)) %>% group_by(sex) %>% summarize(n=n())
d_pre_imp %>% filter(!is.na(sleep_weekday_cq16)) %>% group_by(sex) %>% summarize(n=n())


d <- complete(d_imp) %>%
  merge(., y=age17, by="aid") %>%
  mutate(
    men_date_12y = men_date_12y/365.25,
    avgsleeptime_wd = avgsleeptime_wd/60,
    age7 = c_age_days_COMP_D_QU7Y/365.25,
    age12 = c_age_days_comp_d_qu_12y/365.25,
    age14 = age_days_cq14/365.25,
    age15 = age_days_cq15/365.25,
    age16 = age_days_cq16/365.25,
    age17 = age_years_comp_d_tq17
  )

hist(d$puberty_7y)
hist(d$puberty_10y)
hist(d$puberty_11y)
hist(d$puberty_12y)
hist(d$puberty_14y)

summary(as.factor(d$puberty_12y))

## Create table comparing distributions of pubertal timing measures and covariates by sex

## Make categorical variables factors
varsToFactor <- c("tannerstg_12y")
d[varsToFactor] <- lapply(d[varsToFactor], factor)

## Create a variable list
vars <- c("ses_index_new_03","BMI_7y","age7","age12","age14","age15","age16","age17",
          "age_peak_velocity_years","puberty_12y","tannerstg_12y",
          "avgsleeptime_wd","sleep_schoolday_cq14","sleep_weekday_child_cq15","sleep_weekday_cq16",
          "C5_depr_tscore")

## Create Table 1 stratified by delay
table_by_sex <- CreateTableOne(vars = vars,
                               strata = c("sex"), 
                               includeNA = F, 
                               addOverall = T,
                               data = d)
table_by_sex <- data.frame(print(table_by_sex, missing=F, catDigits = 1, contDigits = 1))
table_by_sex <- tibble::rownames_to_column(table_by_sex, "Characteristic")

# Column labels
names(table_by_sex) <- c("Characteristic", "Overall", "Girls", "Boys", "p-value","test")

dim(table_by_sex)

# Row labels
length(table_by_sex$Characteristic)

table_by_sex <- data.frame(table_by_sex) %>%
  mutate(
    rownames = case_when(
      Characteristic=="ses_index_new_03 (mean (SD))" ~ "SES Index (mean (SD))",
      Characteristic=="BMI_7y (mean (SD))" ~ "Mid-childhood BMI (mean (SD))",
      Characteristic=="age7 (mean (SD))" ~ "Mid-childhood (mean (SD))",
      Characteristic=="age12 (mean (SD))" ~ "Early adolescence (mean (SD))",
      Characteristic=="age14 (mean (SD))" ~ "14-year (mean (SD))",
      Characteristic=="age15 (mean (SD))" ~ "15-year (mean (SD))",
      Characteristic=="age16 (mean (SD))" ~ "16-year (mean (SD))",
      Characteristic=="age17 (mean (SD))" ~ "17-year (mean (SD))",
      Characteristic=="age_peak_velocity_years (mean (SD))" ~ "APHV (mean (SD))",
      Characteristic=="puberty_12y (mean (SD))" ~ "PDS (mean (SD))",
      Characteristic=="tannerstg_12y (%)" ~ "Tanner pubic hair stage (n (%))",
      Characteristic=="avgsleeptime_wd (mean (SD))" ~ "Early adolescence visit, actigraphy, hours (mean (SD))",
      Characteristic=="sleep_schoolday_cq14 (mean (SD))" ~ "14-year visit, self-report, hours (mean (SD))",
      Characteristic=="sleep_weekday_child_cq15 (mean (SD))" ~ "15-year visit, self-report, hours (mean (SD))",
      Characteristic=="sleep_weekday_cq16 (mean (SD))" ~ "16-year visit, self-report, hours (mean (SD))",
      Characteristic=="C5_depr_tscore (mean (SD))" ~ "Mid-adolescence PROMIS t-score (mean (SD))",
      T ~ Characteristic
    )
  ) %>%
  dplyr::select(rownames,Girls,Boys,Overall)

table_by_sex

write.csv(table_by_sex, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/table_1.csv")


check_g <- d %>%
  filter(sex=="Female")
summary(check_g$age17)
check_g <- d %>%
  filter(sex=="Male")
summary(check_g$age17)


cor(d$feeling_score_ET, d$C5_depr_tscore)

summary(lm(data=d, C5_depr_tscore~feeling_score_11y))
