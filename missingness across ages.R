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


#--- Load in unimputed dataset ---#
load("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/da_updated.Rdata")

#-- Create age-residualized Tanner & PDS variables
#-- Note: BMI_7y missing 207/1138 & tannerstg_12y missing 60/1138 - resulting in 889 obs
summary(lm(tannerstg_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da))
da$tanner12_age_res <- da$tannerstg_12y-(-2.8685+0.001*da$c_age_days_comp_d_qu_12y+0.0809*da$BMI_7y) 

#-- Note: BMI_7y missing 207/1138 & puberty_12y missing 2/1138 - resulting in 930 obs
summary(lm(puberty_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da))
da$pds12_age_res <- da$puberty_12y-(-3.255+0.000991*da$c_age_days_comp_d_qu_12y+0.0607*da$BMI_7y) 

names(da)

#-- Rename the age variables
names(da)[71:78] <- c("age7","age9","age10","age11","age12","age14","age15","age16")
names(da)[34:38] <- c("menarch_age9","menarch_age10","menarch_age11","menarch_age12","menarch_age14")
names(da)


#-- Identify age at menarche among girls:
da <- da %>%
  group_by(aid) %>%
  mutate(
    age_menarche = first(na.omit(c(menarch_age9,menarch_age10,menarch_age11,menarch_age12,menarch_age14)))
  ) %>% ungroup()


#-- Split the data by sex
da_b <- da %>% filter(female_d==0) 
da_g <- da %>% filter(female_d==1) 

#-- Create a BMI-residualized APHV variable
summary(lm(age_peak_velocity_years~BMI_7y,data=da_b))
da_b$age_peak_velocity_years_bmi <- da_b$age_peak_velocity_years-(14.80-0.10160*da_b$BMI_7y) 
summary(lm(age_peak_velocity_years~BMI_7y,data=da_g))
da_g$age_peak_velocity_years_bmi <- da_g$age_peak_velocity_years-(13.46-0.12973*da_g$BMI_7y) 

da <- da %>%
  mutate(
    age_peak_velocity_years_bmi = case_when(
      female_d==0 ~ age_peak_velocity_years-(14.80-0.10160*BMI_7y),
      female_d==1 ~ age_peak_velocity_years-(13.46-0.12973*BMI_7y)
    )
  ) 

summary(da_g$age_menarche)


ages <- da %>% 
  select(aid,age11,age12,age14,age15,age16)

ages_missing <- missing_pattern(ages)
