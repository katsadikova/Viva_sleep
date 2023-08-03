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


#----------------------------------------------------------------------------------------#
#-- Crude relationships
#-- Tanner & Dep
g_dep11_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="11Y") 
g_dep12_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y") 
g_dep14_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y") 
g_dep15_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y") 
g_dep16_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y") 
g_dep17_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="17Y") 

b_dep11_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="11Y") 
b_dep12_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y") 
b_dep14_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y") 
b_dep15_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y") 
b_dep16_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y") 
b_dep17_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="17Y")  


#-- Significant interaction between tanner & sex at the ET visit - run repeated cross-sectionals separately by sex
dep12_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_ET ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y")
dep12_tanner

#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_crude <- rbind(g_dep11_tanner,g_dep12_tanner,g_dep14_tanner,g_dep15_tanner,g_dep16_tanner,g_dep17_tanner,
                      b_dep11_tanner,b_dep12_tanner,b_dep14_tanner,b_dep15_tanner,b_dep16_tanner,b_dep17_tanner) %>%
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep11_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="11Y") 
g_dep12_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y") 
g_dep14_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y")
g_dep15_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y") 
g_dep16_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y")
g_dep17_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="17Y")

b_dep11_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="11Y") 
b_dep12_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y") 
b_dep14_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y") 
b_dep15_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y") 
b_dep16_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y") 
b_dep17_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="17Y")

#-- Significant interaction between PDS & sex at the ET visit - run repeated cross-sectionals separately by sex
dep12_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_ET ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y")
dep12_pds

#-- Put all the coefficients for tanner together, format sub-table
pds_dep_crude <- rbind(g_dep11_pds,g_dep12_pds,g_dep14_pds,g_dep15_pds,g_dep16_pds,g_dep17_pds,
                   b_dep11_pds,b_dep12_pds,b_dep14_pds,b_dep15_pds,b_dep16_pds,b_dep17_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- APHV & Dep (adjusting for age at 12Y visit doesn't change associations for APHV)
g_dep11_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="11Y")
g_dep12_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y") 
g_dep14_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y") 
g_dep15_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y")
g_dep16_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y") 
g_dep17_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="17Y")

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="11Y") 
b_dep12_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y")
b_dep14_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y") 
b_dep15_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y") 
b_dep16_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y") 
b_dep17_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="17Y")

#-- No significant interaction between APHV & sex at the ET visit - run repeated cross-sectionals separately by sex despite, since differences for PDS & Tanner
dep12_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_ET ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="12Y")
dep12_aphv

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_crude <- rbind(g_dep11_aphv,g_dep12_aphv,g_dep14_aphv,g_dep15_aphv,g_dep16_aphv,g_dep17_aphv,
                    b_dep11_aphv,b_dep12_aphv,b_dep14_aphv,b_dep15_aphv,b_dep16_aphv,b_dep17_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- Menarche & Dep
g_dep11_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="11Y") 
g_dep12_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y") 
g_dep14_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y")
g_dep15_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y") 
g_dep16_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y")
g_dep17_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="17Y")

#-- Put all the coefficients for tanner together, format sub-table
m_dep_crude <- rbind(g_dep11_m,g_dep12_m,g_dep14_m,g_dep15_m,g_dep16_m,g_dep17_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_crude <- rbind(aphv_dep_crude, tanner_dep_crude, pds_dep_crude, m_dep_crude) 


#----------------------------------------------------------------------------------------#
#-- BMI-adjusted relationships
#-- Tanner & Dep
g_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="11Y") 
g_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y") 
g_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y") 
g_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y") 
g_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y")  
g_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="17Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="11Y")
b_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y") 
b_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y") 
b_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y") 
b_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y") 
b_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="17Y") 


#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_adj <- rbind(g_dep11adj_tanner,g_dep12adj_tanner,g_dep14adj_tanner,g_dep15adj_tanner,g_dep16adj_tanner,g_dep17adj_tanner,
                    b_dep11adj_tanner,b_dep12adj_tanner,b_dep14adj_tanner,b_dep15adj_tanner,b_dep16adj_tanner,b_dep17adj_tanner) %>% 
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="11Y") 
g_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y") 
g_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y")
g_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y") 
g_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y")
g_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="17Y")

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="11Y") 
b_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y") 
b_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y") 
b_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y") 
b_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y") 
b_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="17Y") 


#-- Put all the coefficients for tanner together, format sub-table
pds_dep_adj <- rbind(g_dep11adj_pds,g_dep12adj_pds,g_dep14adj_pds,g_dep15adj_pds,g_dep16adj_pds,g_dep17adj_pds,
                 b_dep11adj_pds,b_dep12adj_pds,b_dep14adj_pds,b_dep15adj_pds,b_dep16adj_pds,b_dep17adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Dep
g_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="11Y")
g_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y") 
g_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y") 
g_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y")
g_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y") 
g_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="17Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="11Y") 
b_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y")
b_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y") 
b_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y") 
b_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y") 
b_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="17Y") 

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_adj <- rbind(g_dep11adj_aphv,g_dep12adj_aphv,g_dep14adj_aphv,g_dep15adj_aphv,g_dep16adj_aphv,g_dep17adj_aphv,
                  b_dep11adj_aphv,b_dep12adj_aphv,b_dep14adj_aphv,b_dep15adj_aphv,b_dep16adj_aphv,b_dep17adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- Menarche & Dep
g_dep11adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="11Y") 
g_dep12adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y") 
g_dep14adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y")
g_dep15adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y") 
g_dep16adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y")
g_dep17adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="17Y")

#-- Put all the coefficients for tanner together, format sub-table
m_dep_adj <- rbind(g_dep11adj_m,g_dep12adj_m,g_dep14adj_m,g_dep15adj_m,g_dep16adj_m,g_dep17adj_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_adj <- rbind(aphv_dep_adj, tanner_dep_adj, pds_dep_adj, m_dep_adj)

#----------------------------------------------------------------------------------------#
#-- BMI and SES adjusted relationships
#-- Tanner & Dep
g_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+race_black+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="11Y") 
g_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y") 
g_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y") 
g_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y") 
g_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y")  
g_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="17Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="11Y")
b_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y") 
b_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y") 
b_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y") 
b_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y") 
b_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="17Y") 

#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_adj <- rbind(g_dep11adj_tanner,g_dep12adj_tanner,g_dep14adj_tanner,g_dep15adj_tanner,g_dep16adj_tanner,g_dep17adj_tanner,
                        b_dep11adj_tanner,b_dep12adj_tanner,b_dep14adj_tanner,b_dep15adj_tanner,b_dep16adj_tanner,b_dep17adj_tanner) %>% 
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="11Y") 
g_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y") 
g_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y")
g_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y") 
g_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y")
g_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="17Y")

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="11Y") 
b_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y") 
b_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y") 
b_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y") 
b_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y") 
b_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="17Y") 


#-- Put all the coefficients for tanner together, format sub-table
pds_dep_adj <- rbind(g_dep11adj_pds,g_dep12adj_pds,g_dep14adj_pds,g_dep15adj_pds,g_dep16adj_pds,g_dep17adj_pds,
                     b_dep11adj_pds,b_dep12adj_pds,b_dep14adj_pds,b_dep15adj_pds,b_dep16adj_pds,b_dep17adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Dep
g_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="11Y")
g_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y") 
g_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y") 
g_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y")
g_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y") 
g_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="17Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="11Y") 
b_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y")
b_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y") 
b_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y") 
b_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y") 
b_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="17Y") 

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_adj <- rbind(g_dep11adj_aphv,g_dep12adj_aphv,g_dep14adj_aphv,g_dep15adj_aphv,g_dep16adj_aphv,g_dep17adj_aphv,
                      b_dep11adj_aphv,b_dep12adj_aphv,b_dep14adj_aphv,b_dep15adj_aphv,b_dep16adj_aphv,b_dep17adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- Menarche & Dep
g_dep11adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="11Y") 
g_dep12adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y") 
g_dep14adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y")
g_dep15adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y") 
g_dep16adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y")
g_dep17adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="17Y")
g_dep17adj_m

#-- Put all the coefficients for tanner together, format sub-table
m_dep_adj <- rbind(g_dep11adj_m,g_dep12adj_m,g_dep14adj_m,g_dep15adj_m,g_dep16adj_m,g_dep17adj_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_adj_bmi_ses <- rbind(aphv_dep_adj, tanner_dep_adj, pds_dep_adj, m_dep_adj)

#-----------------------------------------------------------
#-- Rerun full-sample analysis with PROMIS depression scores

#-- Test interactions with sex first - none significant.
dep17adj_ai <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years*as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="17Y")
dep17adj_ai
dep17adj_pi <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y")
dep17adj_pi
dep17adj_ti <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="17Y")
dep17adj_ti


#-- Rerun without interaction terms (adjusted for sex) - significant associations with puberty_12y and tannerstg_12y
dep17_a <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="17Y", model="Crude")
dep17_a
dep17_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y", model="Crude")
dep17_p
dep17_t <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="17Y", model="Crude")
dep17_t

dep17adj_a <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="17Y", model="BMI-adj")
dep17adj_a
dep17adj_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y", model="BMI-adj")
dep17adj_p
dep17adj_t <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="17Y", model="BMI-adj")
dep17adj_t

dep17adjses_a <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+as.factor(sex)+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="All", x="APHV", y="17Y", model="BMI- & SES-adj")
dep17adjses_a
dep17adjses_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y", model="BMI- & SES-adj")
dep17adjses_p
dep17adjses_t <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="17Y", model="BMI- & SES-adj")
dep17adjses_t

dep17_adj <- rbind(dep17_a,dep17_p,dep17_t,
                   dep17adj_a,dep17adj_p,dep17adj_t,
                   dep17adjses_a,dep17adjses_p,dep17adjses_t) %>%
  filter(term %in% c("age_peak_velocity_years","puberty_12y","tannerstg_12y")) %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,model,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")
dep17_adj

#----------------------------------------------------------------------------------------#
#-- Look at relationships with sleep duration
#----------------------------------------------------------------------------------------#
#-- Crude relationships
#-- Tanner & Sleep

g_sl12_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y (actigraphy)") 
g_sl14_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y (self-report)") 
g_sl15_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y (self-report)") 
g_sl16_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y (self-report)") 

b_sl12_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y (actigraphy)") 
b_sl14_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y (self-report)") 
b_sl15_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y (self-report)") 
b_sl16_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y (self-report)") 

#-- Test for interaction by sex - no significant interactions
sl12_tanner_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y (actigraphy)") 
sl12_tanner_i
sl14_tanner_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="14Y (self-report)") 
sl14_tanner_i
sl15_tanner_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="15Y (self-report)") 
sl15_tanner_i
sl16_tanner_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="16Y (self-report)") 
sl16_tanner_i

#-- Rerun in full sample (adjusted for sex)
sl12_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y (actigraphy)") 
sl12_tanner
sl14_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="14Y (self-report)") 
sl14_tanner
sl15_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="15Y (self-report)") 
sl15_tanner
sl16_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="16Y (self-report)") 
sl16_tanner

#-- Put all the coefficients for tanner together, format sub-table
tanner_sl_crude <- rbind(g_sl12_tanner,g_sl14_tanner,g_sl15_tanner,g_sl16_tanner,
                      b_sl12_tanner,b_sl14_tanner,b_sl15_tanner,b_sl16_tanner,
                      sl12_tanner,sl14_tanner,sl15_tanner,sl16_tanner) %>%
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Sleep

g_sl12_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y (actigraphy)") 
g_sl14_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y (self-report)")
g_sl15_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y (self-report)") 
g_sl16_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y (self-report)")

b_sl12_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y (actigraphy)") 
b_sl14_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y (self-report)") 
b_sl15_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y (self-report)") 
b_sl16_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y (self-report)") 

#-- Test for interaction by sex - none significant.
sl12_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y (actigraphy)") 
sl12_pds_i
sl14_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="14Y (self-report)") 
sl14_pds_i
sl15_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="15Y (self-report)") 
sl15_pds_i
sl16_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y (self-report)")
sl16_pds_i

#-- Rerun in the full sample (adjusted for sex)
sl12_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y (actigraphy)") 
sl12_pds
sl14_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="14Y (self-report)") 
sl14_pds
sl15_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="15Y (self-report)") 
sl15_pds
sl16_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ puberty_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y (self-report)")
sl16_pds

#-- Put all the coefficients for tanner together, format sub-table
pds_sl_crude <- rbind(g_sl12_pds,g_sl14_pds,g_sl15_pds,g_sl16_pds,
                   b_sl12_pds,b_sl14_pds,b_sl15_pds,b_sl16_pds,
                   sl12_pds,sl14_pds,sl15_pds,sl16_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Sleep
g_sl12_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y (actigraphy)") 
g_sl14_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y (self-report)") 
g_sl15_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y (self-report)")
g_sl16_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y (self-report)") 

b_sl12_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y (actigraphy)")
b_sl14_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y (self-report)") 
b_sl15_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y (self-report)") 
b_sl16_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y (self-report)")

#-- Test interaction by sex - none significant.
sl12_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="12Y (actigraphy)")
sl12_aphv_i
sl14_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="14Y (self-report)") 
sl14_aphv_i
sl15_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="15Y (self-report)") 
sl15_aphv_i
sl16_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="16Y (self-report)")
sl16_aphv_i

#-- Rerun in the full sample (adjusting for sex)
sl12_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="12Y (actigraphy)")
sl12_aphv
sl14_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="14Y (self-report)") 
sl14_aphv
sl15_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="15Y (self-report)") 
sl15_aphv
sl16_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="16Y (self-report)")
sl16_aphv


#-- Put all the coefficients for tanner together, format sub-table
aphv_sl_crude <- rbind(g_sl12_aphv,g_sl14_aphv,g_sl15_aphv,g_sl16_aphv,
                    b_sl12_aphv,b_sl14_aphv,b_sl15_aphv,b_sl16_aphv,
                    sl12_aphv,sl14_aphv,sl15_aphv,sl16_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- Menarche & Sleep
g_sl12_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y (actigraphy)") 
g_sl14_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y (self-report)")
g_sl15_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y (self-report)") 
g_sl16_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ men_date_12y))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y (self-report)")

#-- Put all the coefficients for tanner together, format sub-table
m_sl_crude <- rbind(g_sl12_m,g_sl14_m,g_sl15_m,g_sl16_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_SLEEP_crude <- rbind(aphv_sl_crude, tanner_sl_crude, pds_sl_crude, m_sl_crude)



#-- BMI-adjusted relationships
#-- Tanner & Dep
g_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y (actigraphy)") 
g_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y (self-report)") 
g_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y (self-report)") 
g_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y (self-report)") 

b_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y (actigraphy)") 
b_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y (self-report)") 
b_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y (self-report)") 
b_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y (self-report)") 

#-- Test interaction between tanner & sex for ET sleep - not significant
sl12adj_tanner_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y (actigraphy)") 
sl12adj_tanner_i

#-- Rerun in full sample (adjusted for sex)
sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y (actigraphy)") 
sl12adj_tanner
sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="14Y (self-report)") 
sl14adj_tanner
sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="15Y (self-report)") 
sl15adj_tanner
sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="16Y (self-report)") 
sl16_tanner

#-- Put all the coefficients for tanner together, format sub-table
tanner_sl_adj <- rbind(g_sl12adj_tanner,g_sl14adj_tanner,g_sl15adj_tanner,g_sl16adj_tanner,
                    b_sl12adj_tanner,b_sl14adj_tanner,b_sl15adj_tanner,b_sl16adj_tanner,
                    sl12adj_tanner,sl14adj_tanner,sl15adj_tanner,sl16adj_tanner) %>%
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep

g_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y (actigraphy)") 
g_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y (self-report)")
g_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y (self-report)") 
g_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y (self-report)")

b_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y (actigraphy)") 
b_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y (self-report)") 
b_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y (self-report)") 
b_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y (self-report)") 

#-- Test for interaction by sex - none significant.
sl12adj_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y (actigraphy)") 
sl12adj_pds_i
sl14adj_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ puberty_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="14Y (self-report)") 
sl14adj_pds_i
sl15adj_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="15Y (self-report)") 
sl15adj_pds_i
sl16adj_pds_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ puberty_12y*as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y (self-report)")
sl16adj_pds_i

#-- Rerun in the full sample (adjusted for sex)
sl12adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y (actigraphy)") 
sl12adj_pds
sl14adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="14Y (self-report)") 
sl14adj_pds
sl15adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="15Y (self-report)") 
sl15adj_pds
sl16adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ puberty_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y (self-report)")
sl16adj_pds

#-- Put all the coefficients for tanner together, format sub-table
pds_sl_adj <- rbind(g_sl12adj_pds,g_sl14adj_pds,g_sl15adj_pds,g_sl16adj_pds,
                    b_sl12adj_pds,b_sl14adj_pds,b_sl15adj_pds,b_sl16adj_pds,
                    sl12adj_pds,sl14adj_pds,sl15adj_pds,sl16adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Dep
g_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y (actigraphy)") 
g_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y (self-report)") 
g_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y (self-report)")
g_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y (self-report)") 

b_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y (actigraphy)")
b_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y (self-report)") 
b_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y (self-report)") 
b_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y (self-report)") 

#-- Test interaction by sex - none significant.
sl12adj_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years*as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="12Y (actigraphy)")
sl12adj_aphv_i
sl14adj_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years*as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="14Y (self-report)") 
sl14adj_aphv_i
sl15adj_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years*as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="15Y (self-report)") 
sl15adj_aphv_i
sl16adj_aphv_i <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years*as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="16Y (self-report)")
sl16adj_aphv_i

#-- Rerun in the full sample (adjusting for sex)
sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="12Y (actigraphy)")
sl12adj_aphv
sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="14Y (self-report)") 
sl14adj_aphv
sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="15Y (self-report)") 
sl15adj_aphv
sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+as.factor(sex)+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="16Y (self-report)")
sl16adj_aphv

#-- Put all the coefficients for tanner together, format sub-table
aphv_sl_adj <- rbind(g_sl12adj_aphv,g_sl14adj_aphv,g_sl15adj_aphv,g_sl16adj_aphv,
                     b_sl12adj_aphv,b_sl14adj_aphv,b_sl15adj_aphv,b_sl16adj_aphv,
                     sl12adj_aphv,sl14adj_aphv,sl15adj_aphv,sl16adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- Menarche & Sleep
g_sl12adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y (actigraphy)") 
g_sl14adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y (self-report)")
g_sl15adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y (self-report)") 
g_sl16adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ men_date_12y+BMI_7y))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y (self-report)")

#-- Put all the coefficients for tanner together, format sub-table
m_sl_adj <- rbind(g_sl12adj_m,g_sl14adj_m,g_sl15adj_m,g_sl16adj_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_SLEEP_adj <- rbind(aphv_sl_adj, tanner_sl_adj, pds_sl_adj, m_sl_adj)

#-- BMI and SES-adjusted relationships
#-- Tanner & Sleep
g_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="12Y (actigraphy)") 
g_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="14Y (self-report)") 
g_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="15Y (self-report)") 
g_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="16Y (self-report)") 

b_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="12Y (actigraphy)") 
b_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="14Y (self-report)") 
b_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="15Y (self-report)") 
b_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="16Y (self-report)") 

#-- Rerun in full sample (adjusted for sex)
sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="12Y (actigraphy)") 
sl12adj_tanner
sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="14Y (self-report)") 
sl14adj_tanner
sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="15Y (self-report)") 
sl15adj_tanner
sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="16Y (self-report)") 
sl16_tanner

#-- Put all the coefficients for tanner together, format sub-table
tanner_sl_adj <- rbind(g_sl12adj_tanner,g_sl14adj_tanner,g_sl15adj_tanner,g_sl16adj_tanner,
                       b_sl12adj_tanner,b_sl14adj_tanner,b_sl15adj_tanner,b_sl16adj_tanner,
                       sl12adj_tanner,sl14adj_tanner,sl15adj_tanner,sl16adj_tanner) %>%
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Sleep

g_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="12Y (actigraphy)") 
g_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="14Y (self-report)")
g_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="15Y (self-report)") 
g_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="16Y (self-report)")

b_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="12Y (actigraphy)") 
b_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="14Y (self-report)") 
b_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="15Y (self-report)") 
b_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="16Y (self-report)") 

#-- Rerun in the full sample (adjusted for sex)
sl12adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="12Y (actigraphy)") 
sl12adj_pds
sl14adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="14Y (self-report)") 
sl14adj_pds
sl15adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="15Y (self-report)") 
sl15adj_pds
sl16adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ puberty_12y+as.factor(sex)+ses_index_new_10pct+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="16Y (self-report)")
sl16adj_pds


#-- Put all the coefficients for tanner together, format sub-table
pds_sl_adj <- rbind(g_sl12adj_pds,g_sl14adj_pds,g_sl15adj_pds,g_sl16adj_pds,
                    b_sl12adj_pds,b_sl14adj_pds,b_sl15adj_pds,b_sl16adj_pds,
                    sl12adj_pds,sl14adj_pds,sl15adj_pds,sl16adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Sleep
g_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="12Y (actigraphy)") 
g_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="14Y (self-report)") 
g_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="15Y (self-report)")
g_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="16Y (self-report)") 

b_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="12Y (actigraphy)")
b_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="14Y (self-report)") 
b_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="15Y (self-report)") 
b_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="16Y (self-report)") 


#-- Rerun in the full sample (adjusting for sex)
sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+as.factor(sex)+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="12Y (actigraphy)")
sl12adj_aphv
sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+as.factor(sex)+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="14Y (self-report)") 
sl14adj_aphv
sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+as.factor(sex)+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="15Y (self-report)") 
sl15adj_aphv
sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+as.factor(sex)+ses_index_new_10pct+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="16Y (self-report)")
sl16adj_aphv

#-- Put all the coefficients for tanner together, format sub-table
aphv_sl_adj <- rbind(g_sl12adj_aphv,g_sl14adj_aphv,g_sl15adj_aphv,g_sl16adj_aphv,
                     b_sl12adj_aphv,b_sl14adj_aphv,b_sl15adj_aphv,b_sl16adj_aphv,
                     sl12adj_aphv,sl14adj_aphv,sl15adj_aphv,sl16adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- Menarche & Sleep
g_sl12adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="12Y (actigraphy)") 
g_sl14adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="14Y (self-report)")
g_sl15adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="15Y (self-report)") 
g_sl16adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ men_date_12y+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Menarche", y="16Y (self-report)")

#-- Put all the coefficients for tanner together, format sub-table
m_sl_adj <- rbind(g_sl12adj_m,g_sl14adj_m,g_sl15adj_m,g_sl16adj_m) %>% 
  filter(term=="men_date_12y" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_SLEEP_adj_bmi_ses <- rbind(aphv_sl_adj, tanner_sl_adj, pds_sl_adj, m_sl_adj)

#-- Rerun full-sample analysis with PROMIS depression scores

depsleepETadj_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+avgsleeptime_wd_hr+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y")
depsleepETadj_p
depsleep14adj_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+sleep_schoolday_cq14+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y")
depsleep14adj_p
depsleep15adj_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+sleep_weekend_child_cq15+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y")
depsleep15adj_p
depsleep16adj_p <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+sleep_weekday_cq16+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="17Y")
depsleep16adj_p

depsleep17_adj <- rbind(depsleepETadj_p,depsleep14adj_p,depsleep15adj_p,depsleep16adj_p) %>%
  filter( term %in% c("avgsleeptime_wd_hr","sleep_schoolday_cq14","sleep_weekend_child_cq15","sleep_weekday_cq16")) %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,term,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")
depsleep17_adj


write.csv(rep_cross_all_crude, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_depr_crude_age_fixed1.csv")
write.csv(rep_cross_all_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_depr_adjusted_fixed1.csv")
write.csv(rep_cross_all_adj_bmi_ses, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_depr_adjusted_BMI_SES_fixed1.csv")
write.csv(dep17_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/PROMIS_dep_puberty_full_sample.csv")

write.csv(rep_cross_all_SLEEP_crude, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_sleep_crude_age_fixed1.csv")
write.csv(rep_cross_all_SLEEP_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_sleep_adjusted_fixed1.csv")
write.csv(rep_cross_all_SLEEP_adj_bmi_ses, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_sleep_adjusted_BMI_SES_fixed1.csv")

write.csv(depsleep17_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/PROMIS_dep_sleep_puberty_full_sample.csv")


#----------------------------------------------------------
#-- Dep & Sleep, adjusted for PDS

g_dep12adj_actsleep <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Actigraphy sleep", y="12Y")
g_dep12adj_actsleep
g_dep14adj_actsleep <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Actigraphy sleep", y="14Y")
g_dep14adj_actsleep
g_dep15adj_actsleep <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Actigraphy sleep", y="15Y")
g_dep15adj_actsleep
g_dep16adj_actsleep <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Actigraphy sleep", y="16Y")
g_dep16adj_actsleep
g_dep17adj_actsleep <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Actigraphy sleep", y="17Y")
g_dep17adj_actsleep

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep12adj_actsleep <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Actigraphy sleep", y="12Y")
b_dep12adj_actsleep
b_dep14adj_actsleep <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Actigraphy sleep", y="14Y")
b_dep14adj_actsleep
b_dep15adj_actsleep <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Actigraphy sleep", y="15Y")
b_dep15adj_actsleep
b_dep16adj_actsleep <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ avgsleeptime_wd_hr +tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Actigraphy sleep", y="16Y")
b_dep16adj_actsleep
b_dep17adj_actsleep <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~avgsleeptime_wd_hr + tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Actigraphy sleep", y="17Y")
b_dep17adj_actsleep

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in total sample
dep12adj_actsleep <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_ET ~ avgsleeptime_wd_hr +tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Actigraphy sleep", y="12Y")
dep12adj_actsleep
dep14adj_actsleep <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_14y ~ avgsleeptime_wd_hr +tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Actigraphy sleep", y="14Y")
dep14adj_actsleep
dep15adj_actsleep <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_15y ~ avgsleeptime_wd_hr +tannerstg_12y+as.factor(sex)+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Actigraphy sleep", y="15Y")
dep15adj_actsleep
dep16adj_actsleep <- data.frame(summary(pool(with(data = d_imp, exp = lm(feeling_score_16y ~ avgsleeptime_wd_hr +tannerstg_12y+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Actigraphy sleep", y="16Y")
dep16adj_actsleep
dep17adj_actsleep <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~avgsleeptime_wd_hr + tannerstg_12y+as.factor(sex)+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Actigraphy sleep", y="17Y")
dep17adj_actsleep

#-- Put all the coefficients for tanner together, format sub-table
dep_act_sleep_adj <- rbind(g_dep12adj_actsleep,g_dep14adj_actsleep,g_dep15adj_actsleep,g_dep16adj_actsleep,g_dep17adj_actsleep,
                           b_dep12adj_actsleep,b_dep14adj_actsleep,b_dep15adj_actsleep,b_dep16adj_actsleep,b_dep17adj_actsleep,
                           dep12adj_actsleep,dep14adj_actsleep,dep15adj_actsleep,dep16adj_actsleep,dep17adj_actsleep) %>%
  filter(term=="avgsleeptime_wd_hr" & y != "11Y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")
dep_act_sleep_adj
write.csv(dep_act_sleep_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Repeated_cross_sections_dep_act_sleep_adj1.csv")


# b_dep17_sleep_adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ sleep_weekday_cq16+puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="17Y") 
# b_dep17_sleep_adj_pds
# write.csv(b_dep17_sleep_adj_pds, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/b_dep17_sleep_adj_pds.csv")
# 
# 
# b_dep17_sleep_adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ mother_weight_12y_bin +puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="17Y") 
# b_dep17_sleep_adj_pds
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------------------------------------------------
# #-- FIGURES
# #-------------------------------------------------------------------------------------------------------------------------
# 
# #-- Adjusted relationships
# #-- Tanner & Dep
# g_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="Tanner", y="12Y (actigraphy)")
# g_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="Tanner", y="14Y (self-report)")
# g_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="Tanner", y="15Y (self-report)")
# g_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="Tanner", y="16Y (self-report)")
# 
# #-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
# 
# b_sl12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="Tanner", y="12Y (actigraphy)")
# b_sl14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="Tanner", y="14Y (self-report)")
# b_sl15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="Tanner", y="15Y (self-report)")
# b_sl16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="Tanner", y="16Y (self-report)")
# 
# #-- Put all the coefficients for tanner together, format sub-table
# tanner_sl_adj <- rbind(g_sl12adj_tanner,g_sl14adj_tanner,g_sl15adj_tanner,g_sl16adj_tanner,
#                        b_sl12adj_tanner,b_sl14adj_tanner,b_sl15adj_tanner,b_sl16adj_tanner) %>%
#   filter(term=="tannerstg_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
#   ) %>%
#   dplyr::select(sex,x,y,Beta_CI) %>%
#   pivot_wider(
#     names_from = c("y"),
#     values_from = "Beta_CI")
# 
# 
# #-- PDS & Dep
# 
# g_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="PDS", y="12Y (actigraphy)")
# g_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="PDS", y="14Y (self-report)")
# g_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="PDS", y="15Y (self-report)")
# g_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Girls", x="PDS", y="16Y (self-report)")
# 
# #-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
# 
# b_sl12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="12Y (actigraphy)")
# b_sl14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="14Y (self-report)")
# b_sl15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="15Y (self-report)")
# b_sl16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="PDS", y="16Y (self-report)")
# 
# #-- Put all the coefficients for tanner together, format sub-table
# pds_sl_adj <- rbind(g_sl12adj_pds,g_sl14adj_pds,g_sl15adj_pds,g_sl16adj_pds,
#                     b_sl12adj_pds,b_sl14adj_pds,b_sl15adj_pds,b_sl16adj_pds) %>%
#   filter(term=="puberty_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
#   ) %>%
#   dplyr::select(sex,x,y,Beta_CI) %>%
#   pivot_wider(
#     names_from = c("y"),
#     values_from = "Beta_CI")
# 
# 
# #-- APHV & Dep
# g_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Girls", x="APHV", y="12Y (actigraphy)")
# g_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Girls", x="APHV", y="14Y (self-report)")
# g_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Girls", x="APHV", y="15Y (self-report)")
# g_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Girls", x="APHV", y="16Y (self-report)")
# 
# #-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
# b_sl12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(avgsleeptime_wd_hr ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Boys", x="APHV", y="12Y (actigraphy)")
# b_sl14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_schoolday_cq14 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Boys", x="APHV", y="14Y (self-report)")
# b_sl15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_child_cq15 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Boys", x="APHV", y="15Y (self-report)")
# b_sl16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(sleep_weekday_cq16 ~ age_peak_velocity_years+BMI_7y))))) %>%
#   mutate(sex="Boys", x="APHV", y="16Y (self-report)")
# 
# #-- Put all the coefficients for tanner together, format sub-table
# aphv_sl_adj <- rbind(g_sl12adj_aphv,g_sl14adj_aphv,g_sl15adj_aphv,g_sl16adj_aphv,
#                      b_sl12adj_aphv,b_sl14adj_aphv,b_sl15adj_aphv,b_sl16adj_aphv) %>%
#   filter(term=="age_peak_velocity_years" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
#   ) %>%
#   dplyr::select(sex,x,y,Beta_CI) %>%
#   pivot_wider(
#     names_from = c("y"),
#     values_from = "Beta_CI")
# 
# #-- Menarche & Sleep
# g_sl12adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(avgsleeptime_wd_hr ~ men_date_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Menarche", y="12Y (actigraphy)")
# g_sl14adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_schoolday_cq14 ~ men_date_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Menarche", y="14Y (self-report)")
# g_sl15adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_child_cq15 ~ men_date_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Menarche", y="15Y (self-report)")
# g_sl16adj_m <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(sleep_weekday_cq16 ~ men_date_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Menarche", y="16Y (self-report)")
# 
# #-- Put all the coefficients for tanner together, format sub-table
# m_sl_adj <- rbind(g_sl12adj_m,g_sl14adj_m,g_sl15adj_m,g_sl16adj_m) %>%
#   filter(term=="men_date_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
#   ) %>%
#   dplyr::select(sex,x,y,Beta_CI) %>%
#   pivot_wider(
#     names_from = c("y"),
#     values_from = "Beta_CI")
# 
# rep_cross_all_SLEEP_adj <- rbind(aphv_sl_adj, tanner_sl_adj, pds_sl_adj, m_sl_adj)
# 
# #-------------------------------------------------------------------------------------------#
# #-- PROMIS depression scores & sleep (actigraphy @ ET visit & self-reported at 16y visit) --#
# 
# g_dep17act_slp <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ avgsleeptime_wd+c_age_days_comp_d_qu_12y+tannerstg_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Actigraphy sleep", y="17Y")
# b_dep17act_slp <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ avgsleeptime_wd+c_age_days_comp_d_qu_12y+puberty_12y+BMI_7y))))) %>%
#   mutate(sex="Boys", x="Actigraphy sleep", y="17Y")
# 
# g_dep17sr_slp <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ sleep_weekday_cq16+tannerstg_12y+c_age_days_comp_d_qu_12y+BMI_7y))))) %>%
#   mutate(sex="Girls", x="Actigraphy sleep", y="17Y")
# b_dep17sr_slp <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ sleep_weekday_cq16+puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
#   mutate(sex="Boys", x="Actigraphy sleep", y="17Y")
# 
# 
# 
# g_dep17act_slp <- data.frame(summary(pool(with(data = d_imp_g, exp = glm(promis_modsev ~ avgsleeptime_wd+c_age_days_comp_d_qu_12y+tannerstg_12y+BMI_7y, family=binomial))))) %>%
#   mutate(sex="Girls", x="Actigraphy sleep", y="17Y")
# b_dep17act_slp <- data.frame(summary(pool(with(data = d_imp_b, exp = glm(promis_modsev ~ avgsleeptime_wd+c_age_days_comp_d_qu_12y+puberty_12y+BMI_7y, family=binomial))))) %>%
#   mutate(sex="Boys", x="Actigraphy sleep", y="17Y")
# 
# g_dep17sr_slp <- data.frame(summary(pool(with(data = d_imp_g, exp = glm(promis_modsev ~ sleep_weekday_cq16+tannerstg_12y+c_age_days_comp_d_qu_12y+BMI_7y, family=binomial))))) %>%
#   mutate(sex="Girls", x="Actigraphy sleep", y="17Y")
# b_dep17sr_slp <- data.frame(summary(pool(with(data = d_imp_b, exp = glm(promis_modsev ~ sleep_weekday_cq16+puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y, family=binomial))))) %>%
#   mutate(sex="Boys", x="Actigraphy sleep", y="17Y")
# 
# 
# #-----------------------------------------------#
# #-- Organize the tabular results into figures --#
# 
# #-- Depression symptom outcome - adjusted for BMI & age
# tanner_dep_adj <- rbind(g_dep12adj_tanner,g_dep14adj_tanner,g_dep15adj_tanner,g_dep16adj_tanner,
#                         b_dep12adj_tanner,b_dep14adj_tanner,b_dep15adj_tanner,b_dep16adj_tanner) %>%
#   filter(term=="tannerstg_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Est=round(estimate,2),
#     SE=round(std.error,2)) %>%
#   dplyr::select(sex,x,y,Est,SE)
# 
# pds_dep_adj <- rbind(g_dep12adj_pds,g_dep14adj_pds,g_dep15adj_pds,g_dep16adj_pds,
#                      b_dep12adj_pds,b_dep14adj_pds,b_dep15adj_pds,b_dep16adj_pds) %>%
#   filter(term=="puberty_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Est=round(estimate,2),
#     SE=round(std.error,2)) %>%
#   dplyr::select(sex,x,y,Est,SE)
# 
# aphv_dep_adj <- rbind(g_dep12adj_aphv,g_dep14adj_aphv,g_dep15adj_aphv,g_dep16adj_aphv,
#                       b_dep12adj_aphv,b_dep14adj_aphv,b_dep15adj_aphv,b_dep16adj_aphv) %>%
#   filter(term=="age_peak_velocity_years" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Est=round(estimate,2),
#     SE=round(std.error,2)) %>%
#   dplyr::select(sex,x,y,Est,SE)
# 
# m_dep_adj <- rbind(g_dep12adj_m,g_dep14adj_m,g_dep15adj_m,g_dep16adj_m) %>%
#   filter(term=="men_date_12y" & y != "11Y") %>%
#   mutate(
#     sig<-case_when(p.value<0.01 ~ "***",
#                    p.value<0.05 ~ "**",
#                    p.value<0.1 ~ "*",
#                    T~""),
#     Est=round(estimate,2),
#     SE=round(std.error,2)) %>%
#   dplyr::select(sex,x,y,Est,SE)
# 
# fig1dat <- data.frame(rbind(tanner_dep_adj,pds_dep_adj,aphv_dep_adj,m_dep_adj))
# fig1dat_g <- fig1dat %>% filter(sex=="Girls")
# fig1dat_b <- fig1dat %>% filter(sex=="Boys")
# 
# library(ggplot2)
# library(gridExtra)
# library(ggpubr) # Combine figures
# 
# dodge <- position_dodge(width=0.5)
# fig1_g <- ggplot(fig1dat_g, aes(x=y, y=Est, colour=x)) +
#   geom_errorbar(aes(ymin=Est-1.96*SE, ymax=Est+1.96*SE), width=.1,position=dodge) +
#   geom_line(position=dodge) +
#   geom_point(position=dodge) +
#   geom_hline(yintercept=0) +
#   ylab(" Beta") +
#   xlab("Outcome") +
#   scale_y_continuous(limits=c(-1.2,1.7),
#                      breaks = c(-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,
#                                 -0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,
#                                 0.4,0.5,0.6,0.7,0.8,0.9,1,
#                                 1.1,1.2,1.3,1.4,1.5,1.6,1.7)) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) +
#   ggtitle("Girls: Associations between EPT & depression symptoms") +
#   theme_light()
# fig1_g
# 
# fig1_b <- ggplot(fig1dat_b, aes(x=y, y=Est, colour=x)) +
#   geom_errorbar(aes(ymin=Est-1.96*SE, ymax=Est+1.96*SE), width=.1,position=dodge) +
#   geom_line(position=dodge) +
#   geom_point(position=dodge) +
#   geom_hline(yintercept=0) +
#   ylab(" Beta") +
#   xlab("Outcome") +
#   scale_y_continuous(limits=c(-1.2,1.5),
#                      breaks = c(-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,
#                                 -0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,
#                                 0.4,0.5,0.6,0.7,0.8,0.9,1,
#                                 1.1,1.2,1.3,1.4,1.5,1.6,1.7)) +
#   scale_x_discrete(guide = guide_axis(angle = 25)) +
#   ggtitle("Boys: Associations between EPT & depression symptoms") +
#   theme_light()
# fig1_b
# 
# # #-- For girls, there's mediation of the impact of PDS on sleep duration by worsened depression symptoms
# # d_imp_g_1 <- complete(d_imp_g,3)
# # set.seed(23)
# # mediation.rb <- cmest(data = d_imp_g_1, 
# #                       model = "rb", 
# #                       outcome = "avgsleeptime_wd", 
# #                       exposure = "puberty_12y", 
# #                       mediator = c("feeling_score_ET"),  
# #                       EMint = F,
# #                       basec = c("BMI_7y"), 
# #                       mreg = list("linear"), 
# #                       yreg = "linear", 
# #                       a = 3.7, 
# #                       astar = 3, 
# #                       mval = list(1),
# #                       estimation = "imputation", 
# #                       inference = "bootstrap")
# # summary(mediation.rb)$summarydf
# # 
# # 
# # #-- For girls, 
# # d_imp_g_1 <- complete(d_imp_g,6)
# # 
# # with(data=d_imp_g, 
# #      exp=cmest(model = "rb", 
# #                outcome = "feeling_score_14y", 
# #                exposure = "puberty_12y", 
# #                mediator = c("avgsleeptime_wd"),  
# #                EMint = F,
# #                basec = c("BMI_7y"), 
# #                mreg = list("linear"), 
# #                yreg = "linear", 
# #                a = 3.7, 
# #                astar = 3, 
# #                mval = list(438),
# #                estimation = "imputation", 
# #                inference = "bootstrap"))
# # 
# # set.seed(23)
# # mediation.rb <- cmest(model = "rb", 
# #                       outcome = "feeling_score_14y", 
# #                       exposure = "puberty_12y", 
# #                       mediator = c("avgsleeptime_wd"),  
# #                       EMint = F,
# #                       basec = c("BMI_7y"), 
# #                       mreg = list("linear"), 
# #                       yreg = "linear", 
# #                       a = 3.7, 
# #                       astar = 3, 
# #                       mval = list(438),
# #                       estimation = "imputation", 
# #                       inference = "bootstrap")
# # summary(mediation.rb)$summarydf
# # 
# # 
# # 
# # #--- Regression estimated mediation -> THIS IS WHAT's REPORTED!!!!
# # set.seed(23)
# # mediation.rb <- cmest(data = da_b, 
# #                       model = "rb", 
# #                       outcome = "feeling_score_14y", 
# #                       exposure = "age_peak_velocity_years", 
# #                       mediator = c("avgsleeptime_wd"),  
# #                       EMint = F,
# #                       basec = c("BMI_7y"), 
# #                       mreg = list("linear"), 
# #                       yreg = "linear", 
# #                       a = 1, 
# #                       astar = 0, 
# #                       mval = list(1),
# #                       estimation = "imputation", 
# #                       inference = "bootstrap")
# # summary(mediation.rb)$summarydf
# # 
# # 
# # 
# # 
# # #-- Create a BMI-residualized APHV variable
# # 
# # 
# # Factors <- names(
# #   da %>% select(csection,
# #                 race_amind,race_asian,race_black,
# #                 race_hisp,race_morethan1,race_other,
# #                 female_d,
# #                 pcoll_grad,
# #                 MARITAL_QU7Y,WELFARE_QU7Y,
# #                 tannerstg_10y,tannerstg_11y,tannerstg_12y,tannerstg_14y,tannerstg_15y,tannerstg_16y,
# #                 FIRSTMENST_9y,FIRSTMENST_10y,firstmenst_11y))
# # da[Factors] <- lapply(da[Factors], as.numeric)
# # 
# # #-- Create age-residualized Tanner & PDS variables
# # summary(lm(tannerstg_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da))
# # da$tanner12_age_res <- da$tannerstg_12y-predict(lm(tannerstg_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da),da) 
# # summary(lm(puberty_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da))
# # da$pds12_age_res <- da$puberty_12y-predict(lm(puberty_12y ~ c_age_days_comp_d_qu_12y+BMI_7y,data=da),data=da)
# # 
# # #-- Split the data by sex
# # da_b <- da %>% filter(sex=="Male") 
# # da_g <- da %>% filter(sex=="Female") 
# # 
# # 
# # summary(lm(age_peak_velocity_years~BMI_7y,data=da_b))
# # da_b$age_peak_velocity_years_bmi <- da_b$age_peak_velocity_years-predict(lm(age_peak_velocity_years~BMI_7y,data=da_b),data=da_b)
# # summary(lm(age_peak_velocity_years~BMI_7y,data=da_g))
# # da_g$age_peak_velocity_years_bmi <- da_g$age_peak_velocity_years-predict(lm(age_peak_velocity_years~BMI_7y,data=da_g),data=da_g)
# # 
# # da <- da %>%
# #   mutate(
# #     age_peak_velocity_years_bmi = case_when(
# #       female_d==1 ~ age_peak_velocity_years-(14.35-0.08*BMI_7y),
# #       female_d==2 ~ age_peak_velocity_years-(12.87-0.09*BMI_7y)
# #     )
# #   )
# 
# 
# plot(complete(d_imp_b)$feeling_score_16y, complete(d_imp_b)$C5_depr_tscore)
# 
# plot(complete(d_imp_g)$feeling_score_16y, complete(d_imp_g)$C5_depr_tscore)
