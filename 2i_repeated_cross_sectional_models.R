#--------------------------------------------------------------------------------------------------------------------#
#--- 2i_repeated_cross_sectional_models.R
#--- Date updated: 2/7/2024
#--- Description: Run linear models with one-time (dep) & repeated cross-sectional (sleep) outcomes (for supplement)
#--------------------------------------------------------------------------------------------------------------------#

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

#--------------------------------------------------------------------------------------#
#-- Crude relationships between one-time measure of depression symptoms (Y17) and PT --#
#----------------------------------------------------------------------------------------#

#-- Tanner & Dep
g_dep17_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="Y17") 
b_dep17_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="Y17")  

#-- Test interaction between tanner & sex at the 17Y visit (not signif, p=0.6693) - can run together, adjusting for sex
dep17_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="Y17")
dep17_tanner

dep17_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="Y17")
dep17_tanner


#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_crude <- rbind(g_dep17_tanner,b_dep17_tanner,dep17_tanner) %>%
  filter(term=="tannerstg_12y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    #-- If APHV, multiply estimate & CI bounds by -1 (for Younger APHV, for congruence with the other 2 measures)
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep17_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="Y17")
b_dep17_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="Y17")

#-- Check interaction between PDS & sex at the 17y visit - not significant (p-value 0.8198)
dep17_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y*as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="Y17")
dep17_pds
dep17_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y + as.factor(sex)+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="Y17")
dep17_pds

#-- Put all the coefficients for tanner together, format sub-table
pds_dep_crude <- rbind(g_dep17_pds,b_dep17_pds,dep17_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    #-- If APHV, multiply estimate & CI bounds by -1 (for Younger APHV, for congruence with the other 2 measures)
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

#-- APHV & Dep (adjusting for age at 12Y visit doesn't change associations for APHV)
g_dep17_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years))))) %>%
  mutate(sex="Girls", x="APHV", y="Y17")
b_dep17_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years))))) %>%
  mutate(sex="Boys", x="APHV", y="Y17")

#--Check interaction between APHV & sex at the 17 visit - not significant (p-vale =0.8024)
dep17_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years*as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="Y17")
dep17_aphv
dep17_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+as.factor(sex)))))) %>%
  mutate(sex="All", x="APHV", y="Y17")
dep17_aphv

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_crude <- rbind(g_dep17_aphv,b_dep17_aphv,dep17_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    #-- If APHV, multiply estimate & CI bounds by -1 (for Younger APHV, for congruence with the other 2 measures)
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


rep_cross_all_crude <- rbind(aphv_dep_crude, tanner_dep_crude, pds_dep_crude) %>%
  group_by(x) %>%
    pivot_wider(
      names_from=c(sex),
      values_from=c(Y17)
    ) %>%
  mutate(
    model="Model 1"
  ) %>%
  arrange(x)


#----------------------------------------------------------------------------------------#
#-- Model 2: BMI-adjusted relationships
#-- Tanner & Dep
g_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="Y17") 
b_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="Y17") 
dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+sex+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="Y17") 

#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_adj <- rbind(g_dep17adj_tanner,b_dep17adj_tanner,dep17adj_tanner) %>% 
  filter(term=="tannerstg_12y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="Y17")
b_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="Y17") 
dep17adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+sex+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="Y17") 

#-- Put all the coefficients for tanner together, format sub-table
pds_dep_adj <- rbind(g_dep17adj_pds,b_dep17adj_pds,dep17adj_pds) %>% 
  filter(term=="puberty_12y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Dep
g_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Girls", x="APHV", y="Y17") 
b_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y))))) %>%
  mutate(sex="Boys", x="APHV", y="Y17") 
dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+sex+BMI_7y))))) %>%
  mutate(sex="All", x="APHV", y="Y17")

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_adj <- rbind(g_dep17adj_aphv,b_dep17adj_aphv,dep17adj_aphv) %>% 
  filter(term=="age_peak_velocity_years") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_adj <- rbind(aphv_dep_adj, tanner_dep_adj, pds_dep_adj) %>%
  group_by(x) %>%
  pivot_wider(
    names_from=c(sex),
    values_from=c(Y17)
  ) %>%
  mutate(
    model="Model 2"
  ) %>%
  arrange(x)

#----------------------------------------------------------------------------------------#
#-- Model 3: BMI and SES adjusted relationships
#-- Tanner & Dep
g_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Tanner", y="Y17") 
b_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Tanner", y="Y17") 
dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ tannerstg_12y+sex+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="Tanner", y="Y17") 

#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_adj <- rbind(g_dep17adj_tanner,b_dep17adj_tanner,dep17adj_tanner) %>% 
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- PDS & Dep
g_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS", y="Y17")
b_dep17adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS", y="Y17") 
dep17adj_pds <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ puberty_12y+sex+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="All", x="PDS", y="Y17") 

#-- Put all the coefficients for tanner together, format sub-table
pds_dep_adj <- rbind(g_dep17adj_pds,b_dep17adj_pds,dep17adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- APHV & Dep
g_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="APHV", y="Y17") 
b_dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="APHV", y="Y17")
dep17adj_aphv <- data.frame(summary(pool(with(data = d_imp, exp = lm(C5_depr_tscore ~ age_peak_velocity_years+sex+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="All", x="APHV", y="Y17")

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_adj <- rbind(g_dep17adj_aphv,b_dep17adj_aphv,dep17adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                  p.value<0.05 ~ "*",
                  p.value<0.1 ~ ".",
                  T~""),
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_adj_bmi_ses <- rbind(aphv_dep_adj, tanner_dep_adj, pds_dep_adj) %>%
  group_by(x) %>%
  pivot_wider(
    names_from=c(sex),
    values_from=c(Y17)
  ) %>%
  mutate(
    model="Model 3"
  ) %>%
  arrange(x)

table2 <- rbind(rep_cross_all_crude, rep_cross_all_adj, rep_cross_all_adj_bmi_ses) %>%
  arrange(x)
write.csv(table2, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/table2.csv")

#---------------------------------------------------------------------------------------------------#
#-- Look at relationships with sleep duration
#-- Table 3 includes results for 12Y actigraphy for the full sample only, by sex to the supplement
#---------------------------------------------------------------------------------------------------#

#-------------------------#
#-- Crude relationships --#
#-------------------------#

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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_SLEEP_crude <- rbind(aphv_sl_crude, tanner_sl_crude, pds_sl_crude) %>% 
  mutate(
    model="Model 1"
  ) %>%
  arrange(x)



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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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

#-- Put all the coefficients for APHV together, format sub-table
aphv_sl_adj <- rbind(g_sl12adj_aphv,g_sl14adj_aphv,g_sl15adj_aphv,g_sl16adj_aphv,
                     b_sl12adj_aphv,b_sl14adj_aphv,b_sl15adj_aphv,b_sl16adj_aphv,
                     sl12adj_aphv,sl14adj_aphv,sl15adj_aphv,sl16adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")

rep_cross_all_SLEEP_adj <- rbind(aphv_sl_adj, tanner_sl_adj, pds_sl_adj) %>% 
  mutate(
    model="Model 2"
  ) %>%
  arrange(x)

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
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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


#-- Put all the coefficients for PDS together, format sub-table
pds_sl_adj <- rbind(g_sl12adj_pds,g_sl14adj_pds,g_sl15adj_pds,g_sl16adj_pds,
                    b_sl12adj_pds,b_sl14adj_pds,b_sl15adj_pds,b_sl16adj_pds,
                    sl12adj_pds,sl14adj_pds,sl15adj_pds,sl16adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
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

#-- Put all the coefficients for APHV together, format sub-table
aphv_sl_adj <- rbind(g_sl12adj_aphv,g_sl14adj_aphv,g_sl15adj_aphv,g_sl16adj_aphv,
                     b_sl12adj_aphv,b_sl14adj_aphv,b_sl15adj_aphv,b_sl16adj_aphv,
                     sl12adj_aphv,sl14adj_aphv,sl15adj_aphv,sl16adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,y,Beta_CI) %>%
  pivot_wider(
    names_from = c("y"),
    values_from = "Beta_CI")


#-- Stack all the results for M3 results linking sleep across time with PT measures at ET
rep_cross_all_SLEEP_adj_bmi_ses <- rbind(aphv_sl_adj, tanner_sl_adj, pds_sl_adj) %>% 
  mutate(
    model="Model 3"
  ) %>%
  arrange(x)


table3_aphv <- rbind(rep_cross_all_SLEEP_crude,rep_cross_all_SLEEP_adj,rep_cross_all_SLEEP_adj_bmi_ses)
names(table3_aphv) <- c("sex","x","Y12_act","Y14_SR","Y15_SR","Y16_SR","model")

#- For the supplement (Supplemental Table S.1), keep sex-specific estimates
sup_table3_aphv <- table3_aphv %>%
  #-- Drop the repeated cross-sectional results for self-reported sleep (analyzed properly with linear mixed effects models - this was just to look at patterns)
  dplyr::select(sex,x,model,Y12_act) %>%
  arrange(x,desc(sex))

#- For main text Table 3, only keep full-sample estimates for APHV
table3_aphv <- table3_aphv %>%
  filter(sex=="All") %>%
  dplyr::select(sex,x,model,Y12_act) %>%
  arrange(x)

write.csv(sup_table3_aphv, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/sup_table3_aphv.csv")
write.csv(table3_aphv, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/table3_aphv.csv")
