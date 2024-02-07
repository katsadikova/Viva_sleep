#-----------------------------------------------------------------------------------------------------------#
#--- 3i_linear_mixed_effects.R
#--- Date: 4/30/2023 (original)
#--- Date updated: 2/7/2024
#--- Description: Look at longitudinal trends in self-reported sleep duration by time (relative to ET)
#-----------------------------------------------------------------------------------------------------------#

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

summary(d_pre_imp$tannerstg_12y)

summary(d_pre_imp$sleep_schoolday_cq14) #missing 755/1138
summary(d_pre_imp$sleep_weekday_child_cq15) #missing 684/1138
summary(d_pre_imp$sleep_weekday_cq16) #missing 454/1138

#----------------------------------------------------------------------------------------#
#-- Make edits to the mids object d_imp - might need later

#--- convert to long-from data set to add the "max_problems" variable
long <- complete(d_imp, action="long", include=TRUE)
long$ses_index_new_10pct <- ifelse(long$ses_index_new_03<=7,1,0)
long$miss_actigraphy <-ifelse(is.na(d_pre_imp$avgsleeptime),1,0)
long$miss_ET <-ifelse(is.na(d_pre_imp$c_age_days_comp_d_qu_12y),1,0)
long$tannerstg_12y <- as.numeric(long$tannerstg_12y)
long$men_date_12y <- long$men_date_12y/365.25
long$avgsleeptime_wd <- long$avgsleeptime_wd/60
long$avgsleeptime <- long$avgsleeptime/60
long$overweight <- ifelse(long$BMI_7y>=25,1,0)
long$obese <- ifelse(long$BMI_7y>=30,1,0)
long$BMI_7y_z <- (long$BMI_7y-17.18)/2.85

names(long)

#-- Rename the age & sleep variables
names(long)
long$age7 = long$c_age_days_COMP_D_QU7Y
long$age9 = long$c_age_days_comp_d_qu_9y
long$age10 = long$c_age_days_comp_d_qu_10y
long$age11 = long$c_age_days_comp_d_qu_11y
long$age12 = long$c_age_days_comp_d_qu_12y
long$age14 = long$age_days_cq14
long$age15 = long$age_days_cq15
long$age16 = long$age_days_cq16
long$sr_wd_sleep14 = long$sleep_schoolday_cq14
long$sr_wd_sleep15 = long$sleep_weekday_child_cq15
long$sr_wd_sleep16 = long$sleep_schoolday_cq14
long$sr_wd_sleep14_8h = ifelse(long$sleep_schoolday_cq14<8,1,0)
long$sr_wd_sleep15_8h = ifelse(long$sleep_weekday_child_cq15<8,1,0)
long$sr_wd_sleep16_8h = ifelse(long$sleep_schoolday_cq14<8,1,0)
names(long)


#######################################################################################################
#-- Run LME models with repeated sleep measures
#-- Method resource: https://stats.stackexchange.com/questions/117605/lmer-with-multiply-imputed-data 

#-- Pivot data such that there's a row for every visit
ages <- long %>% 
  select(.imp,.id,aid,sex,age_peak_velocity_years,age11,age12,age14,age15,age16) %>%
  pivot_longer(
    cols = age11:age16,
    names_to = c("visit"),
    values_to = "age") %>%
  mutate(
    age = age/365.25,
    time_rt_aphv = round(age-age_peak_velocity_years,1),
    time_rt_avg_aphv = case_when(
      sex=="Female" ~ round(age-11.2,1),
      sex=="Male" ~ round(age-13.1,1)
    ),
    time_rt_avg_12y = age-13.3,
    time_rt_avg_14y = age-14.5,
    time_rt_avg_15y = age-15.6,
    sq_time_rt_avg_12y = time_rt_avg_12y*time_rt_avg_12y,
    visit = case_when(visit=="age11" ~ 11,
                     visit=="age12" ~ 12,
                     visit=="age14" ~ 14,
                     visit=="age15" ~ 15,
                     visit=="age16" ~ 16)) %>%
  select(-c(age_peak_velocity_years,sex))

# 
# pds <- long %>% 
#   select(.imp,.id,aid,puberty_11y,puberty_12y,puberty_14y) %>%
#   pivot_longer(
#     cols = puberty_11y:puberty_14y,
#     names_to = c("visit"),
#     values_to = "pds") %>%
#   mutate(
#     visit = case_when(visit=="puberty_11y" ~ 11,
#                       visit=="puberty_12y" ~ 12,
#                       visit=="puberty_14y" ~ 14))
# 
# feelings <- long %>% 
#   select(.imp,.id,aid,feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y) %>%
#   pivot_longer(
#     cols = feeling_score_11y:feeling_score_16y,
#     names_to = c("visit"),
#     values_to = "feeling_score") %>%
#   mutate(
#     visit = case_when(visit=="feeling_score_11y" ~ 11,
#                      visit=="feeling_score_ET" ~ 12,
#                      visit=="feeling_score_14y" ~ 14,
#                      visit=="feeling_score_15y" ~ 15,
#                      visit=="feeling_score_16y" ~ 16)) 
# 
# ages_feelings <- merge(ages,feelings,by=c(".imp",".id","aid","visit"))
# 
# other_da <- long %>%
#   select(-c(age11,age12,age14,age15,age16,
#             feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y,
#             sr_wd_sleep14,sr_wd_sleep15,sr_wd_sleep16)) 
# 
# all_da_long <- merge(other_da,ages_feelings,by=c(".imp",".id","aid")) %>%
#   mutate(.id = as.numeric(paste0(.id,visit)))
# 
# d_imp_long <- as.mids(all_da_long) 
# #-- 5/16 - remove the age 11 timepoint
# b_da_long <- d_imp_long %>% filter(sex=="Male" & visit>11)
# g_da_long <- d_imp_long %>% filter(sex=="Female" & visit>11)
# 
# summary(complete(b_da_long)$time_rt_avg_12y)
# 
# save(b_da_long, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/b_da_long.Rda")
# save(g_da_long, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/g_da_long.Rda")

sleep <- long %>%
  select(.imp,.id,aid,avgsleeptime_wd,sr_wd_sleep14,sr_wd_sleep15,sr_wd_sleep16) %>%
  mutate(wdsl12=avgsleeptime_wd,
         wdsl14=sr_wd_sleep14,
         wdsl15=sr_wd_sleep15,
         wdsl16=sr_wd_sleep16) %>%
  select(.imp,.id,aid,wdsl12,wdsl14,wdsl15,wdsl16) %>%
  pivot_longer(
    cols = wdsl12:wdsl16,
    names_to = c("visit"),
    values_to = "sleep") %>%
  mutate(
    visit = case_when(visit=="wdsl12" ~ 12,
                      visit=="wdsl14" ~ 14,
                      visit=="wdsl15" ~ 15,
                      visit=="wdsl16" ~ 16),
    sleep8 = ifelse(sleep<8,1,0)) %>%
  group_by(.imp,.id) %>%
  mutate(
    sleep_lag = lag(sleep)
  ) %>%
  ungroup()

ages_sleep<-merge(ages,sleep,by=c(".imp",".id","aid","visit"),all.x=TRUE)

other_da <- long %>%
  select(-c(age11,age12,age14,age15,age16,
            feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y,
            sr_wd_sleep14,sr_wd_sleep15,sr_wd_sleep16)) 

all_da_long <- merge(other_da,ages_sleep,by=c(".imp",".id","aid")) %>%
  mutate(.id = as.numeric(paste0(.id,visit))) 

d_imp_long <- as.mids(all_da_long) 

#-- Limit the long data set to 14,15,16Y visits (self-reported sleep)
b_da_long <- d_imp_long %>% filter(sex=="Male" & visit>12)
g_da_long <- d_imp_long %>% filter(sex=="Female" & visit>12)
da_long <- d_imp_long %>% filter(visit>12)

#--------------------------------------------------------------------------------------------------------------
#-- Model 1: CRUDE sleep & pubertal timing (14,15,16 - self-reported only)

#----------#
#-- APHV --#

#-- Girls
lmm1a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm1a_gi
lmm1a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm1a_g
#-- Boys
lmm1a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm1a_bi
lmm1a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm1a_b
#-- All
lmm1ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm1ai
lmm1a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + as.factor(sex) +  c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm1a

#------------#
#-- Tanner --#

#-- Girls
lmm1t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm1t_gi
lmm1t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm1t_g
#-- Boys
lmm1t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm1t_bi
lmm1t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm1t_b
#-- All
lmm1ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm1ti
lmm1t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 + tannerstg_12y| aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm1t

#---------#
#-- PDS --#

#-- Girls
lmm1p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm1p_gi
lmm1p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm1p_g
#-- Boys
lmm1p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm1p_bi
lmm1p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm1p_b
#-- All
lmm1pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm1pi
lmm1p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 + puberty_12y | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm1p


lme_sleep_SR_crude <- rbind(lmm1a_g, lmm1a_b, lmm1a, 
                            lmm1t_g, lmm1t_b, lmm1t,
                            lmm1p_g, lmm1p_b, lmm1p) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
         term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = case_when(
      x=="APHV" ~ paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig),
      T ~ paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
    ),
    model="Model 1"
  ) %>%
  dplyr::select(sex,x,term,Beta_CI,model) %>%
  arrange(sex,x)


#------------------------------------------------------------------------------
#-- Model 2: ADJUSTED sleep & pubertal timing (14,15,16 - self-reported only)

lmm2a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm2a_gi
lmm2a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm2a_g
lmm2a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm2a_bi
lmm2a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm2a_b
lmm2ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm2ai
lmm2a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm2a

lmm2t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm2t_gi
lmm2t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 + tannerstg_12y| aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm2t_g
lmm2t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm2t_bi
lmm2t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm2t_b
lmm2ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm2ti
lmm2t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm2t

lmm2p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm2p_gi
lmm2p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm2p_g
lmm2p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm2p_bi
lmm2p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm2p_b
lmm2pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm2pi
lmm2p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm2p


lme_sleep_SR_adj <- rbind(lmm2a_g, lmm2a_b, lmm2a,
                          lmm2t_g, lmm2t_b, lmm2t,
                          lmm2p_g, lmm2p_b, lmm2p) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
           term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = case_when(
      x=="APHV" ~ paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig),
      T ~ paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
    ),
    model="Model 2"
  ) %>%
  dplyr::select(sex,x,term,Beta_CI,model) %>%
  arrange(sex,x)

#--------------------------------------------------------------------------------------
#-- Model 3: BMI & SES ADJUSTED sleep & pubertal timing (14,15,16 - self-reported only)

lmm3a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm3a_gi
lmm3a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm3a_g
lmm3a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm3a_bi
lmm3a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm3a_b
lmm3ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm3ai
lmm3a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm3a

lmm3t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm3t_gi
lmm3t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm3t_g
lmm3t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm3t_bi
lmm3t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm3t_b
lmm3ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm3ti
lmm3t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm3t

lmm3p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm3p_gi
lmm3p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm3p_g
lmm3p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm3p_bi
lmm3p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm3p_b
lmm3pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+  (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm3pi
lmm3p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm3p

lme_sleep_SR_adj_bmi_ses <- rbind(lmm3a_g, lmm3a_b, lmm3a,
                                  lmm3t_g, lmm3t_b, lmm3t,
                                  lmm3p_g, lmm3p_b, lmm3p) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
           term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male" & term!="ses_index_new_10pct") %>%
  mutate(
    sig=case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   p.value<0.1 ~ ".",
                   T~""),
    Beta_CI = case_when(
      x=="APHV" ~ paste0(-1*round(estimate,2),"(",-1*round((estimate+1.96*std.error),2),",",-1*round((estimate-1.96*std.error),2),")",sig),
      T ~ paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
    ),
    model="Model 3"
  ) %>%
  dplyr::select(sex,x,term,Beta_CI,model) %>%
  arrange(sex,x)

#-- For Table 3 in manuscript
table3_SR <- rbind(lme_sleep_SR_crude,lme_sleep_SR_adj,lme_sleep_SR_adj_bmi_ses) %>%
  filter(sex=="All") %>%
  dplyr::select(-term) %>%
  arrange(x)

#-- For Supplement Table 1
sup_table3_SR <- rbind(lme_sleep_SR_crude,lme_sleep_SR_adj,lme_sleep_SR_adj_bmi_ses) %>%
  dplyr::select(-term) %>%
  arrange(x)

write.csv(table3_SR, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/table3_SR.csv")
write.csv(sup_table3_SR, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/sup_table3_SR.csv")
