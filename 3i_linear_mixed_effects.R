#-----------------------------------------------------------------------------------------------------------#
#--- 5i_mixed_effects.R
#--- Date:    4/30/2023 (original)
#--- Updated: 7/15/2023
#--- Description: Look at longitudinal trends in dep & self-reported sleep duration by time (relative to ET)
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

#-- Pivot da such that there's a row for every repeated dep entry
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


pds <- long %>% 
  select(.imp,.id,aid,puberty_11y,puberty_12y,puberty_14y) %>%
  pivot_longer(
    cols = puberty_11y:puberty_14y,
    names_to = c("visit"),
    values_to = "pds") %>%
  mutate(
    visit = case_when(visit=="puberty_11y" ~ 11,
                      visit=="puberty_12y" ~ 12,
                      visit=="puberty_14y" ~ 14))

feelings <- long %>% 
  select(.imp,.id,aid,feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y) %>%
  pivot_longer(
    cols = feeling_score_11y:feeling_score_16y,
    names_to = c("visit"),
    values_to = "feeling_score") %>%
  mutate(
    visit = case_when(visit=="feeling_score_11y" ~ 11,
                     visit=="feeling_score_ET" ~ 12,
                     visit=="feeling_score_14y" ~ 14,
                     visit=="feeling_score_15y" ~ 15,
                     visit=="feeling_score_16y" ~ 16)) 

ages_feelings <- merge(ages,feelings,by=c(".imp",".id","aid","visit"))

other_da <- long %>%
  select(-c(age11,age12,age14,age15,age16,
            feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y,
            sr_wd_sleep14,sr_wd_sleep15,sr_wd_sleep16)) 

all_da_long <- merge(other_da,ages_feelings,by=c(".imp",".id","aid")) %>%
  mutate(.id = as.numeric(paste0(.id,visit)))

d_imp_long <- as.mids(all_da_long) 
#-- 5/16 - remove the age 11 timepoint
b_da_long <- d_imp_long %>% filter(sex=="Male" & visit>11)
g_da_long <- d_imp_long %>% filter(sex=="Female" & visit>11)

summary(complete(b_da_long)$time_rt_avg_12y)

save(b_da_long, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/b_da_long.Rda")
save(g_da_long, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/data/g_da_long.Rda")

names(complete(b_da_long))
#-- Look at age relative to APHV as the time scale
#-- Linear mixed effects models 
#-- Good resource: https://stats.stackexchange.com/questions/117605/lmer-with-multiply-imputed-data 

#-- Model 1: APHV and age 
lmm1a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm1a_gi
lmm1a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm1a_g
lmm1a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm1a_bi
lmm1a_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm1a_b

#-- Model 1: Tanner and age 
lmm1t_gi = data.frame(summary(est <- pool( with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm1t_gi
lmm1t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm1t_g
lmm1t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm1t_bi
lmm1t_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm1t_b

#-- Model 1: PDS and age 
lmm1p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm1p_gi
lmm1p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm1p_g
lmm1p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm1p_bi
lmm1p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm1p_b

#-- Model 1: Age at menarche and age 
lmm1m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+time_rt_avg_12y+men_date_12y*time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="Yes")
lmm1m_gi
lmm1m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="No")
lmm1m_g

lme_crude <- rbind(lmm1a_g,lmm1t_g,lmm1p_g,lmm1m_g,lmm1a_b,lmm1t_b,lmm1p_bi) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_crude, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Dep_crude_fixed.csv")

#-- For the alternative version of Table 2, export coefficients for main & interaction terms
lme_crude_ints <- rbind(lmm1a_gi,lmm1t_gi,lmm1p_gi,lmm1m_gi,
                        lmm1a_bi,lmm1t_bi,lmm1p_bi) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="c_age_days_comp_d_qu_12y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_crude_ints, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Dep_crude_INTERACTIONS.csv")


#-- Model 2: APHV, age, BMI  (no interaction with time)
lmm2a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm2a_gi
lmm2a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm2a_g
lmm2a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm2a_bi
lmm2a_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm2a_b

#-- Model 2: Tanner, age, BMI (interaction not significant)
lmm2t_gi = data.frame(summary(est <- pool( with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm2t_gi
lmm2t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm2t_g
lmm2t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm2t_bi
lmm2t_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm2t_b

#-- Model 2: PDS, age, BMI (interaction significant among boys (p-value 0.04)!!!!)
lmm2p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm2p_gi
lmm2p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm2p_g
lmm2p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm2p_bi
lmm2p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm2p_b

#-- Model 2: Age at menarche and age 
lmm2m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+BMI_7y+time_rt_avg_12y+men_date_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="Yes")
lmm2m_gi
lmm2m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="No")
lmm2m_g

lme_bmi_adj <- rbind(lmm2a_g,lmm2t_g,lmm2p_g,lmm2m_g,lmm2a_b,lmm2t_b,lmm2p_bi) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_bmi_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Dep_bmi_adj_fixed.csv")

#-- For the alternative version of Table 2, export coefficients for main & interaction terms
lme_bmi_adj_ints <- rbind(lmm2a_gi,lmm2t_gi,lmm2p_gi,lmm2m_gi,
                     lmm2a_bi,lmm2t_bi,lmm2p_bi) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y"
         & term!="c_age_days_comp_d_qu_12y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_bmi_adj_ints, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Dep_bmi_adj_INTERACTIONS.csv")


#-- Check if adding SES adds to the fit of the models - check among girls
test<-with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.pooled.AIC<-pool(test)$glanced$AIC
test.pooled.AIC
test.SES<-with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+ses_index_new_10pct+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.SES.pooled.AIC<-pool(test.SES)$glanced$AIC
test.SES.pooled.AIC 
aic_girls_aphv <- data.frame(cbind(test.pooled.AIC,test.SES.pooled.AIC))
write.csv(aic_girls_aphv, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/aic_girls_aphv.csv")

test<-with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.pooled.AIC<-pool(test)$glanced$AIC
test.pooled.AIC
test.SES<-with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+ses_index_new_10pct+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.SES.pooled.AIC<-pool(test.SES)$glanced$AIC
test.SES.pooled.AIC 
aic_girls_tanner <- data.frame(cbind(test.pooled.AIC,test.SES.pooled.AIC))
write.csv(aic_girls_tanner, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/aic_girls_tanner.csv")

test<-with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.pooled.AIC<-pool(test)$glanced$AIC
test.pooled.AIC
test.SES<-with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+ses_index_new_10pct+BMI_7y+time_rt_avg_12y + (1 | aid)))
test.SES.pooled.AIC<-pool(test.SES)$glanced$AIC
test.SES.pooled.AIC 
aic_girls_pds <- data.frame(cbind(test.pooled.AIC,test.SES.pooled.AIC))
write.csv(aic_girls_pds, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/aic_girls_pds.csv")

#--- Adding SES increases AIC values (poorer model fit) across all three pubertal timing measures (i.e. will only widen CIs, but not add to prediction accuracy)

#-- Model 3: APHV, age, BMI + SES
lmm3a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm3a_gi
lmm3a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm3a_g
lmm3a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm3a_bi
lmm3a_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm3a_b

#-- Model 3: Tanner, age, BMI, SES (interaction not significant)
lmm3t_gi = data.frame(summary(est <- pool( with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm3t_gi
lmm3t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm3t_g
lmm3t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit)))))%>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm3t_bi
lmm3t_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm3t_b

#-- Model 3: PDS, age, BMI, SES
lmm3p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm3p_gi
lmm3p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm3p_g
lmm3p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+puberty_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm3p_bi
#Check if confounded by mother making comments about weight
#No evidence of confounding, but strong independent predictor of greater depression symptoms
lmm3p_bi_mat = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ puberty_12y+mother_weight_12y_bin+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+puberty_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm3p_bi_mat

#-- Model 3: Age at menarche, BMI, age, SES 
lmm3m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+men_date_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="Yes")
lmm3m_gi
lmm3m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ men_date_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="No")
lmm3m_g


lme_bmi_ses_adj <- rbind(lmm3a_g,lmm3t_g,lmm3p_g,lmm3m_g,lmm3a_b,lmm3t_b,lmm3p_bi) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" 
         & term!="ses_index_new_10pct" & term!="c_age_days_comp_d_qu_12y") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_bmi_ses_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Dep_bmi_ses_aadj.csv")


#######################################################################################################
#-- Run LME models with repeated sleep measures


#-- Construct an alternative long data set

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

ages_feelings_sleep<-merge(ages_feelings,sleep,by=c(".imp",".id","aid","visit"),all.x=TRUE)

other_da <- long %>%
  select(-c(age11,age12,age14,age15,age16,
            feeling_score_11y,feeling_score_ET,feeling_score_14y,feeling_score_15y,feeling_score_16y,
            sr_wd_sleep14,sr_wd_sleep15,sr_wd_sleep16)) 

all_da_long <- merge(other_da,ages_feelings_sleep,by=c(".imp",".id","aid")) %>%
  mutate(.id = as.numeric(paste0(.id,visit))) 

d_imp_long <- as.mids(all_da_long) 

#-- Limit the long data set to 14,15,16Y visits (self-reported sleep)
b_da_long <- d_imp_long %>% filter(sex=="Male" & visit>12)
g_da_long <- d_imp_long %>% filter(sex=="Female" & visit>12)
da_long <- d_imp_long %>% filter(visit>12)

#--------------------------------------------------------------------------------------------------------------
#-- Model 4: CRUDE sleep & pubertal timing (14,15,16 - self-reported only)
#-- Girls
lmm4a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm4a_gi
lmm4a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm4a_g
#-- Boys
lmm4a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm4a_bi
lmm4a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm4a_b
#-- All
lmm4ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm4ai
lmm4a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y + as.factor(sex) +  c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm4a

#-- Girls
lmm4t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm4t_gi
lmm4t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm4t_g
#-- Boys
lmm4t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm4t_bi
lmm4t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm4t_b
#-- All
lmm4ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm4ti
lmm4t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 + tannerstg_12y| aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm4t

#-- Girls
lmm4p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm4p_gi
lmm4p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm4p_g
#-- Boys
lmm4p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm4p_bi
lmm4p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm4p_b
#-- All
lmm4pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm4pi
lmm4p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 + puberty_12y | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm4p

lmm4m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+time_rt_avg_12y+men_date_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm4m_gi
lmm4m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm4m_g


lme_sleep_SR_crude <- rbind(lmm4a_g, lmm4a_b, lmm4a, 
                            lmm4t_g, lmm4t_b, lmm4t,
                            lmm4p_g, lmm4p_b, lmm4p,
                            lmm4m_g) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
         term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_sleep_SR_crude, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Sleep_all_puberty_CRUDE_fixed.csv")


#------------------------------------------------------------------------------
#-- Model 5: ADJUSTED sleep & pubertal timing (14,15,16 - self-reported only)

lmm5a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm5a_gi
lmm5a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm5a_g
lmm5a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm5a_bi
lmm5a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm5a_b
lmm5ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm5ai
lmm5a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm5a

lmm5t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm5t_gi
lmm5t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 + tannerstg_12y| aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm5t_g
lmm5t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm5t_bi
lmm5t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm5t_b
lmm5ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm5ti
lmm5t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm5t

lmm5p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm5p_gi
lmm5p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm5p_g
lmm5p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm5p_bi
lmm5p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm5p_b
lmm5pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm5pi
lmm5p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm5p


lmm5m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+time_rt_avg_12y+BMI_7y+men_date_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="Yes")
lmm5m_gi
lmm5m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="No")
lmm5m_g

lme_sleep_SR_adj <- rbind(lmm5a_g, lmm5a_b, lmm5a,
                          lmm5t_g, lmm5t_b, lmm5t,
                          lmm5p_g, lmm5p_b, lmm5p,
                          lmm5m_g) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
           term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_sleep_SR_adj, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Sleep_SR_adj_fixed.csv")


#--------------------------------------------------------------------------------------
#-- Model 6: BMI & SES ADJUSTED sleep & pubertal timing (14,15,16 - self-reported only)

lmm6a_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="Yes")
lmm6a_gi
lmm6a_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="APHV",interaction="No")
lmm6a_g
lmm6a_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="Yes")
lmm6a_bi
lmm6a_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="APHV",interaction="No")
lmm6a_b
lmm6ai = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+age_peak_velocity_years*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="Yes")
lmm6ai
lmm6a =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="APHV",interaction="No")
lmm6a

lmm6t_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="Yes")
lmm6t_gi
lmm6t_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Tanner",interaction="No")
lmm6t_g
lmm6t_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="Yes")
lmm6t_bi
lmm6t_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Tanner",interaction="No")
lmm6t_b
lmm6ti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+tannerstg_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="Yes")
lmm6ti
lmm6t =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Tanner",interaction="No")
lmm6t

lmm6p_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm6p_gi
lmm6p_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm6p_g
lmm6p_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm6p_bi
lmm6p_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm6p_b
lmm6pi = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(sleep ~ puberty_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+puberty_12y*time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+  (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm6pi
lmm6p =data.frame(summary(est <- pool( with(data=da_long, exp=lmer(sleep ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + as.factor(sex) + c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="No")
lmm6p


lmm6m_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+men_date_12y*time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="Yes")
lmm6m_gi
lmm6m_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(sleep ~ men_date_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Menarche",interaction="No")
lmm6m_g

lme_sleep_SR_adj_bmi_ses <- rbind(lmm6a_g, lmm6a_b, lmm6a,
                                  lmm6t_g, lmm6t_b, lmm6t,
                                  lmm6p_g, lmm6p_b, lmm6p,
                                  lmm6m_g) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" &
           term!="c_age_days_comp_d_qu_12y" & term!="as.factor(sex)Male" & term!="ses_index_new_10pct") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
write.csv(lme_sleep_SR_adj_bmi_ses, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_Sleep_SR_adj_bmi_ses_fixed.csv")

#------------------------------------------------------------------------


#-- Model 7: Associations between dep (14,15,16Y) & 12Y actigraphy sleep
lmm7act_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ avgsleeptime_wd*time_rt_avg_12y+tannerstg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="Yes")
lmm7act_gi
lmm7act_g =data.frame(summary(est <- pool( with(data=g_da_long, exp=lmer(feeling_score ~ avgsleeptime_wd+tannerstg_12y+time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="PDS",interaction="No")
lmm7act_g
# test <- complete(g_da_long,3)
# est <- lmer(feeling_score ~ avgsleeptime_wd+tannerstg_12y+time_rt_avg_12y + BMI_7y +
#                             (1 + avgsleeptime_wd | aid), na.action=na.omit, data=test)
# mean(ranef(est)$aid[,1])

lmm7act_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ avgsleeptime_wd*time_rt_avg_12y+tannerstg_12y  + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="Yes")
lmm7act_bi
lmm7act_b =data.frame(summary(est <- pool( with(data=b_da_long, exp=lmer(feeling_score ~ avgsleeptime_wd+tannerstg_12y+time_rt_avg_12y  + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1| aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="PDS",interaction="No")
lmm7act_b

#-- Both boys & girls
lmm7acti = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(feeling_score ~ avgsleeptime_wd*time_rt_avg_12y+tannerstg_12y + as.factor(sex) + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm7acti

lmm7act = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(feeling_score ~ avgsleeptime_wd+time_rt_avg_12y+tannerstg_12y + as.factor(sex) + BMI_7y + ses_index_new_10pct +ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="PDS",interaction="Yes")
lmm7act


lme_dep_actig_sleep <- rbind(lmm7act_g, lmm7act_b, lmm7act) %>%
  filter(term!="(Intercept)" & term!="time_rt_avg_12y" & term!="BMI_7y" & term!="ses_index_new_10pct" 
         & term!="c_age_days_comp_d_qu_12y" & term!="tannerstg_12y" & term!="as.factor(sex)Male") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "***",
                   p.value<0.05 ~ "**",
                   p.value<0.1 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
lme_dep_actig_sleep
write.csv(lme_dep_actig_sleep, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_lme_dep_actig_sleep.csv")

#-- Model 7: Associations between dep (14,15,16Y) & repeated contemporaneous sleep (14,15,16Y)
lmm7_gi = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y+sleep*time_rt_avg_12y+time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Sleep",interaction="Yes")
lmm7_gi
lmm7_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y + time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Girls",x="Sleep",interaction="No")
lmm7_g
lmm7_bi = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y + sleep*time_rt_avg_12y+time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Sleep",interaction="Yes")
lmm7_bi
lmm7_b = data.frame(summary(est <- pool(with(data=b_da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y + time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="Boys",x="Sleep",interaction="No")
lmm7_b
lmm7i = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y + sleep*time_rt_avg_12y+time_rt_avg_12y + BMI_7y + ses_index_new_10pct + c_age_days_comp_d_qu_12y +as.factor(sex)+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Sleep",interaction="Yes")
lmm7i
#-- No interaction by sex (p-value=0.8998)
lmm7i = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(feeling_score ~ sleep*as.factor(sex)+tannerstg_12y+time_rt_avg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Sleep",interaction="Yes")
lmm7i
lmm7 = data.frame(summary(est <- pool(with(data=da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+as.factor(sex)+(1 | aid), na.action=na.omit))))) %>%
  mutate(sex="All",x="Sleep",interaction="No")
lmm7

lme_dep_sr_sleep <- rbind(lmm7_g, lmm7_b, lmm7) %>%
  filter(term=="sleep") %>%
  mutate(
    sig<-case_when(p.value<0.01 ~ "**",
                   p.value<0.05 ~ "*",
                   T~""),
    Beta_CI = paste0(round(estimate,2),"(",round((estimate-1.96*std.error),2),",",round((estimate+1.96*std.error),2),")",sig)
  ) %>%
  dplyr::select(sex,x,term,Beta_CI)
lme_dep_sr_sleep

write.csv(lme_dep_sr_sleep, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/Longitudinal_lme_dep_sr_sleep.csv")

# 
# #-- Mediation?
# 
# longg <- complete(g_da_long, action="long", include=TRUE)
# 
# ests <- list()
# for(i in seq(1:20)){
#   est <- lmer(sleep ~ tannerstg_12y+BMI_7y+time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 + tannerstg_12y| aid), data=complete(g_da_long,i))
#   ests[[i]] <- ranef(est)$aid$tannerstg_12y
# }
# g<-complete(g_da_long) %>% filter(visit==14)
# ests_d <- cbind(do.call(cbind.data.frame, ests),g$aid)
# names(ests_d) <- c("pred_sl_1","pred_sl_2","pred_sl_3","pred_sl_4","pred_sl_5",
#                    "pred_sl_6","pred_sl_7","pred_sl_8","pred_sl_9","pred_sl_10",
#                    "pred_sl_11","pred_sl_12","pred_sl_13","pred_sl_14","pred_sl_15",
#                    "pred_sl_16","pred_sl_17","pred_sl_18","pred_sl_19","pred_sl_20","aid")
# ests_d <-ests_d %>%
#   pivot_longer(
#     cols = pred_sl_1:pred_sl_20,
#     names_to = c("imp"),
#     values_to = "pred_sl") %>%
#   mutate(
#     .imp = as.numeric(sub(".*pred_sl_", "", imp))
#   ) %>%
#   dplyr::select(aid,pred_sl,.imp)
# str(ests_d)
# str(longg)
# longg1 <- merge(longg,ests_d,by=c(".imp","aid"), all.x = TRUE)
# g_da_long <- as.mids(longg1)
# 
# lmm7_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ sleep+tannerstg_12y+BMI_7y+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit))))) %>%
#   mutate(sex="Girls",x="Sleep",interaction="No")
# lmm7_g
# 
# lmm8_g = data.frame(summary(est <- pool(with(data=g_da_long, exp=lmer(feeling_score ~ pred_sl*time_rt_avg_12y + (1 | aid), na.action=na.omit))))) %>%
#   mutate(sex="Girls",x="Sleep",interaction="No")
# lmm8_g
# 
# g<-complete(g_da_long) %>% filter(visit==14)
# cor(g$sleep,g$pred_sl)
# g<-complete(g_da_long) %>% filter(visit==15)
# cor(g$sleep,g$pred_sl)
# g<-complete(g_da_long) %>% filter(visit==16)
# cor(g$sleep,g$pred_sl)
# 
# 
# #-- Stronger associations between contemporaneous sleep & dep at 14,15y, weaker at 16y
# #-- No indication of non-linearity
# dg14<-complete(g_da_long) %>% filter(visit==14)
# plot(dg14$sleep,dg14$feeling_score)
# cor(dg14$sleep,dg14$feeling_score)
# dg15<-complete(g_da_long) %>% filter(visit==15)
# plot(dg15$sleep,dg15$feeling_score)
# cor(dg15$sleep,dg15$feeling_score)
# dg16<-complete(g_da_long) %>% filter(visit==16)
# plot(dg16$sleep,dg16$feeling_score)
# cor(dg16$sleep,dg16$feeling_score)
# 
# dg14<-complete(b_da_long) %>% filter(visit==14)
# plot(dg14$sleep,dg14$feeling_score)
# cor(dg14$sleep,dg14$feeling_score)
# dg15<-complete(b_da_long) %>% filter(visit==15)
# plot(dg15$sleep,dg15$feeling_score)
# cor(dg15$sleep,dg15$feeling_score)
# dg16<-complete(b_da_long) %>% filter(visit==16)
# plot(dg16$sleep,dg16$feeling_score)
# cor(dg16$sleep,dg16$feeling_score)
# 
# 
# 
# dg15<-complete(b_da_long) %>% filter(visit==15)
# plot(dg15$avgsleeptime_wd,dg15$C5_depr_tscore)
# cor(dg15$avgsleeptime_wd,dg15$C5_depr_tscore)
# 
# summary(lm(data=dg15, avgefficiency_wd ~ tannerstg_12y+time_rt_avg_12y+BMI_7y))
# 
