#-----------------------------------------------------------------------------------------------------------#
#--- 4i_gformula_sleep_PROMIS.R
#--- Date: 7/16/2023 (updated 7/17/2023)
#--- Description: Estimate the effect of increasing weeknight self-reported sleep on PROMIS dep symptoms at age 17
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

#--------------------------------------------------------------------------------------------------

ests <- list()

for(i in seq(1:20)){
  print(i)

  d <- complete(da_long,i) %>% mutate(
    visit=visit-14,
    boy=ifelse(sex=="Male",1,0)
  )
  
  library(gfoRmula)
  id <- "aid"
  time_points <- 3
  time_name <- "visit"
  covnames <- c("sleep")
  basecovs <- c("tannerstg_12y","boy","BMI_7y","ses_index_new_10pct","age")
  outcome_name <- "C5_depr_tscore"
  outcome_type <- "continuous_eof"
  histories <- c(lagged)
  histvars <- list(c("sleep"))
  covtypes <- c("bounded normal")
  covparams <- list(covmodels = c(sleep ~ lag1_sleep + tannerstg_12y + boy + BMI_7y + ses_index_new_10pct + age))
  ymodel <- C5_depr_tscore ~ sleep + tannerstg_12y + boy + BMI_7y + age
  intvars <- list("sleep","sleep")
  interventions <- list(list(c(static, rep(6, 3))),
                        list(c(static, rep(7, 3))))
  interventions
  int_descript <- c("sleep=7","sleep=8")
  
  gform_cont_eof <- gformula_continuous_eof(obs_data = d,
                                            id = id,
                                            time_name = time_name,
                                            covnames = covnames,
                                            outcome_name = outcome_name,
                                            covtypes = covtypes,
                                            covparams = covparams, 
                                            ymodel = ymodel,
                                            intvars = intvars,
                                            interventions = interventions,
                                            int_descript = int_descript,
                                            histories = histories, 
                                            histvars = histvars,
                                            basecovs = basecovs,
                                            nsimul = 100000, 
                                            nsamples = 1000,
                                            boot_diag = T,
                                            seed = 1234)
  est <- round(gform_cont_eof$coeffs$C5_depr_tscore[2],2)
  SE <- round(gform_cont_eof$stderrs$C5_depr_tscore[2],2)
  ests[[i]] <- cbind(est,SE)
}

ests_d <- do.call(rbind.data.frame, ests) 
ests_mean <-ests_d %>%
  summarise(n_imp=n(),
            mean=mean(est),
            Vw = mean(SE*SE),
            Vb = sum((est - mean)^2) / 20,
            SE_MEAN = mean(SE),
            SE = sqrt(Vw+Vb+Vb/n_imp),
            p = 2*pnorm(abs(mean / SE), lower.tail=F))
ests_mean
write.csv(ests_d, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/sleep_PROMIS_with_SES_gformula_ests_by_imputation.csv")
write.csv(ests_mean, "/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results_imp_chronological/sleep_PROMIS_with_SES_gformula_MEAN_est.csv")

