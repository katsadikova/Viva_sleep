#-----------------------------------------------------------------------------------------------------------#
#--- Figures for slides - LMEs.R
#--- Date:    8/11/2023 (original)
#--- Description: Plot linear mixed effects models for depression outcomes across early-mid adolescence

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
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions


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

rep_cross <- read.csv("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/Viva/results/Table_S1.csv") %>%
  mutate(
    Est1 = case_when(PT=="APHV" ~ -1*Est,
                     T~Est)
  )

#--- Adding SES increases AIC values (poorer model fit) across all three pubertal timing measures (i.e. will only widen CIs, but not add to prediction accuracy)

#-- Model 3: APHV, age, BMI + SES (main effects and interaction with time not significant)
d<-data.frame(complete(g_da_long))
lmm3a <- lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  #--recode to younger age at PHV
  mutate(
    estimate=estimate*-1,
    conf.low1=conf.high*-1,
    conf.high1=conf.low*-1
  ) %>%
  filter(term=="age_peak_velocity_years") %>%
  dplyr::select(term, estimate, conf.low1, conf.high1)
#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0,0,0,0), ncol = 2))
names(d) <- c("visit","int")
dtg <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time
  )
aphv_rep <- rep_cross %>% filter(PT=="APHV" & sex=="Female")
dtg1 <- merge(dtg,aphv_rep,by="visit")

d<-data.frame(complete(b_da_long))
lmm3a <- lmer(feeling_score ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  #--recode to younger age at PHV
  mutate(
    estimate=estimate*-1,
    conf.low1=conf.high*-1,
    conf.high1=conf.low*-1
  ) %>%
  filter(term=="age_peak_velocity_years") %>%
  dplyr::select(term, estimate, conf.low1, conf.high1)
#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0,0,0,0), ncol = 2))
names(d) <- c("visit","int")
dtb <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time
  )
aphv_rep <- rep_cross %>% filter(PT=="APHV" & sex=="Male")
dtb1 <- merge(dtb,aphv_rep,by="visit")

a_plot <- ggplot() + 
  geom_point(data=dtg1, aes(visit, Est1, color="Female")) + 
  geom_line(data=dtg1, aes(x=visit, y=est,  color="Female")) +
  geom_ribbon(data= dtg1, aes(x=visit, ymin=conf.low1, ymax=conf.high1, color="Female"), alpha= 0.3, fill="red") + 
  geom_point(data=dtb1, aes(visit, Est1, color="Male")) + 
  geom_line(data=dtb1, aes(x=visit, y=est, color="Male")) +
  geom_ribbon(data= dtb1, aes(x=visit, ymin=conf.low1, ymax=conf.high1, color="Male"), alpha= 0.3, fill="blue") +
  labs(x="Study visit", y="Beta coefficient") +
  scale_y_continuous(limits=c(-0.4,1.4),
                     breaks = c(-0.4,0,0.4,0.8,1.2)) +
  scale_color_manual(name='Biological sex',
                     breaks=c('Female', 'Male'),
                     values=c('Female'='red', 'Male'='blue')) + 
  ggtitle("Younger APHV") +
  theme_light() +
  theme(text=element_text(size=12,  family="serif"))
a_plot



#-- Model 3: Tanner, age, BMI, SES (interaction with time not significant)
d<-data.frame(complete(g_da_long))
lmm3a <- lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  filter(term=="tannerstg_12y") %>%
  dplyr::select(term, estimate, conf.low, conf.high)
#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0,0,0,0), ncol = 2))
names(d) <- c("visit","int")
dtg <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time
  )
tanner_rep <- rep_cross %>% filter(PT=="Adrenarche" & sex=="Female")
dtg1 <- merge(dtg,tanner_rep,by="visit")

d<-data.frame(complete(b_da_long,2))
lmm3a <- lmer(feeling_score ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  filter(term=="tannerstg_12y") %>%
  dplyr::select(term, estimate, conf.low, conf.high)
#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0,0,0,0), ncol = 2))
names(d) <- c("visit","int")
dtb <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time
  )
tanner_rep <- rep_cross %>% filter(PT=="Adrenarche" & sex=="Male")
dtb1 <- merge(dtb,tanner_rep,by="visit")

t_plot <- ggplot() + 
  geom_point(data=dtg1, aes(visit, Est1, color="Female")) + 
  geom_line(data=dtg1, aes(x=visit, y=est,  color="Female")) +
  geom_ribbon(data= dtg1, aes(x=visit, ymin=conf.low, ymax=conf.high, color="Female"), alpha= 0.3, fill="red") + 
  geom_point(data=dtb1, aes(visit, Est1, color="Male")) + 
  geom_line(data=dtb1, aes(x=visit, y=est, color="Male")) +
  geom_ribbon(data= dtb1, aes(x=visit, ymin=conf.low, ymax=conf.high, color="Male"), alpha= 0.3, fill="blue") +
  labs(x="Study visit", y="Beta coefficient") +
  scale_y_continuous(limits=c(-0.4,1.4),
                     breaks = c(-0.4,0,0.4,0.8,1.2)) +
  scale_color_manual(name='Biological sex',
                     breaks=c('Female', 'Male'),
                     values=c('Female'='red', 'Male'='blue')) + 
  ggtitle("Self-reported Adrenarche") +
  theme_light() +
  theme(text=element_text(size=12,  family="serif"))
t_plot


#-- Model 3: PDS, age, BMI, SES (interaction with time significant for boys)
d<-data.frame(complete(g_da_long))
lmm3a <- lmer(feeling_score ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y + c_age_days_comp_d_qu_12y+(1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  filter(term=="puberty_12y") %>%
  dplyr::select(term, estimate, conf.low, conf.high)
#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0,0,0,0), ncol = 2))
names(d) <- c("visit","int")
dtg <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time
  )
pds_rep <- rep_cross %>% filter(PT=="PDS score" & sex=="Female")
dtg1 <- merge(dtg,pds_rep,by="visit")

d<-data.frame(complete(b_da_long))
lmm3a <- lmer(feeling_score ~ puberty_12y+BMI_7y+ses_index_new_10pct+time_rt_avg_12y+puberty_12y*time_rt_avg_12y +c_age_days_comp_d_qu_12y+ (1 | aid), na.action=na.omit, data=d)
summary(lmm3a)
d$fit = predict(lmm3a)

#-- Pull coefficient & 95%CI for aphv
t <- data.frame(get_model_data(lmm3a, type="est")) %>%
  filter(term=="puberty_12y") %>%
  dplyr::select(term,estimate,std.error)

#-- Create a scaffolding data set with study visit ages and interaction, which is 0 for aphv
d <- data.frame(matrix(c(13,14,15,16,0.162,0.162,0.162,0.162), ncol = 2))
names(d) <- c("visit","int")
dtb <- cbind(d,t) %>%
  mutate(
    time=visit-13,
    est=estimate+int*time,
    conf.low=est-1.96*std.error,
    conf.high=est+1.96*std.error
  )
pds_rep <- rep_cross %>% filter(PT=="PDS score" & sex=="Male")
dtb1 <- merge(dtb,pds_rep,by="visit")



p_plot <- ggplot() + 
  geom_point(data=dtg1, aes(visit, Est1, color="Female")) + 
  geom_line(data=dtg1, aes(x=visit, y=est,  color="Female")) +
  geom_ribbon(data= dtg1, aes(x=visit, ymin=conf.low, ymax=conf.high, color="Female"), alpha= 0.3,fill="red") + 
  geom_point(data=dtb1, aes(visit, Est1, color="Male")) + 
  geom_line(data=dtb1, aes(x=visit, y=est, color="Male")) +
  geom_ribbon(data= dtb1, aes(x=visit, ymin=conf.low, ymax=conf.high, color="Male"), alpha= 0.3,fill="blue") +
  labs(x="Study visit", y="Beta coefficient") +
  scale_y_continuous(limits=c(-0.4,1.4),
                     breaks = c(-0.4,0,0.4,0.8,1.2)) +
  scale_color_manual(name='Biological sex',
                     breaks=c('Female', 'Male'),
                     values=c('Female'='red', 'Male'='blue')) + 
  ggtitle("Parent-reported PDS Score") +
  theme_light() +
  theme(text=element_text(size=12,  family="serif"))
p_plot

ggarrange(a_plot,t_plot,p_plot,
          ncol=3, 
          nrow=1, common.legend = TRUE,
          legend="right")
