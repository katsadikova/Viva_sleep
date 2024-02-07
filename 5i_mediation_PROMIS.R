#----------------------------------------------------------------------#
#--- 5i_mediation_PROMIS.R
#--- Date: 7/3/2023
#--- Date updated: 2/7/2024
#----------------------------------------------------------------------#

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
long$sr_sl16_8hr <- ifelse(long$sleep_weekday_cq16<7,1,0)
long$promis_modsev <- ifelse(long$C5_depr_tscore>=65,1,0)

#--- Convert back to mids object
d_imp <- as.mids(long)
#----------------------------------------------------------------------------------------#

#-- Subset the imputed mids object by sex
d_imp_b <- filter(d_imp, sex == "Male")
nrow(complete(d_imp_b)) #585
d_imp_g <- filter(d_imp, sex == "Female")
nrow(complete(d_imp_g)) #553



#-- Mediation by longitudinally self-reported weeknight sleep IS significant (Tanner --- Promis Dep).
ests_17 <- list()
for(i in seq(1:20)){
  print(i)
  d_imp_i <- complete(d_imp,i)
  set.seed(23)
  mediation.rb <- cmest(data = d_imp_i,
                        model = "gformula",
                        outcome = "C5_depr_tscore",
                        exposure = "tannerstg_12y",
                        mediator = c("sleep_schoolday_cq14","sleep_weekday_child_cq15","sleep_weekday_cq16"),
                        EMint = F,
                        basec = c("BMI_7y","ses_index_new_10pct","c_age_days_comp_d_qu_12y","sex"),
                        mreg = list("linear","linear","linear"),
                        yreg = "linear",
                        a = 2,
                        astar = 1,
                        mval = list(8,8,8),
                        estimation = "imputation",
                        inference = "bootstrap")
  ests_17[[i]]<-tibble::rownames_to_column(data.frame(summary(mediation.rb)$summarydf),"Effect")
}
ests_17_d <- do.call(rbind.data.frame, ests_17) 
names(ests_17_d) <- c("Effect","est","se","low","high","pval")
ests_17_mean <-ests_17_d %>%
  group_by(Effect) %>%
  summarise(n_imp=n(),
            est_mean=mean(est),
            lowCI_mean=mean(low),
            highCI_mean=mean(high))
ests_17_mean



#----------------------------------------------------------------------------------------------------------------------------------------#
#---- Testing other possibilities - not real chance of them panning out due to individual relationships not really being substantiated --#
#-- Mediation by actigraphy-measured ET sleep is not significant
ests_17a <- list()
for(i in seq(1:20)){
  print(i)
  d_imp_i <- complete(d_imp,i)
  set.seed(23)
  mediation.rb <- cmest(data = d_imp_i,
                        model = "gformula",
                        outcome = "C5_depr_tscore",
                        exposure = "tannerstg_12y",
                        mediator = c("avgsleeptime_wd_hr"),
                        EMint = F,
                        basec = c("BMI_7y","ses_index_new_10pct","c_age_days_comp_d_qu_12y","sex"),
                        mreg = list("linear"),
                        yreg = "linear",
                        a = 2,
                        astar = 1,
                        mval = list(8),
                        estimation = "imputation",
                        inference = "bootstrap")
  ests_17a[[i]]<-tibble::rownames_to_column(data.frame(summary(mediation.rb)$summarydf),"Effect")
}
ests_17a_d <- do.call(rbind.data.frame, ests_17a) 
names(ests_17a_d) <- c("Effect","est","se","low","high","pval")

ests_17a_mean <-ests_17a_d %>%
  group_by(Effect) %>%
  summarise(n_imp=n(),
            est_mean=mean(est),
            lowCI_mean=mean(low),
            highCI_mean=mean(high))
ests_17a_mean


#-- Mediation adjusted for post-baseline confounder of med-outcome relationship --> feeling_score_ET.
ests_17 <- list()
for(i in seq(1:20)){
  print(i)
  d_imp_i <- complete(d_imp,i)
  set.seed(23)
  mediation.rb <- cmest(data = d_imp_i,
                        model = "gformula",
                        outcome = "C5_depr_tscore",
                        exposure = "tannerstg_12y",
                        mediator = c("sleep_schoolday_cq14","sleep_weekday_child_cq15","sleep_weekday_cq16"),
                        EMint = F,
                        basec = c("BMI_7y","ses_index_new_10pct","c_age_days_comp_d_qu_12y","sex"),
                        postc = c("feeling_score_ET"),
                        postcreg = list("linear"),
                        mreg = list("linear","linear","linear"),
                        yreg = "linear",
                        a = 2,
                        astar = 1,
                        mval = list(8,8,8),
                        estimation = "imputation",
                        inference = "bootstrap")
  ests_17[[i]]<-tibble::rownames_to_column(data.frame(summary(mediation.rb)$summarydf),"Effect")
}
ests_17_d <- do.call(rbind.data.frame, ests_17) 
names(ests_17_d) <- c("Effect","est","se","low","high","pval")

ests_17_mean <-ests_17_d %>%
  group_by(Effect) %>%
  summarise(n_imp=n(),
            est_mean=mean(est),
            lowCI_mean=mean(low),
            highCI_mean=mean(high))
ests_17_mean


#-- Mediation by longitudinally self-reported weeknight sleep is not significant (PDS ---> Promis Dep).
ests_17p <- list()
for(i in seq(1:20)){
  print(i)
  d_imp_i <- complete(d_imp,i)
  set.seed(23)
  mediation.rb <- cmest(data = d_imp_i,
                        model = "gformula",
                        outcome = "C5_depr_tscore",
                        exposure = "puberty_12y",
                        mediator = c("sleep_schoolday_cq14","sleep_weekday_child_cq15","sleep_weekday_cq16"),
                        EMint = F,
                        basec = c("BMI_7y","ses_index_new_10pct","c_age_days_comp_d_qu_12y","sex"),
                        mreg = list("linear","linear","linear"),
                        yreg = "linear",
                        a = 2,
                        astar = 1,
                        mval = list(8,8,8),
                        estimation = "imputation",
                        inference = "bootstrap")
  ests_17p[[i]]<-tibble::rownames_to_column(data.frame(summary(mediation.rb)$summarydf),"Effect")
}
ests_17p_d <- do.call(rbind.data.frame, ests_17p) 
names(ests_17p_d) <- c("Effect","est","se","low","high","pval")

ests_17p_mean <-ests_17p_d %>%
  group_by(Effect) %>%
  summarise(n_imp=n(),
            est_mean=mean(est),
            lowCI_mean=mean(low),
            highCI_mean=mean(high))
ests_17p_mean
