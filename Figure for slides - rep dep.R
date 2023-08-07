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
#-- BMI and SES adjusted relationships
#-- Tanner & Dep
g_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+race_black+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="11Y") 
g_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="12Y") 
g_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="14Y") 
g_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="15Y") 
g_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="16Y")  
g_dep17adj_tanner <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(C5_depr_tscore ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="Adrenarche", y="17Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Adrenarche", y="11Y")
b_dep12adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Adrenarche", y="12Y") 
b_dep14adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Adrenarche", y="14Y") 
b_dep15adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Adrenarche", y="15Y") 
b_dep16adj_tanner <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ tannerstg_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="Adrenarche", y="16Y") 

#-- Put all the coefficients for tanner together, format sub-table
tanner_dep_adj <- rbind(g_dep11adj_tanner,g_dep12adj_tanner,g_dep14adj_tanner,g_dep15adj_tanner,g_dep16adj_tanner,
                        b_dep11adj_tanner,b_dep12adj_tanner,b_dep14adj_tanner,b_dep15adj_tanner,b_dep16adj_tanner) %>% 
  filter(term=="tannerstg_12y" & y != "11Y") %>%
  mutate(
    Beta = round(estimate,2),
    SE = round(std.error,2)
  ) %>%
  dplyr::select(sex,x,y,Beta,SE) 


#-- PDS & Dep
g_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS score", y="11Y") 
g_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS score", y="12Y") 
g_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS score", y="14Y")
g_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS score", y="15Y") 
g_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Girls", x="PDS score", y="16Y")

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS score", y="11Y") 
b_dep12adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS score", y="12Y") 
b_dep14adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS score", y="14Y") 
b_dep15adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS score", y="15Y") 
b_dep16adj_pds <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ puberty_12y+BMI_7y+ses_index_new_10pct+c_age_days_comp_d_qu_12y))))) %>%
  mutate(sex="Boys", x="PDS score", y="16Y") 

#-- Put all the coefficients for tanner together, format sub-table
pds_dep_adj <- rbind(g_dep11adj_pds,g_dep12adj_pds,g_dep14adj_pds,g_dep15adj_pds,g_dep16adj_pds,
                     b_dep11adj_pds,b_dep12adj_pds,b_dep14adj_pds,b_dep15adj_pds,b_dep16adj_pds) %>% 
  filter(term=="puberty_12y" & y != "11Y") %>%
  mutate(
    Beta = round(estimate,2),
    SE = round(std.error,2)
  ) %>%
  dplyr::select(sex,x,y,Beta,SE) 


#-- APHV & Dep
g_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Younger APHV", y="11Y")
g_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Younger APHV", y="12Y") 
g_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Younger APHV", y="14Y") 
g_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Younger APHV", y="15Y")
g_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_g, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Girls", x="Younger APHV", y="16Y") 

#-- Look at relationships between dep, Tanner, and actigraphy-measured sleep in boys
b_dep11adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_11y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="Younger APHV", y="11Y") 
b_dep12adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_ET ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="Younger APHV", y="12Y")
b_dep14adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_14y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="Younger APHV", y="14Y") 
b_dep15adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_15y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="Younger APHV", y="15Y") 
b_dep16adj_aphv <- data.frame(summary(pool(with(data = d_imp_b, exp = lm(feeling_score_16y ~ age_peak_velocity_years+BMI_7y+ses_index_new_10pct))))) %>%
  mutate(sex="Boys", x="Younger APHV", y="16Y") 

#-- Put all the coefficients for tanner together, format sub-table
aphv_dep_adj <- rbind(g_dep11adj_aphv,g_dep12adj_aphv,g_dep14adj_aphv,g_dep15adj_aphv,g_dep16adj_aphv,
                      b_dep11adj_aphv,b_dep12adj_aphv,b_dep14adj_aphv,b_dep15adj_aphv,b_dep16adj_aphv) %>% 
  filter(term=="age_peak_velocity_years" & y != "11Y") %>%
  mutate(
    Beta = -1*round(estimate,2),
    SE = round(std.error,2)
  ) %>%
  dplyr::select(sex,x,y,Beta,SE) 



fig1dat <- data.frame(rbind(aphv_dep_adj, tanner_dep_adj, pds_dep_adj)) %>%
  mutate(
    PT=x
  )
fig1dat_g <- fig1dat %>% filter(sex=="Girls")
fig1dat_b <- fig1dat %>% filter(sex=="Boys")

library(ggplot2)
library(gridExtra)
library(ggpubr) # Combine figures
library(extrafont)
loadfonts(device = "all")

dodge <- position_dodge(width=0.5)
fig1_g <- ggplot(fig1dat_g, aes(x=y, y=Beta, colour=PT)) +
  geom_errorbar(aes(ymin=Beta-1.96*SE, ymax=Beta+1.96*SE), width=.1,position=dodge) +
  geom_line(size=2,position=dodge) +
  geom_point(position=dodge) +
  geom_hline(yintercept=0) +
  ylab("Linear regression coefficient, adjusted") +
  xlab("Outcome") +
  scale_y_continuous(limits=c(-0.5,1.7),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,
                                0.4,0.5,0.6,0.7,0.8,0.9,1,
                                1.1,1.2,1.3,1.4,1.5,1.6,1.7)) +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  ggtitle("Females: Associations between PT & repeatedly measured depression symptoms") +
  theme_light() +
  theme(text=element_text(size=12,  family="serif"))
fig1_g

fig1_b <- ggplot(fig1dat_b, aes(x=y, y=Beta, colour=PT)) +
  geom_errorbar(aes(ymin=Beta-1.96*SE, ymax=Beta+1.96*SE), width=.1,position=dodge) +
  geom_line(position=dodge, width=0.5) +
  geom_point(position=dodge) +
  geom_hline(yintercept=0) +
  ylab("Linear regression coefficient, adjusted") +
  xlab("Outcome") +
  scale_y_continuous(limits=c(-0.5,1.5),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,
                                0.4,0.5,0.6,0.7,0.8,0.9,1,
                                1.1,1.2,1.3,1.4,1.5,1.6,1.7)) +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  ggtitle("Males: Associations between PT & repeatedly measured depression symptoms") +
  theme_light() +
  theme(text=element_text(size=12,  family="serif"))
fig1_b

ggarrange(fig1_g,fig1_b,
          ncol=1, 
          nrow=2, common.legend = TRUE,
          legend="right")
