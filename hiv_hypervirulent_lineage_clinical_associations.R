# Author: Chris Wymant
# Acknowledgment: I wrote this while funded by ERC Advanced Grant PBDR-339251
# and a Li Ka Shing Foundation grant, both awarded to Christophe Fraser.
# Abbreviations: VL = viral load (log10 copies per ml), LMM = linear mixed model,
# df = dataframe

library(tidyverse)
library(lme4)

# WARNING: deletes all objects in your R session's memory
rm(list = ls()) 

setwd("~/Dropbox (Infectious Disease)/BEEHIVE/LineageEffect/AssociationTesting_PaperVersion/")

# Read in the input csvs. For you, lucky reader, this is after much data
# wrangling.
# "individual" characteristics (one row per individual):
df_ind <- read_csv("individual_summary.csv", col_types = cols(
  id_paper = col_factor(),
  in_lineage = col_factor(),
  diagnosis_period = col_factor(),
  sex = col_factor(),
  age_diagnosed = col_factor(),
  death_at_censoring = col_logical()
))
# And CD4 counts (many per individual):
df_cd4_decline <- read_csv("cd4_counts_pretreatment.csv", col_types = cols(
  id_paper = col_factor()
))

# Inspect the mean VL in each time period for lineage and not lineage (what we
# plot in Figure 1a).
df_ind %>%
  group_by(diagnosis_period, in_lineage) %>%
  summarise(vl_group_mean = mean(vl_mean, na.rm = TRUE),
            number_of_vl_measurements = sum(!is.na(vl_mean)))

# Merge individual-level data into the CD4 counts.
stopifnot(all(df_cd4_decline$id_paper %in% df_ind$id_paper))
stopifnot(! anyNA(df_ind$in_lineage))
df_cd4_decline <- left_join(df_cd4_decline, df_ind, by = "id_paper") 

# Calculate the mean VL of the reference category: not-lineage males diagnosed 
# in their thirties. We'll measure VLs offset by this value, i.e. relative to 
# it, for the VL-adjusted LMM, for more easily interpretable regression
# coefficients. 
mean_vl_reference_category <- df_ind %>%
  filter(sex %in% "male",
         age_diagnosed %in% "[30, 40)",
         in_lineage %in% "not.lineage") %>%
  pull(vl_mean) %>%
  mean(na.rm = TRUE)
df_cd4_decline <- df_cd4_decline %>%
  mutate(vl_ref_cat_diff = vl_mean - mean_vl_reference_category)

# Define the reference category
df_cd4_decline <- df_cd4_decline %>%
  mutate(sex = relevel(sex, ref = "male"),
         age_diagnosed = relevel(age_diagnosed, ref = "[30, 40)"),
         in_lineage = relevel(df_cd4_decline$in_lineage, ref = "not.lineage"))

lmm <- lmer(data = df_cd4_decline, 
            # Model CD4 counts as a linear function of time,
            cd4_count ~ years_since_diagnosis +  
              # with a fixed effect of age on the intercept,
              age_diagnosed + 
              # a fixed effect of sex on both intercept and slope,
              years_since_diagnosis * sex + 
              # a fixed effect of the lineage on both intercept and slope,
              years_since_diagnosis * in_lineage + 
              # and a random effect of the individual on both intercept and slope.
              (years_since_diagnosis | id_paper))  
summary(lmm)

lmm_adjusted_for_vl <- lmer(data = df_cd4_decline, 
                            # Model CD4 counts as a linear function of time,
                            cd4_count ~ years_since_diagnosis +  
                              # with a fixed effect of age on the intercept,
                              age_diagnosed + 
                              # a fixed effect of sex on both intercept and slope,
                              years_since_diagnosis * sex + 
                              # a fixed effect of the lineage on both intercept and slope,
                              years_since_diagnosis * in_lineage + 
                              # a fixed effect of VL on both intercept and slope,
                              years_since_diagnosis * vl_ref_cat_diff + 
                              # and a random effect of the individual on both intercept and slope.
                              (years_since_diagnosis | id_paper))  
summary(lmm_adjusted_for_vl)

# Calculate LMM confidence intervals. Takes a few minutes each.
confints <- confint(lmm)
confints_adjusted_for_vl <- confint(lmm_adjusted_for_vl)
confints
confints_adjusted_for_vl

# Cox model for the survival analysis.
# Test association between the hazard for death and age, sex and the lineage.
fit.coxph <- coxph(data = df_survival,
                   Surv(time = years_to_censoring, event = death_at_censoring) ~
                      age_diagnosed + sex + in_lineage)
summary(fit.coxph)
