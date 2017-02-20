## LOAD PACKAGES ####
library(mgcv) # For GAMs
library(lme4)     # For LMERs


## READ IN DATA ####
data = read.table("data_gamm.txt", header = T, sep = "\t")


## LINEAR REGRESSION ####
# Run model
percip.lm = lm(Percipitation_Median ~ Month_Num_Real, data = data)

# Look at results
percip.lm_sum = summary(percip.lm)
percip.lm_sum


## GENERAL ADDITIVE MODEL ####
# Run model
percip.gam.10 = gam(Percipitation_Median ~ s(Month_Num_Real, k = 10, bs="tp"),
                    data = data)

# Look at results
percip.gam.10_sum = summary(percip.gam.10)
percip.gam.10_sum


## LINEAR MIXED EFFECTS MODEL ####
# Run model with random intercept only
percip.lmer.int = lmer(Percipitation_Median ~ Month_Num_Real +
                         (1|River_Basin),
                       data = data)

# Look at results
percip.lmer.int_sum = summary(percip.lmer.int)
percip.lmer.int_sum

# Run model with random intercept and random slope
percip.lmer.intslp = lmer(Percipitation_Median ~ Month_Num_Real +
                         (1+Month_Num_Real|River_Basin),
                       data = data)

# Look at results
percip.lmer.intslp_sum = summary(percip.lmer.intslp)
percip.lmer.intslp_sum


## GENERAL ADDITIVE MIXED EFFECTS MODEL ####
# Run model with random intercept only
percip.gamm.10.int = gam(Percipitation_Median ~ s(Month_Num_Real, k = 10, bs = "tp") +
                           s(River_Basin, bs = "re"),
                         data = data)

# Look at results
percip.gamm.10.int_sum = summary(percip.gamm.10.int)
percip.gamm.10.int_sum

# Run model with random intercept and random slope
percip.gamm.10.intslp = gam(Percipitation_Median ~ s(Month_Num_Real, k = 10, bs = "tp") +
                              s(River_Basin, bs = "re") +
                              s(River_Basin, Month_Num_Real, bs = "re"),
                            data = data)

# Look at results
percip.gamm.10.intslp_sum = summary(percip.gamm.10.intslp)
percip.gamm.10.intslp_sum

# Run model with random smooth
# NOTE: gives a warning, may need to add more parameters to random smooth
percip.gamm.10.smt = gam(Percipitation_Median ~ s(Month_Num_Real, k = 10, bs = "tp") +
                           s(Month_Num_Real, River_Basin, bs = "fs", m = 1),
                         data = data)

# Look at results
percip.gamm.10.int_smt = summary(percip.gamm.10.smt)
percip.gamm.10.int_smt

