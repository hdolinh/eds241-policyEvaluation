library(MASS)
library(ggplot2)
library(vtable)
library(stargazer)
library(estimatr)
library(dplyr)
library(tidyr)
library(here)

### Directory

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set's directory where script is located

# getwd()


# IMPORT CSV DATA
CGL <- read.csv(here::here("data/cgl_collapse_data_extract.csv"))


# SUMMARY STATISTICS
# always first check the summary statistics
# DAPever this is our treatment for a 120 fisheries 
# LME large mor... ecosystem
stargazer(CGL, type="text", digits=2)


# EXAMINE BALANCE IN COVARIATES
# COVARIATE MEAN DIFFERENCES by DAPever
m1 <- lm(formula = LME ~ DAPever, data=CGL)
m2 <- lm(formula = genus ~ DAPever, data=CGL)
m3 <- lm(formula = species ~ DAPever, data=CGL)
se_models = starprep(m1, m2, m3, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(m1, m2, m3, se = se_models, type="text")
# use stargazer to compare the 3 models
# we should expect to see no difference in the average of the covariants bc this has been randomized
# LME has a larger average by 5.406

# BOXPLOTS TO EXAMINE BALANCE IN COVARIATES
# one box for values where treatment is 1
# one box for values where treatment is 0
ggplot(CGL, aes(x=as.factor(DAPever), y=LME)) + 
  geom_boxplot(fill="cyan") + xlab("ITQ Yes/No")
# LME is really a categorical variable
# median value is around 35 which is similar to the ITQ's are treated
# the overlap is a bit problematic bc if we want to compare the means's...bc we can't compare the 75th percentile and there's a large difference between the 25th percentiles

ggplot(CGL, aes(x=as.factor(DAPever), y=genus)) + 
  geom_boxplot(fill="cyan") + xlab("ITQ Yes/No")
# this looks pretty good!

ggplot(CGL, aes(x=as.factor(DAPever), y=species)) + 
  geom_boxplot(fill="cyan") + xlab("ITQ Yes/No")
# species also looks good


# BASIC OLS by DAPever -- THEN ADD INDICATORS FOR OTHER COVARIATES 
# NOTE DO NOT INCLUDE SPECIES IN MODELS TO KEEP RUNNING TIME FAST
mA <- lm(formula = collapse ~ DAPever, data=CGL) # collapse on DAPever
mB <- lm(formula = collapse ~ DAPever + as.factor(LME), data=CGL) # add a dummy
mC <- lm(formula = collapse ~ DAPever + as.factor(LME) + as.factor(genus), data=CGL) # another dummy
se_models = starprep(mA, mB, mC, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(mA, mB, mC, se = se_models, type="text", omit = "(LME)|(genus)|(species)")
# first assumption of our treatment


# BASIC PROPENSITY SCORE --- THIS IS A TOY MODEL
# ESTIMATE PROPENSITY SCORE MODEL AND PREDICT (EPS)
ps_model <- glm(DAPever ~ LME + genus, family = binomial(), data = CGL)
summary(ps_model)
EPS <- predict(ps_model, type = "response")
PS_WGT <- (CGL$DAPever/EPS) + ((1-CGL$DAPever)/(1-EPS))


# COLLECT ALL RELEVANT VARIABLES IN DATAFRAME
DF <- data.frame(years = CGL$years, collapse = CGL$collapse, DAPever = CGL$DAPever, 
                 LME = CGL$LME, genus = CGL$genus, species = CGL$species, EPS, PS_WGT)


# BOXPLOTS TO EXAMINE OVERLAP IN P-SCORE DISTRIBUTIONS
ggplot(DF, aes(x=as.factor(DAPever), y=EPS)) + 
  geom_boxplot(fill="cyan") + xlab("Ever Collapsed")
# pretty good overlap


# WLS USING EPS WEIGHTS
wls1 <- lm(formula = collapse ~ DAPever, data=DF, weights=PS_WGT)
wls2 <- lm(formula = collapse ~ DAPever + LME + genus, data=DF, weights=PS_WGT)
se_models = starprep(wls1, wls2, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(wls1, wls2, se = se_models, type="text", omit = "(LME)|(genus)|(species)")
# if you don't do any adjustments we got -0.142
# if we did reweighing, etc we got -0.148
# this tells us there is no bias or there it didn't remove / control bias
