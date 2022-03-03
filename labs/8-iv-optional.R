library(stargazer)
library(estimatr)
library(AER)
library(ggplot2)
library(cowplot)
library(sandwich)
library(lmtest)
library(dplyr)
library(lfe)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set's directory where script is located
getwd()

# IMPORT CSV DATA
FULTON <- read.csv("FULTON.csv")%>%
  mutate(log_tots = log(tots),
         log_price = log(pricelevel))


# SUMMARY STATISTICS
stargazer(FULTON, type="text", digits=2)


# Other approach to TSLS using the lfe package (OPTIONAL) #####################
# Focus on just-identified model only

# Estimate the OLS equation first-stage regression
ols_felm <- felm(formula = log_tots ~ log_price, data=FULTON)
fs1_felm <- felm(formula = log_price ~ windspd, data=FULTON)

# Estimate 2SLS
# "log_tots ~ 1" is not the first stage, it is the all variables in the first stage, BUT the endogenous one
# | 0 | means that we are not including fixed effects here.

tsls1_felm <- felm(formula = log_tots ~ 1 | 0 | (log_price ~ windspd),  data=FULTON)

# The robust standard errors are calculated (not reported) by default in felm(), so here we can fetch and combine them
# It might be HC1, but the documentation is not great. 

se_models_felm <- list(ols_felm$rse,fs1_felm$rse, tsls1_felm$rse)

stargazer(ols_felm, fs1_felm, tsls1_felm, se = se_models_felm, type="text")

