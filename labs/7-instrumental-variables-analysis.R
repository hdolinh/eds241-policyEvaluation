library(stargazer) # tables
library(estimatr) # lm_robust()
library(AER) 
library(ggplot2)
library(cowplot)
library(sandwich)
library(lmtest)
library(dplyr)
library(lfe) # felm()
library(here)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set's directory where script is located
# getwd()

# IMPORT CSV DATA
FULTON <- read.csv(here("data/FULTON.csv")) %>%
  mutate(log_tots = log(tots),
         log_price = log(pricelevel))


# SUMMARY STATISTICS
stargazer(FULTON, type="text", digits=2)


# BASIC OLS REGRESSION
ols <- lm(formula = log_tots ~ log_price, data=FULTON)
summary(ols)


# FIRST_STAGE REGRESSION - JUST-IDENTIFIED MODEL
fs1 <- lm(formula = log_price ~ windspd, data=FULTON) # what is the effect of windspeed on log price of fish
summary(fs1) # found that windspeed increases price by 0.7181
# reject null hypothesis with small p-value
# high wind could make it difficult to fish


# TSLS - JUST-IDENTIFIED MODEL
tsls1 <- ivreg(log_tots ~ log_price | windspd, data = FULTON) # vertical bar specifies list of instrument 
summary(tsls1)


# Calculate robust standard errors for OLS and FS1 using starprep()
se_ols_fs1 <- starprep(ols,fs1, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

# Calculate robust standard errors sandwich and lmtest packages (starprep() does not like ivreg() objects)
se_tsls11 <- coeftest(tsls1, vcov = vcovHC(tsls1, type = "HC2"))[, "Std. Error"]

# Combine standard errors and output results with stargazer()
se_models <- append(se_ols_fs1,list(se_tsls11))
stargazer(ols, fs1, tsls1, se = se_models, type="text")



# Other approach using the lfe package #####################

# Estimate the first two models
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

