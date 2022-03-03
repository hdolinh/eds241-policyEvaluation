library(stargazer)
library(estimatr)
library(ggplot2)
library(plm)
library(sandwich)
library(lmtest)
library(dplyr)
library(tidyr)


# SET WORKING DIRECTORY

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Automatically sets directory where script is saved
getwd()

# IMPORT CSV DATA 
QWI <- read.csv("qwi_employment_wind.csv")%>%
  mutate(pop1000 = population/1000)

# SUMMARY STATISTICS
stargazer(QWI, type="text", digits=2)


# BASIC OLS REGRESSION
ols1 <- lm(formula = log_emp ~ capacity_gwwind + pop1000 + gdp_state_millions_usd, data=QWI)
se_ols1 <- starprep(ols1, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

# OLS REGRESSION, adding year fixed effects
ols2 <- lm(formula = log_emp ~ capacity_gwwind + pop1000 + gdp_state_millions_usd + as.factor(year), data=QWI)
se_ols2 <- starprep(ols2, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

# OLS REGRESSION, adding year fixed effects, and cluster-robust std errors (state)
ols3 <- lm(formula = log_emp ~ capacity_gwwind + pop1000 + gdp_state_millions_usd + as.factor(year), data=QWI)
se_ols3 <- starprep(ols3, stat = c("std.error"), se_type = "CR2", clusters=QWI$state_name, alpha = 0.05) 

se_models <- list(se_ols1[[1]], se_ols2[[1]], se_ols3[[1]])
stargazer(ols1, ols2, ols3, se = se_models, keep=c("capacity_gwwind"), type="text")


# ESTIMATE THE BASIC WITHIN FIXED EFFECTS REGRESSION
within1 <- plm(log_emp ~ capacity_gwwind + pop1000 + gdp_state_millions_usd, 
               index = c("county_id", "year"), model = "within", effect = "twoways", data = QWI)

# Calculate standard errors (note slightly different procedure with plm package)
se_within1 <- coeftest(within1, vcov = vcovHC(within1, type = "HC2"))[, "Std. Error"]
# Reformat standard errors for stargazer()
se_within1 <- list(se_within1)
# Output results with stargazer
stargazer(within1, keep=c("capacity_gwwind"), se = se_within1, type="text")



