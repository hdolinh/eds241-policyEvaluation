library(stargazer)
library(estimatr)
library(ggplot2)
library(plm)
library(sandwich)
library(lmtest)
library(dplyr)
library(tidyr)
library(here)
library(tidyverse)

# SET WORKING DIRECTORY
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Automatically sets directory where script is saved
# getwd()

# IMPORT CSV DATA 
RPS <- read_csv(here::here("data/RPS_data.csv"))

# SUMMARY STATISTICS
stargazer(RPS, type="text", digits=2)

# LIST SOME VARIABLES FOR CALIFORNIA
RPS%>%
  filter(state_name == "California")%>%
  select(state_name, year, rps_D, rps_ever_adopter, rps_implementation_year)%>%
  View

# DD REGRESSION, Y = Wind+Solar installed capacity (MW), using lm package
DD_cap1 <- lm(formula = cap_WS_mw ~ rps_D + as.factor(state_name) + as.factor(year), data=RPS)
se_DD_cap1 <- starprep(DD_cap1, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

DD_cap2 <- lm(formula = cap_WS_mw ~ rps_D + as.factor(state_name) + as.factor(year), data=RPS)
se_DD_cap2 <- starprep(DD_cap2, stat = c("std.error"), se_type = "CR2", clusters=RPS$state_name, alpha = 0.05) 

se_models <- list(se_DD_cap1[[1]], se_DD_cap2[[1]])
stargazer(DD_cap1, DD_cap2, se = se_models, keep=c("rps_D"), type="text")


# DD REGRESSION, Y = Wind+Solar generation (GWh), using plm package
DD_gen1 <- plm(gen_WS_gwh ~ rps_D, 
               index = c("state_name", "year"), model = "within", effect = "twoways", data = RPS)

# Calculate standard errors (note slightly different procedure with plm package)
se_DD_gen1 <- coeftest(DD_gen1, vcov = vcovHC(DD_gen1, type = "HC2"))[, "Std. Error"]
# Reformat standard errors for stargazer()
se_DD_gen1 <- list(se_DD_gen1)
# Output results with stargazer
stargazer(DD_gen1, keep=c("rps_D"), se = se_DD_gen1, type="text")

