library(stargazer)
library(estimatr)
library(AER)
library(ggplot2)
library(cowplot)
library(sandwich)
library(lmtest)
library(dplyr)
library(lfe)
library(here)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set's directory where script is located
# # getwd()

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
# Lecture 8, slide 12
fs1 <- lm(formula = log_price ~ windspd, data=FULTON)
summary(fs1)
# F-test for non-weak and relevant instruments (Lecture 9, slides 13-14)
linearHypothesis(fs1, c("windspd=0"), white.adjust = "hc2")
# don't use the f statistic here because all characteristics are being tested here

# TSLS - JUST-IDENTIFIED MODEL
# Lecture 8, slide 13
tsls1 <- ivreg(log_tots ~ log_price | windspd, data = FULTON)
summary(tsls1)


# Calculate robust standard errors for OLS and FS1 using starprep()
se_ols_fs1 <- starprep(ols,fs1, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

# Calculate robust standard errors for TSLS1 using sandwich and lmtest packages (starprep() does not like ivreg() objects)
se_tsls1 <- coeftest(tsls1, vcov = vcovHC(tsls1, type = "HC2"))[, "Std. Error"]

# Combine standard errors and output results with stargazer()
se_models <- append(se_ols_fs1,list(se_tsls1))
stargazer(ols, fs1, tsls1, se = se_models, type="text")






# FIRST_STAGE REGRESSION - OVER-IDENTIFIED MODEL
# Lecture 9, slide 6
fs2 <- lm(formula = log_price ~ day1 + day2 + day3 + day4 + windspd + cold, data=FULTON)
summary(fs2)

# F-test for non-weak and relevant instruments (Lecture 9, slides 13-14)
linearHypothesis(fs2, c("windspd=0", "cold=0"), white.adjust = "hc2")
# even tho p value is significant we don't move forward because the f statistic is below 10

# TSLS - OVER-IDENTIFIED MODEL
# Lecture 9, Slide 7
tsls2 <- ivreg(log_tots ~ log_price + day1 + day2 + day3 + day4 | day1 + day2 + day3 + day4 + windspd + cold, data = FULTON)
summary(tsls2)
# right side of bar are the instruments 

# FINAL TSLS - REMOVED WEAKED INSTRUMENT
tsls3 <- ivreg(log_tots ~ log_price + day1 + day2 + day3 + day4 | day1 + day2 + day3 + day4 + windspd, data = FULTON)
summary(tsls3)


# Calculate robust standard errors for FS2 using starprep()
se_fs2 <- starprep(fs2, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 

# Calculate robust standard errors for TSLS2 using sandwich and lmtest packages (starprep() does not like ivreg() objects)
se_tsls2 <- coeftest(tsls2, vcov = vcovHC(tsls2, type = "HC2"))[, "Std. Error"]

# Calculate robust standard errors for TSLS3 using sandwich and lmtest packages (starprep() does not like ivreg() objects)
se_tsls3 <- coeftest(tsls3, vcov = vcovHC(tsls3, type = "HC2"))[, "Std. Error"]

# Combine standard errors and output results with stargazer()
se_models <- append(se_fs2,list(se_tsls2))
stargazer(fs2, tsls2, tsls3, se = se_models, type="text")

## when we removed cold we got a high enough f statistic (around 15)
