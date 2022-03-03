library(estimatr)
library(stargazer)
library(here)
library(tidyverse)

# IMPORT CSV DATA
HPRICE2 <- read_csv(here::here("data", "HPRICE2.csv"))

# SUMMARY STATISTICS
stargazer(HPRICE2, type="text", digits=1)

# BIVARIATE REGRESSION WITH ROBUST STD ERRORS
model1 <- lm_robust(formula = price ~ nox , data = HPRICE2)
summary(model1)

# PREDICTING PRICE WHEN NOX =7 AND DERIVING STD ERROR + 95% CI 

PredPrice=data.frame(nox=c(7))
predict(model1, newdata=PredPrice, se.fit=TRUE, interval='confidence')



# MULTIPLE REGRESSION WITH ROBUST STD ERRORS
model2 <- lm_robust(formula = price ~ nox + rooms, data = HPRICE2)
summary(model2)

# PREDICTING PRICE WHEN NOX=5 & ROOMS=6 AND DERIVING STD ERROR + 95% CI 

PredPrice=data.frame(nox=c(5), rooms=c(6))
predict(model2, newdata=PredPrice, se.fit=TRUE, interval='confidence')
