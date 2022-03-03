library(stargazer)
library(estimatr)
library(car)

# IMPORT CSV DATA
HPRICE2 <- read.csv("HPRICE2.csv")

# SUMMARY STATISTICS
stargazer(HPRICE2, type="text", digits=1)

# MULTIPLE REGRESSION WITH ROBUST STD ERRORS
model1 <- lm_robust(formula = price ~ nox + rooms + stratio, data = HPRICE2)
summary(model1)

# TEST THAT NOX=0
linearHypothesis(model1,c("nox=0"), white.adjust = "hc2")

# TEST THAT NOX=-1000
linearHypothesis(model1,c("nox=-1000"), white.adjust = "hc2")

# TEST THAT NOX=0 & STRATIO=0
linearHypothesis(model1,c("nox=0", "stratio=0"), white.adjust = "hc2")

# TEST THAT 4*NOX=-ROOMS (just as example)
linearHypothesis(model1,c("4*nox=-rooms"), white.adjust = "hc2")

