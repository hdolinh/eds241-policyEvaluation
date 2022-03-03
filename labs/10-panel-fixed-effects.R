library(stargazer)
library(estimatr)
library(ggplot2)
library(plm)
library(sandwich)
library(lmtest)

# SET WORKING DIRECTORY

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Automatically sets directory where script is saved
getwd()

# IMPORT CSV DATA and CREATE SUBSET for 1988
SWTF <- read.csv("SW_TrafficFatality.csv")
SWTF88 <- subset(SWTF, year == "1988")
SWTF82 <- subset(SWTF, year == "1982")


# SUMMARY STATISTICS
stargazer(SWTF, type="text", digits=2)


# BASIC OLS REGRESSION
m88 <- lm(formula = frate ~ beertax, data=SWTF88)

se_m88 <- starprep(m88, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 
stargazer(m88, se = se_m88, type="text")


# SCATTERPLOT with REGRESSION LINE
ggplot(data=SWTF88, aes(x=beertax, y=frate)) +
  geom_smooth(method=lm) +
  geom_text(aes(label=state), size = 3) 


# FIRST DIFFERENCE MODEL FOR 1988 AND 1982
dfrate <- SWTF88$frate - SWTF82$frate
dbtax <- SWTF88$beertax - SWTF82$beertax
fd8882 <- lm(formula = dfrate ~ dbtax)

se_fd8882 <- starprep(fd8882, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 
stargazer(fd8882, se = se_fd8882, type="text")




# ESTIMATE THE BASIC WITHIN FIXED EFFECTS REGRESSION
# NOTE "plm" ONLY PRODUCES CLUSTER-ROBUST STANDARD ERRORS

within1 <- plm(frate ~ beertax, index = c("state", "year"), model = "within", effect = "twoways", data = SWTF)

# Calculate standard errors (note slightly different procedure with plm package)
se_within1 <- coeftest(within1, vcov = vcovHC(within1, type = "HC2"))[, "Std. Error"]
# Reformat standard errors for stargazer()
se_within1 <- list(se_within1)
# Output results with stargazer
stargazer(within1, keep=c("beertax"), se = se_within1, type="text")



# LEAST SQUARES DUMMY VARIABLES REGRESSION
# COMPARE HETEROSKEDASTICITY AND CLUSTER-ROBUST STANDARD ERRORS
lsdv1 <- lm(formula = frate ~ beertax + as.factor(state) + as.factor(year), data=SWTF)
lsdv2 <- lm(formula = frate ~ beertax + as.factor(state) + as.factor(year), data=SWTF)

se_lsdv_Hrobust <- starprep(lsdv1, stat = c("std.error"), se_type = "HC2", alpha = 0.05) 
se_lsdv_Crobust <- starprep(lsdv2, stat = c("std.error"), se_type = "CR2", clusters=SWTF$state, alpha = 0.05) 
se_list <- list(se_lsdv_Hrobust[[1]], se_lsdv_Crobust[[1]])

stargazer(lsdv1, lsdv2, keep=c("beertax"), se = se_list, type="text")

