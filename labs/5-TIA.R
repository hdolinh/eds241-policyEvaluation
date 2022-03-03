library(MASS)
library(ggplot2)
library(vtable)
library(stargazer)
library(estimatr)
library(dplyr)
library(tidyr)

set.seed(7307)

bigN=20000

W <- runif(bigN,0,5)
X=as.integer(W+1)
X1 <- as.numeric(X==1)
X2 <- as.numeric(X==2)
X3 <- as.numeric(X==3)
X4 <- as.numeric(X==4)
X5 <- as.numeric(X==5)



# GENERATE MEAN COMPONENT OF POTENTIAL OUTCOMES
MU0=(1/2)*X1 + (2/2)*X2 + (3/2)*X3 + (4/2)*X4 + (5/2)*X5
mean(MU0)
MU1=1*X1 + 2*X2 + 3*X3 + 4*X4 + 5*X5
mean(MU1)


# GENERATE ERROR COMPONENT OF POTENTIAL OUTCOMES
Sigma <- matrix(c(1,0.75,0.75,1),2,2)
Sigma
e <- (mvrnorm(n=bigN, rep(0, 2), Sigma))
e0 <- e[,c(1)]
mean(e0)
e1 <- e[,c(2)]  
mean(e1)



# GENERATE POTENTIAL OUTCOMES
Y0 <- MU0 + e0
mean(Y0)
Y1 <- MU1 + e1
mean(Y1)

ATE <- mean(Y1)-mean(Y0)
print(ATE) # verify treatment effect, so ATE should be 1.5

# ASSIGN A TREATMENT IGNORABLE COND ON X
v <- rnorm(bigN, 0, 2) # random variable v, random erro from treatment to control
D <- as.numeric((-2*X1+ -1*X2 + 0.5*X3 + 1*X4 + 2*X5 + v)>0) # if all is greater than zero then assign to the treatment group; if not assign to the control group
mean(D)


# USE SUTVA TO MAP POTENTIAL OUTCOMES INTO OBSERVED OUTCOMES
Y = D*Y1 + (1-D)*Y0 
# the means are no longer the same

# COLLECT ALL RELEVANT VARIABLES IN A DATAFRAME
TIA_DATA <- data.frame(D, Y0, Y1, X, X1, X2, X3, X4, X5)


# SHOW THAT D IS NOT INDEPENDENT OF X, Y0, Y1 (RECALL Y0,Y1 NOT OBSERVED IN REALITY)
sumtable(TIA_DATA, vars=c('Y0','Y1', 'X'), group='D', group.test=TRUE)
# we see that the means have changed; those with the higher means will be assigned the treatment
# will get a biased regression we see that the dependent variable Y has a D of ~2.5 when it should be 1.5

# SIMPLE OLS ESTIMATOR NOT CONSISTENT FOR ATE
ols <- lm(formula = Y ~ D, data=TIA_DATA)
se_models = starprep(ols, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(ols, se = se_models, type="text")


# DIFFERENCE IN MEAN OF OUTCOME Y FOR D=1 AND D=0, BY CELL OF X
sumtable(TIA_DATA, vars=c('Y', 'X1', 'X2', 'X3', 'X4', 'X5'), group='D', group.test=TRUE)


TIA_table <- TIA_DATA %>%
  mutate(Y = ifelse(D==1, Y1,Y0))%>% #Create observed Y variable
  group_by(X,D)%>% 
  summarise(n_obs = n(),
            Y_mean= mean(Y, na.rm = T))%>% #Calculate number of observations and Y mean by X by treatment cells
  gather(variables, values, n_obs:Y_mean)%>% #Reshape data
  mutate(variables = paste0(variables,"_",D, sep=""))%>% #Combine the treatment and variables for reshaping
  pivot_wider(id_cols = X, names_from = variables,values_from = values)%>% #Reshape data by treatment and X cell
  ungroup()%>%  #Ungroup from X values
  mutate(Y_diff = Y_mean_1 - Y_mean_0, #calculate Y_diff
         w_ATE = (n_obs_0+n_obs_1)/(sum(n_obs_0)+sum(n_obs_1)),
         w_ATT = n_obs_1/sum(n_obs_1))%>% #calculate weights
  mutate_if(is.numeric, round, 2) #Round data

# y diff is delta of x
# delta of x1 is 0.062 (mean)
# the means are what are helpful
# notice that x_obs_0 and n_obs_1 are inverse of each other
# y diffs are our four deltas
# average treatment effect is Y_diff...look below

stargazer(TIA_table, type= "text", summary = FALSE, digits = 2)

# MULTIVARIATE MATCHING ESTIMATES OF ATE AND ATT
ATE=sum((TIA_table$w_ATE)*(TIA_table$Y_diff))
ATE # this should be 1.5, it is about 1.5 it's 1.48
ATT=sum((TIA_table$w_ATT)*(TIA_table$Y_diff))
ATT # 1.83 didn't have a theoretical 

# MULTIVARIATE MATCHING AS REGRESSION ESTIMATOR
reg_ate <- lm(formula = Y ~ D + X2 + X3 + X4 + X5, data=TIA_DATA)
se_models = starprep(reg_ate, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(reg_ate, se = se_models, type="text")
# D is 1.49
# regression here is interpreted as difference of means



