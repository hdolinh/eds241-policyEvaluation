library(MASS) # random number generator for bivariate or multivariate random variables
library(ggplot2)
library(vtable) # useful for summary statistics? 
library(stargazer) # tables
library(estimatr) # lm_robust 

set.seed(7307) # don't have to specify R will use the time or date if you don't set this

bigN=20000 # number of obs we will use

W <- runif(bigN,0,5) # w is random uniform variable so the probability of drawing a number between 0 and 5 is equal
X=as.integer(W+1)
X1 <- as.numeric(X==1) # defining numeric variables 
X2 <- as.numeric(X==2) # double equal sign is an indicator function / logical 
X3 <- as.numeric(X==3) # return value of 1 when x is true and 0 when it's false
X4 <- as.numeric(X==4) # these can be considered "binary variables"
X5 <- as.numeric(X==5)



# GENERATE MEAN COMPONENT OF POTENTIAL OUTCOMES
MU0=(1/2)*X1 + (2/2)*X2 + (3/2)*X3 + (4/2)*X4 + (5/2)*X5
mean(MU0)
MU1=1*X1 + 2*X2 + 3*X3 + 4*X4 + 5*X5
mean(MU1)


# GENERATE ERROR COMPONENT OF POTENTIAL OUTCOMES
Sigma <- matrix(c(1,0.75,0.75,1),2,2) # created a covariance matrix 
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

ATE <- mean(Y1)-mean(Y0) # we are controlling the average treatment effect here
print(ATE)

PO_DF <- data.frame(Y0,Y1,X) # potential outcome dataframe

# PLOT POTENTIAL OUTCOMES AGAINST X, JITTER TO CLARIFY VISUAL
jitter <- position_jitter(width = 0.05, height = 0)
ggplot(PO_DF, aes(x=X, y=Y0)) +
  geom_point(color = "blue") +
  geom_point(position = jitter, color = "red", shape=1, aes(x=X, y=Y1))
# x only takes 5 discrete values 
# we have symbols piled on top of each other - use jitter to distinguish between y(1) and y(0)
# y(1) in red and y(0) in blue
# jitter spreads out the y(1) data points a little more so it's easier to see
## graph takeaway - wanted to highlight relationship between x and the potential outcomes but in reality we never observe potential outcomes so no big takeaway here

## We want to get an ATE of 1.5 ##
 
### SUBJECTS ARE ARRIVING ###

# RANDOMLY ASSIGN A TREATMENT INDICATOR
D <- as.numeric((runif(bigN,0,1)) > 0.5) # runif says draw a random number between 0 and 1, if it's less than 5 you're not in control group aka 0 and vice versa 
mean(D)
# we expect it to be 0.5 so half our obs received the a value of 1 (treatment) and received a value of 0 (no treatment) 


# USE SUTVA TO MAP POTENTIAL OUTCOMES INTO OBSERVED OUTCOMES
Y = D*Y1 + (1-D)*Y0


# COLLECT ALL RELEVANT VARIABLES IN A DATAFRAME
RCT_DATA <- data.frame(Y, D, Y0, Y1, X, X1, X2, X3, X4, X5) # added Y here



# CHECK THAT D IS INDEPENDENT OF X, Y0, Y1 (RECALL Y0,Y1 NOT OBSERVED IN REALITY) 
# "TEST" OF COVARIATE BALANCE - slide 16 in Lecture 5
sumtable(RCT_DATA, vars=c('Y0','Y1', 'Y', 'X1', 'X2', 'X3', 'X4', 'X5'), group='D', group.test=TRUE) # sumtable tests covariate balance
# sumtable is a vtable function
# roughly half the sample with D(0) and D(1) we see from N that it's almost split 50/50
# get an F test value returned we reject when we have a large test value
# we don't reject the null hypothesis here bc none of the F test are significant 
# we can see the balance property here holds (means match and so do N)
# we fail to reject the null hypothesis here in all cases
# we want this and no significance 

## added Y to sumtable and see that value was significant ##

mA <- lm(formula = X ~ D, data=RCT_DATA) # we should fail to reject the null hypothesis again
mB <- lm(formula = Y0 ~ D, data=RCT_DATA)
mC <- lm(formula = Y1 ~ D, data=RCT_DATA)
se_models = starprep(mA, mB, mC, stat = c("std.error"), se_type = "HC2", alpha = 0.05) # use starprep to get heteroskad and robust SDs
stargazer(mA, mB, mC, se = se_models, type="text")
# regress x (1-5) on D
# difference in means is almost 0


# ESTIMATE ATE USING SIMPLE OLS REGRESSION OF Y on D
ate1 <- lm(formula = Y ~ D, data=RCT_DATA)
ate2 <- lm(formula = Y ~ D + X, data=RCT_DATA)
se_models = starprep(ate1, ate2, stat = c("std.error"), se_type = "HC2", alpha = 0.05)
stargazer(ate1, ate2, se = se_models, type="text")

# key regressor is D
# coefficient from model1 to model2 changes but not in a meaningful way
# SD error decreases in the second model because 
# see that adjusted r2 also changes quite the jump

