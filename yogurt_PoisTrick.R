###########################################################################
## Name:         yogurt_PoisTrick.R
## Article:      On the "Poisson Trick" and its Extensions for Fitting 
##               Multinomial Regression Models
## Author:       Jarod Y. L. Lee 
## Last updated: 1 MARCH 2019
## Purpose:      Script for fitting various multinomial models to the 
##               yogurt data, as described in Section 4 of the article.
## Input data:   yogurt.dat                                                  
## Packages:     tidyverse
##               mlogit - fit mlogit as a comparison
##               lme4 - fit models proposed by Chen & Kuo (2001) 
###########################################################################


# Clear R memory ----------------------------------------------------------
rm(list=ls()) 


# Load the required packages ----------------------------------------------
libname <- c("tidyverse","mlogit","lme4")
lapply(libname, require, character.only=T)


# Read raw data -----------------------------------------------------------
yogurt = read.table("yogurt.dat")
colnames(yogurt) = c("id",
                     "yoplait","dannon","weight","hiland",
                     "fy","fd","fw","fh",
                     "py","pd","pw","ph")
head(yogurt)


# Data in wide format -----------------------------------------------------
obs = factor(rep(1:nrow(yogurt)))  # indicator variable: log(delta)

yogurt.wide = cbind(yogurt$id,obs,
                    yogurt$yoplait,yogurt$dannon,yogurt$weight,yogurt$hiland,
                    yogurt$fy,yogurt$fd,yogurt$fw,yogurt$fh,
                    yogurt$py,yogurt$pd,yogurt$pw,yogurt$ph)
colnames(yogurt.wide) = c("id","obs",
                          "yoplait","dannon","weight","control",
                          "fy","fd","fw","fc",
                          "py","pd","pw","pc")
head(yogurt.wide)


# Data in long format - required for “Poisson trick" ----------------------
yogurt.long <- data.frame()
for(i in 1:nrow(yogurt.wide))
{
  new.row1 = yogurt.wide[i,c(1,2,7,11,3)]
  new.row2 = yogurt.wide[i,c(1,2,8,12,4)]
  new.row3 = yogurt.wide[i,c(1,2,9,13,5)]
  new.row4 = yogurt.wide[i,c(1,2,10,14,6)]
  yogurt.long = rbind(yogurt.long,new.row1,new.row2,new.row3,new.row4)
}
brand = rep(c("yoplait","dannon","weight","control"),nrow(yogurt.wide))
yogurt.long = cbind(yogurt.long,brand)
colnames(yogurt.long) = c("id","obs","feature","price","count","brand")

head(yogurt.long)


# Model 1: mlogit ---------------------------------------------------------
yogurt.mlogit = mlogit.data(yogurt.long, choice="count", shape="long", 
                            varying=3:4, alt.var="brand", id.var="id")

# "Hack" to allow for category random intercepts in mlogit
yogurt.mlogit$yop = ifelse(yogurt.mlogit$brand=="yoplait", 1, 0)  
yogurt.mlogit$dan = ifelse(yogurt.mlogit$brand=="dannon", 1, 0)  
yogurt.mlogit$wei = ifelse(yogurt.mlogit$brand=="weight", 1, 0)

# Check the model matrix
head(model.matrix(mFormula(count~yop+dan+wei+feature+price|0),yogurt.mlogit))

m = mlogit(count~yop+dan+wei+feature+price|0, 
           rpar=c(yop="n", dan="n", wei="n"), 
           R=1000, correlation=TRUE, halton=NA, panel=TRUE, 
           data=yogurt.mlogit)
summary(m) # 2 min 22 sec 


# Model 2: gamma-Poisson model --------------------------------------------
data = yogurt.long
id = data$id
obs = data$obs
feature = data$feature
price = data$price
brand = data$brand
y = data$count
data$category = as.numeric(brand)   # create numeric cat var: 1 to Q
I = length(unique(id))              # total number of groups
N = length(unique(obs))             # total number of observations
Q = length(unique(data$category))   # total number of categories

# Create the model matrix
ff = y~-1+factor(obs)+brand+feature+price
utils::str(M <- model.frame(ff, data))
x = model.matrix(ff, M)
head(x[,-c(1:N)])


######### !!! THIS MODEL TOOK AROUND 1.5 DAYS TO CONVERGE !!! #########

############### ECM algorithm ###############
ptm <- proc.time();

# Initialize gamma (result$coef) and alpha
result = glm(y~-1+factor(obs)+brand+feature+price,family="poisson",data=data)  
beta = c(1e-100,rep(1,Q-1))  
#beta = c(1e-100,2.178,6.078,1.941)  

PRINT = TRUE
delta = 1
iteration = 1 
estimates = vector() 
while (delta > .01) 
{
  if(PRINT) print(paste("------------- iteration number ", 
                        as.character(iteration)," -------------"))
  
  result.predicted = exp(x%*%result$coef)
  
  # For each q:
  # (a) calculate pieces needed to estimate the lambda.tilde and loglambda.tilde
  # (b) update beta 
  offset = vector()
  beta.new = vector()
  for (q in 1:Q)
  {
    OK = data$category==q
    Yidot = tapply(y[OK],id[OK],sum)
    SumPredOverj = tapply(result.predicted[OK],id[OK],sum)
    
    # calculate lambda.tilde
    lambda.tilde = (1/beta[q]+Yidot)/(1/beta[q]+SumPredOverj)
    # replicate lambda.tildes and update the offset variable for category q
    lambda.tilde.expanded = rep(lambda.tilde,tapply(id[OK],id[OK],length))
    offset[OK] = log(lambda.tilde.expanded)
    
    # calculate loglambda.tilde
    loglambda.tilde = -log(SumPredOverj+1/beta[q]) + digamma(Yidot+1/beta[q])
    
    # update beta estimate
    beta.objective = function(beta)
    {
      y = (1/beta-1)*sum(loglambda.tilde) - sum(lambda.tilde)/beta -
        I*log(beta)/beta - I*lgamma(1/beta)
      return(-y)
    }
    beta.new[q] = optim(1,beta.objective,method="L-BFGS-B",lower=1e-10,upper=Inf)$par
    if(PRINT) print(paste("lambda.tilde ", as.character(mean(lambda.tilde))))
  }
  beta.new[1] = 1e-100 
  
  # run a new glm with the offset
  result.new = glm(y~-1+factor(obs)+brand+feature+price,
                   offset=offset,family=poisson,data=data)#,maxit=1) 
  
  # test convergence
  delta.delta = (result.new$coef[1:N]-result$coef[1:N])/result$coef[1:N]
  delta.gamma = (result.new$coef[-c(1:N)]-result$coef[-c(1:N)])/result$coef[-c(1:N)]
  delta.beta = (beta.new-beta)/beta
  delta = sum(abs(delta.gamma))+sum(abs(delta.beta))+sum(abs(delta.delta))
  
  # update current estimates
  result = result.new
  beta = beta.new 
  iteration = iteration + 1 
  
  # print current status
  if(PRINT) print(result$coefficients[-c(1:N)])
  if(PRINT) print(paste("beta ", as.character(beta)))
  if(PRINT) print(paste("delta ", as.character(delta)))
  estimates = rbind(estimates,c(result$coefficients[-c(1:N)],beta))
}
proc.time() - ptm;            # 38.8 hours

iteration = iteration - 1
iteration                     # number of iterations
c(result$coef[-c(1:N)],beta)  # estimated coefficients
# branddannon    brandweight   brandyoplait   feature       price                               
# 4.614435e+00   3.677556e+00   5.274746e+00   7.845186e-01  -4.088184e+01   
# b1             b2            b3             b4
# 1.000000e-100   2.200877e+00   6.067862e+00   1.920390e+00


############### Check gradient = 0 ###############
## Score equation
choices = c("control","dannon","weight","yoplait")
yPLUS.qi = matrix(0,Q,I)        # yPLUS.qi[q,i] is y_{i+q}
for(q in 1:Q)
{
  yPLUS.qi[q,] = tapply(y[brand==choices[q]],id[brand==choices[q]],sum)
} 

grad <- function(theta)
{
  log.delta = theta[1:N]
  gamma = theta[(N+1):(N+5)]   # Yoplait, Dannon, Weight, Feature, Price
  b = theta[(N+6):(N+9)]  
  a = 1/b        
  
  # temp.qi[q,i] is \sum_j e^{x_{ijq} \gamma}
  temp1 = exp(x%*%c(log.delta,gamma))
  temp.qi = matrix(0,Q,I)  
  for(q in 1:Q)
  {
    temp.qi[q,] = tapply(temp1[brand==choices[q]],id[brand==choices[q]],sum)
  }
  
  comp.g1 = vector(length=(N+5))
  for(i in 1:I)
  {
    for(q in 1:Q)
    {
      # temp2 is sum_j x_{ijq} e(x_{ijq}\gamma)
      condition = id==i&brand==choices[q]
      temp2 = apply(x[condition,]*temp1[condition],2,sum)
      comp.g1 = comp.g1+(a[q]+yPLUS.qi[q,i])*(temp2)/(1/b[q]+temp.qi[q,i])
    }
  }
  
  comp.g2 = apply(x*y,2,sum)
  comp.g = -comp.g1 + comp.g2
  
  comp.b = -apply(digamma(1/b+yPLUS.qi),1,sum)/b^2 + 
    apply(((1/b^3+yPLUS.qi/b^2)/(1/b+temp.qi)),1,sum) +
    apply(log(1/b+temp.qi),1,sum)/b^2 +
    I*digamma(1/b)/b^2 - I/b^2 + I*log(b)/b^2
  
  grad = c(comp.g,comp.b)
  return(grad)
}

theta=c(result$coef,beta)
summary(grad(theta))
# Min.       1st Qu.     Median     Mean       3rd Qu.    Max. 
# -3.853e-04 -4.571e-06  1.659e-06  1.618e-06  5.436e-06  1.290e-03 


############### Standard errors ###############
# Hessian matrix
comp.gg <- function(theta)
{
  log.delta = theta[1:N]
  gamma = theta[(N+1):(N+5)] # Yoplait, Dannon, Weight, Feature, Price
  b = theta[(N+6):(N+9)]  
  a = 1/b        
  
  # temp.qi[q,i] is \sum_j \delta_{ij}\zeta_{ijq}
  temp = exp(x%*%c(log.delta,gamma))
  temp.qi = matrix(0,Q,I)  
  for(q in 1:Q)
  {
    temp.qi[q,] = tapply(temp[brand==choices[q]],id[brand==choices[q]],sum)
  }
  
  comp.gg = matrix(0,nrow=N+5,ncol=N+5)
  for(i in 1:I)
  {
    for(q in 1:Q)
    {
      condition = id==i&brand==choices[q]
      
      comp1 = 1/b[q]+temp.qi[q,i]
      
      # comp2 is sum_j x_{ijq} e(x_{ijq}\gamma)
      comp2 = apply(x[condition,]*temp[condition],2,sum)
      comp2 = comp2 %*% t(comp2)
      
      J = sum(condition)
      comp3 = matrix(0,nrow=N+5,ncol=N+5)
      for(j in 1:J)
      {
        comp3 = comp3 + x[condition,][j,]%*%t(x[condition,][j,])*temp[condition][j]
      }
      
      comp.gg = comp.gg - (1/b[q]+yPLUS.qi[q,i])*(comp1*comp3-comp2)/(comp1)^2
    }
  }
  
  return(comp.gg)
}


comp.bb <- function(theta)
{
  log.delta = theta[1:N]
  gamma = theta[(N+1):(N+5)]   # Yoplait, Dannon, Weight, Feature, Price
  b = theta[(N+6):(N+9)]  
  a = 1/b        
  
  # temp.qi[q,i] is \sum_j e^{x_{ijq} \gamma}
  temp = exp(x%*%c(log.delta,gamma))
  temp.qi = matrix(0,Q,I)  
  for(q in 1:Q)
  {
    temp.qi[q,] = tapply(temp[brand==choices[q]],id[brand==choices[q]],sum)
  }
  
  comp.bb = vector()
  for(q in 2:Q)
  {
    comp1 = (1/b[q]+temp.qi[q,])*(-3/b[q]-2*yPLUS.qi[q,])+(1/b[q]^2+yPLUS.qi[q,])
    comp2 = b[q]^3*(1/b[q]+temp.qi[q,])^2
    
    comp.bb[q-1] = sum(trigamma(1/b[q]+yPLUS.qi[q,]))/b[q]^4 +
      2*sum(digamma(1/b[q]+yPLUS.qi[q,]))/b[q]^3 +
      sum(comp1/comp2) -
      2*sum(log(1/b[q]+temp.qi[q,]))/b[q]^2 -
      (sum(1/b[q]+temp.qi[q,])*b[q]^4)^-1 - 
      I*trigamma(1/b[q])/b[q]^4 - 
      2*I*digamma(1/b[q])/b[q]^3 +
      3*I/b[q]^3 - 
      2*I*log(b[q])/b[q]^3
  }
  
  return(diag(comp.bb))
}


comp.gb <- function(theta)
{
  log.delta = theta[1:N]
  gamma = theta[(N+1):(N+5)]   # Yoplait, Dannon, Weight, Feature, Price
  b = theta[(N+6):(N+9)]  
  a = 1/b        
  
  # temp.qi[q,i] is \sum_j e^{x_{ijq} \gamma}
  temp = exp(x%*%c(log.delta,gamma))
  temp.qi = matrix(0,Q,I)  
  for(q in 1:Q)
  {
    temp.qi[q,] = tapply(temp[brand==choices[q]],id[brand==choices[q]],sum)
  }
  
  comp.gb = matrix(0,nrow=N+5,ncol=Q-1)
  for(q in 2:Q)
  {
    for(i in 1:I)
    {
      condition = id==i&brand==choices[q]
      # comp1 is sum_j x_{ijq} e(x_{ijq}\gamma)
      comp1 = apply(x[condition,]*temp[condition],2,sum)
      comp2 = comp1*(yPLUS.qi[q,i]-temp.qi[q,i])
      comp3 = b[q]^2*(temp.qi[q,i]+1/b[q])^2
      comp.gb[,q-1] = comp.gb[,q-1] - comp2/comp3
    }
  }
  
  return(comp.gb)
}


hess <- function(theta)
{
  comp_gg = comp.gg(theta)
  comp_gb = comp.gb(theta)
  comp_bg = t(comp_gb)
  comp_bb = comp.bb(theta)
  
  hess1 = cbind(comp_gg,comp_gb)
  hess2 = cbind(comp_bg,comp_bb)
  hess = rbind(hess1,hess2)
  return(hess)
}

theta=c(result$coef,beta)
se = diag(sqrt(-solve(hess(theta))))
se[-c(1:N)]
# Yoplait   Dannon    Weight    Feature   Price
# 0.3088661 0.3918856 0.3423029 0.1777258 3.7779590
# b2        b3        b4
# 0.1335338 0.3739458 0.1349498


# Model 3: Chen & Kuo (2001) ----------------------------------------------

########## !!! THESE MODELS DID NOT CONVERGE WITHIN A MONTH !!! ##########

# Simplified Poisson log-linear model proposed by Chen & Kuo (2001), with
# just a random intercept per household.
ptm <- proc.time()   
poisglmm1 = glmer(count~-1+factor(obs)+brand+feature+price+(1|id),
                  family="poisson",
                  data=yogurt.long)
proc.time() - ptm 
summary(poisglmm1)


### Fit model by Chen and Kuo (2001), treating brand as nested within household
ptm <- proc.time()   
poisglmm2 = glmer(count~-1+factor(obs)+brand+feature+price+(1|brand/id),family="poisson",data=yogurt.long)
proc.time() - ptm  
summary(poisglmm2)



