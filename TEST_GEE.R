setwd("C:/Users/lutang/Dropbox/RESEARCH_METHOD/MissingPattern/Rwd")

##########################################
## dependent packages
##########################################
library(Matrix)
library(geepack)   # function geese, not as good, DON'T USE
library(gee)       # function gee, better than geepack for HD data, naive and robust variance
library(PGEE)      # Lan Wang's PGEE implementation
library(glmnet)
source("datagenerator.R")
source("mygee.R")


##########################################
## testing penalized GEE
##########################################
# gaussian, binomial, poisson
family_test <- "gaussian"
family_test <- "binomial"
family_test <- "poisson"
set.seed(101)
df <- datagenerator(n=100, m=4, p=5, type=family_test, beta0=c(0.5,0.5,0,0,0), 
                    corst_x="ar1", rho_x=0.3, corst_y="cs", rho_y=0.5)

######## No Penalty ########
summary(gee(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="independence")) # library(gee)
summary(geeglm(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="independence")) # library(geepack)
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence")
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence", weights=rep(1, 400))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence", weights=rep(2, 400))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence", weights=rep(c(2,1), each=200))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence", weights=rep(c(4,2), each=200))


summary(gee(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="exchangeable"))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="exchangeable")

summary(gee(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="AR-M", Mv=1))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], waves=df["visit"], family=family_test, corstr="ar1")

summary(gee(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="unstructured"))
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="unstructured")

# # Syntax for geepack function geeglm (DO NOT USE)
# summary(geeglm(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="independence"))
# summary(geeglm(y~X1+X2+X3+X4+X5-1, data=df, id=id, family=family_test, corstr="exchangeable"))
# summary(geeglm(y~X1+X2+X3+X4+X5-1, data=df, id=id, waves=visit, family=family_test, corstr="ar1"))
# summary(geeglm(y~X1+X2+X3+X4+X5-1, data=df, id=id, waves=visit, family=family_test, corstr="unstructured"))


######## With Penalty, compare with PGee ########
# compare under independence case, no penalty, nuisance parameter might be estimated differently
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, corstr="independence", penalty="SCAD", lambda=0)

pgeefit <- PGee(formula = y ~.-id-visit-1, id = df[,1], data = df, na.action = NULL, 
                family = gaussian(link = "identity"), corstr = "independence", Mv = NULL, scale.fix = FALSE,
                lambda = 0, eps = 10^-6, maxiter = 30, tol = 10^-3, silent = FALSE)
summary(pgeefit)

# compare under CS correlation, no penalty, tau might be estimated differently
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test,
      corstr="exchangeable", penalty="SCAD", lambda=0)

pgeefit <- PGee(formula = y ~.-id-visit-1, id = df[,1], data = df, na.action = NULL, 
                family = gaussian(link = "identity"), corstr = "exchangeable", Mv = NULL,  scale.fix = FALSE,
                lambda = 0, eps = 10^-6, maxiter = 30, tol = 10^-3, silent = FALSE)
summary(pgeefit)

# compare under independence with glmnet, with lasso penalty
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, 
      corstr="independence", penalty="lasso", lambda=1)

glmnetfit <- glmnet(x=as.matrix(df[c("X1","X2","X3","X4","X5")]), y=df[,"y"], intercept=FALSE, 
                    family=family_test, lambda=0.25)
coef(glmnetfit)

# compare under independence with glmnet, with SCAD penalty
mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, 
      corstr="independence", penalty="SCAD", lambda=1.4)

pgeefit <- PGee(formula = y ~.-id-visit-1, id = df[,1], data = df, na.action = NULL, 
                family = gaussian(link = "identity"), corstr = "independence", Mv = NULL, scale.fix = FALSE,
                lambda = 2, eps = 10^-6, maxiter = 30, tol = 10^-3, silent = FALSE)
summary(pgeefit)


# create a solution path for lasso penalty
l_seq <- seq(0, 4, 0.05)
r_seq <- c()
b_seq <- c()
for(l in l_seq){
  print(l)
  fit <- mygee(X=df[c("X1","X2","X3","X4","X5")], y=df["y"], id=df["id"], family=family_test, 
               corstr="independence", penalty="lasso", lambda=l)
  b_seq <- cbind(b_seq, fit$beta[,1])
  r_seq <- c(r_seq ,fit$iteration)
}

plot(b_seq[1,], type="l", ylim=c(0,2))
for(i in 2:nrow(b_seq)){
  lines(b_seq[i,], type="l")
}