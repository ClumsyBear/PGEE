library(MASS)
library(SimCorMultRes)
# library(bindata)

##########################################
## data generator
##########################################
workingCor <- function(p, corst, rho = 0) {
  if (corst == "ind") {Sigma <- diag(rep(1, p))}
  if (corst == "cs") {Sigma <- matrix(rho, p, p); diag(Sigma) <- 1}
  if (corst == "ar1") {Sigma <- rho^(abs(outer(1:p, 1:p, "-")))}
  return(Sigma)
}

# to generate a longitudinal data set
# test input: 
# n = 100; m = 4; p = 5; type = "gaussian"; beta0=c(0.5,0.3,0.2,0,0); 
# intercept=FALSE; corst_x="ind"; rho_x=0; corst_y="cs"; rho_y=0.6; seed=1989;
dataGenerator <- function(n, m, p, type, intercept = FALSE, beta0, 
                          corst_x = "ind", rho_x = 0, corst_y = "ind", rho_y = 0, 
                          seed = NULL){
  n          # sample size
  m          # number of visits
  p          # number of covariates (not including intercept)
  type       # c("gaussian", "binomial", "poisson")
  beta0      # coefficients, length equals to p if intercept=FALSE, and p+1 if intercept=TRUE
  intercept  # if TRUE, the beta0[1] is used for the coefficient of intercept
  corst_x    # c("ind", "cs", "ar1")
  rho_x      # parameter of correlation
  corst_y    # c("ind", "cs", "ar1")
  rho_y      # parameter of correlation
  seed       # random seed
  if(!is.null(seed)) {set.seed(seed)}
  if(length(m) != 1) {stop("m must be an integer")}
  
  Sigma_x <- workingCor(p = p, corst = corst_x, rho = rho_x)
  Sigma_y <- workingCor(p = m, corst = corst_y, rho = rho_y)
  id <- rep(1:n, each = m)
  visit <- rep(1:m, times = n)
  
  X <- data.frame(mvrnorm(n = n * m, mu = rep(0, p), Sigma = Sigma_x))
  X[, 1] <- (X[, 1] < 0.0) + 0
  X[, 2] <- (X[, 2] < 0.0) + 0
  # X[, p - 1] <- rep(rnorm(n, 0, 1), each = m )# create a time-invariant continuous covariate
  # X[, p] <- rep(rbinom(n, 1, 0.5), each = m) # create a time-invariant binary covaraite
  
  if(intercept == TRUE) {X <- cbind(X0 = 1, X)}

  Xbeta0 <- drop(as.matrix(X) %*% beta0)
  if(type == "gaussian") {
    epsilon <- mvrnorm(n = n, mu = rep(0, m), Sigma = Sigma_y)
    epsilon <- c(t(epsilon))
    y <- Xbeta0 + epsilon
  }
  if(type == "binomial") {
    y <- c(t(rbin(clsize = m, intercepts = 0, betas = beta0,
                  xformula = ~ X1 + X2 + X3 + X4 + X5, xdata = X,
                  cor.matrix = Sigma_y, link = "logit")$Ysim))
  }
  # if(type=="poisson") {y <- rpois(n*m, lambda=exp(Xbeta0 + epsilon))}
  
  return(data.frame(id = id, visit = visit, y = y, X))
}