##########################################
## mygee function
##########################################
# description: my gee function
# test case: 
# X=df[c("X1","X2","X3","X4","X5")]; y=df["y"]; id=df["id"]; family=family_test; corstr="independence"; intercept=FALSE; penalty=NULL; lambda=0; penalty_factor=NULL; weights = NULL; beta_warm = NULL
# X=df[c("X1","X2","X3","X4","X5")]; y=df["y"]; id=df["id"]; family=family_test; corstr="exchangeable"; intercept=FALSE; penalty=NULL; lambda=NULL; penalty_factor=NULL
# X=df[c("X1","X2","X3","X4","X5")]; y=df["y"]; id=df["id"]; family=family_test; corstr="ar1"; waves=df["visit"]; intercept=FALSE; penalty=NULL; lambda=NULL; penalty_factor=NULL

##########################################
library(MASS)
# library(Rcpp)

mygee <- function (X, y, id, family, intercept = FALSE, corstr = "independence", 
                   waves = NULL, penalty = NULL, lambda = 0, 
                   penalty_factor = NULL, weights = NULL, beta_warm = NULL){
  if(!(family %in% c("gaussian", "binomial", "poisson"))) {stop("'family' value not supported")}
  if(!(corstr %in% c("independence", "exchangeable", "ar1", "unstructured"))) {stop("'corstr' value not supported")}
  if(lambda != 0) {
    if(is.null(penalty)) {stop("need to specify 'penalty' when lambda is not 0")
    } else{if(!(penalty %in% c("lasso", "SCAD", "MCP"))) {stop("'penalty' value not supported")}}
  }
  if(is.null(weights)) {weights <- rep(1, nrow(X))} else {if(length(weights) != nrow(X)) {stop("'weights' length should equal nrow(X)")}}
  if(class(y) == "data.frame") {y <- y[,1]}
  if(class(id) == "data.frame") {id <- id[,1]}
  if(class(X) == "data.frame") {X <- as.matrix(X)}
  if(intercept) {X <- cbind(Int = 1, X)}
  n <- length(unique(id))
  p <- ncol(X)
  N <- nrow(X)
  maxit <- 400 # maximum iteration number, 30 in PGEE
  tol <- 1e-04 # tolerance for convergence, 1e-04 in PGEE
  count <- 0 # iteration count
  stop_flag <- FALSE
  converged <- FALSE
  epsilon <- 1e-06 # epsilon for MM algorithm, 1e-6 in PGEE
  if(is.null(beta_warm)) {beta_old <- coef(glm(y ~ . - 1, data = data.frame(y = y, X), family = family, weights = weights))
  } else {beta_old = beta_warm}
  # cat("glm initial regression estimate \n"); print(beta_old)
  InvLink <- getLinkFuns(family = family)$InvLink
  VarFunc <- getLinkFuns(family = family)$VarFunc
  if(lambda != 0) {
    if(is.null(penalty_factor)) {penalty_factor <- rep(1, p)}
    lambda_long <- lambda * penalty_factor
    QFunc <- getQFunc(penalty = penalty, lambda_long = lambda_long)
  } else {
    QFunc <- 0
  }
  while (!stop_flag) {
    count <- count + 1
    #### update \phi ####
    eta <- drop(X %*% beta_old)
    residual <- (y - InvLink(eta)) / sqrt(VarFunc(eta))
    phi <- sum(residual^2) / (N - p) # same as in gee, different from geese which divides by N
    #### update \tau ####
    tau <- updateTau(residual = residual, id = id, N = N, p = p, n = n, phi = phi, corstr = corstr)
    #### update \beta ####
    S_accum <- 0
    H_accum <- 0
    G_accum <- 0
    # SS_accum <- 0
    for(i in unique(id)) {
      selector <- id == i
      mi <- sum(selector)
      Xi <- subset(X, selector)
      yi <- y[selector]
      mui <- InvLink(eta[selector])   
      nui <- sqrt(VarFunc(eta[selector]))   # THIS NEEDS TO MULTIPLY sqrt(phi)
      if(mi > 1){
        Ai_half <- diag(nui)
        Ai_inv_half <- diag(1/nui)
        wi <- diag(weights[selector])
      } else{
        Ai_half <- nui
        Ai_inv_half <- 1/nui
        wi <- weights[selector]
      }
      Ri <- getRi(tau = tau, mi = mi, corstr = corstr)
      invRi <- ginv(Ri)
      varyi <- (yi - mui) %*% t(yi - mui)
      S_accum <- S_accum + t(Xi) %*% Ai_half %*% invRi %*% Ai_inv_half %*% wi %*% (yi - mui)
      # SS_accum <- SS_accum + outer(as.vector(t(Xi) %*% Ai_half %*% invRi %*% Ai_inv_half %*% wi %*% (yi - mui)), as.vector(t(Xi) %*% Ai_half %*% invRi %*% Ai_inv_half %*% wi %*% (yi - mui)))
      H_accum <- H_accum + t(Xi) %*% Ai_half %*% invRi %*% Ai_half %*% wi %*% Xi
      G_accum <- G_accum + t(Xi) %*% wi %*% Ai_half %*% invRi %*% Ai_inv_half %*% varyi %*% Ai_inv_half %*% invRi %*% Ai_half %*% wi %*% Xi 
    }
    if(lambda == 0) {
      beta_new <- beta_old + drop(ginv(H_accum) %*% (S_accum))
    } else {
      E_term <- diag(drop(QFunc(pmax(0, abs(beta_old))) / (epsilon + abs(beta_old))))
      beta_new <- beta_old + drop(ginv(H_accum + n * E_term) %*% (S_accum - n * E_term %*% beta_old)) 
    }
    # if(max(abs((beta_new - beta_old)/(beta_old + .Machine$double.eps))) < tol) {converged <- TRUE; stop <- TRUE}
    if(sum(abs((beta_new - beta_old))) < tol) {converged <- TRUE; stop_flag <- TRUE}
    if(count >= maxit) {stop_flag <- TRUE}
    beta_old <- beta_new
  } # end of while loop
  if(converged == FALSE) {warning("algorithm reached 'maxit' but did not converge")}
  if(lambda == 0) {
    Godambe <- H_accum %*% ginv(G_accum) %*% H_accum
    varb <- sqrt(diag(ginv(Godambe)))
    out_beta <- cbind(Estimate = beta_new, StdErr = varb, Zscore = beta_new / varb, pvalue = 2 * pnorm(-abs(beta_new / varb)))
  } else {
    beta_new[abs(beta_new) < 1e-2] <- 0
    out_beta <- cbind(Estimate = beta_new)
  }
  out <- list()
  out$beta <- out_beta
  out$phi <- phi
  if (count > 50) {print(paste0("number of PGEE iteration = ", count))}
  out$iteration <- count
  out$corstr <- corstr
  out$Ri <- Ri
  out$converged <- converged
  if (lambda == 0) {out$varbinv <- Godambe}
  return(out)
}

getLinkFuns <- function(family) {
  if (family == "gaussian") { 
    InvLink <- function(x) {x} 
    VarFunc <- function(x) {rep(1, length(x))} 
  }
  if (family == "binomial") { 
    InvLink <- function(x) {exp(x) / (1 + exp(x))}
    VarFunc <- function(x) {exp(x) / (1 + exp(x))^2} 
  }
  if (family == "poisson") { 
    InvLink <- VarFunc <- function(x) {exp(x)} 
  }
  return(list(InvLink = InvLink, VarFunc = VarFunc))
}

getQFunc <- function(penalty, lambda_long) {
  if(penalty == "lasso") { 
    QFunc <- function(x) {lambda_long * rep(1,length(x))} 
  }
  if(penalty == "SCAD") { 
    a = 3.7 # parameter in SCAD
    QFunc <- function(x) { 
      output <- lambda_long * ((x <= lambda_long) + (pmax(0, a * lambda_long - x) / ((a - 1) * lambda_long)) * (x > lambda_long)) 
      output[is.na(output)] <- 0
      output
    } 
  }
  if(penalty == "MCP") {
    a = 1.5 # parameter in MCP
    QFunc <- function(x) {
      output <- lambda_long * pmax(0, a * lambda_long - x) / (a * lambda_long)
      output[is.na(output)] <- 0
      output
    }
  }
  return(QFunc)
}

updateTau <- function(residual, id, N, p, n, phi, corstr) {
  if(corstr == "independence") {tau <- 0}
  if(corstr == "exchangeable") { 
    temp_r_accum <- c()
    for (i in unique(id)) {
      selector <- id == i
      residuali <- residual[selector] %*% t(residual[selector])
      temp_r_accum <- c(temp_r_accum, residuali[upper.tri(residuali)])
    }
    # need to standardize by phi to get true correlation parameter
    tau <- sum(temp_r_accum) / (length(temp_r_accum) - p) / phi # should this be minus p or what???
  }
  if(corstr == "ar1") { # requires 'waves' to work better
    temp_r_accum <- c()
    for (i in unique(id)) {
      selector <- id == i
      residuali <- residual[selector] %*% t(residual[selector])
      temp_r_accum <- c(temp_r_accum, diag(as.matrix(residuali[-nrow(residuali), -1]))) # only using subdiagonal
      # distancei <- as.matrix(dist(waves[selector], method = "manhattan")) # did not use waves, can improve
      # temp_t_accum <- c(temp_t_accum, diag(distancei[-nrow(distancei), -1]))
    }
    tau <- sum(temp_r_accum) * N / (N - p) / (length(temp_r_accum)) / phi # should this be minus p or what???
  }
  if(corstr == "unstructured") {
    temp_r_accum <- c()
    for (i in unique(id)) {
      selector  <- id == i
      residuali <- residual[selector] %*% t(residual[selector])
      temp_r_accum <- rbind(temp_r_accum, residuali[upper.tri(residuali)]) # there is error here, dimension how correct.
    }
    tau <- colSums(temp_r_accum) * N / (N - p) / n / phi # should this be minus p or what???
  }
  return(tau)
}

getRi <- function(tau, mi, corstr){
  if (corstr == "independence") { Ri <- diag(rep(1, mi)) }
  if (corstr == "exchangeable") { Ri <- matrix(tau, mi, mi); diag(Ri) <- 1 }
  if (corstr == "ar1") { Ri <- tau^(abs(outer(1:mi, 1:mi, "-"))) }
  if (corstr == "unstructured") { 
    Ri <- diag(rep(1, mi))
    ut <- upper.tri(Ri, diag=FALSE)
    lt <- lower.tri(Ri, diag=FALSE)
    Ri[ut] <- tau[order(col(ut)[ut], col(ut)[ut])]
    Ri[lt] <- tau[order(row(ut)[ut], col(ut)[ut])]
  }
  return(Ri)
}


# improve loops by writing C++
# compare initial value with PGEE
# don't use t(X) V X
