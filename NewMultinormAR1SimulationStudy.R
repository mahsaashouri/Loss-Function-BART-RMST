library(survival)
library(BART)
library(glmnet)
library(mboost)
library(BARTTrial)
library(AFTrees)
source("DrawIPCW.R")

set.seed(1234)

ndraws <- 1000
burnIn <- 100
sigma <- 1.0
n <- 2000 # 250 or 2000 # number of training observation
n.test <- 4000 # 2000 or 4000 # number of test observation
num_covar <- 10 # 10 or 100 # total number of predictors
coef <- c(runif(5, 0, .5), rep(0, num_covar-5))
Rho <- 0.5
nreps <- 100 # number of simulation replications
## choosing this big tau value cause warning
gam_alpha <- 20


## function to simulate linear model AR1
sim.reg <- function(nobs, coef, mu, sd, Rho){
  num.var = length(coef)
  beta = as.matrix(coef)
  ## generate data from multivariate normal with AR(1) - using the above function and MASS package
  #H = mvrnorm(n = nobs, mu, Sigma = AR1Cor(num.var, Rho))
  ## generate data from multivariate normal with AR(1) - using 'arima.sim' function
  H = t(replicate(nobs, c(arima.sim(length(coef), model = list(ar = Rho)))))
  Y = exp(H %*% beta)
  return(list(Z = H, Y = c(Y), coeff = coef))
}

CoxExpectedSurv <- function(X, beta_val, time, H0.vals, tau) {
  ## This function computes E( min(T_i, tau) |x_i) for
  ## a cox ph model
  mu.x <- colMeans(X)
  tpoints <- c(time[time < tau], tau)
  Hpoints <- H0.vals[time < tau]
  nn <- nrow(X)
  
  fitted_vals <- rep(NA, nn)
  for(k in 1:nn) {
    nu <- sum((X[k,] - mu.x)*beta_val)
    fitted_vals[k] <- sum(exp(-Hpoints*exp(nu))*diff(tpoints))
  }
  return(fitted_vals)
}



cens_rate <- 1.5 # Use 0.25 (high censoring) or 1.5 (low censoring)
tau <- 25
sgrid <- seq(0, tau, by=.1)

cens_prop <- rep(NA, nreps)
CorCT <- rep(NA, nreps)

cens_prop <- rep(NA, nreps)
rmse_bcart <- rmse_bart <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rep(NA, nreps)
rmse_aft <- rmse_aft_bart <- rmse_aft_null <- rmse_ipcw <- rep(NA, nreps)

coverage_bcart <- coverage_bart <- coverage_aft_bart <- rep(NA, nreps)

for(j in 1:nreps) {
  ## training set
  DataSim <- sim.reg(n, coef = coef, mu = mu, Rho = Rho)
  X.train <- DataSim$Z
  colnames(X.train) <- paste0('X', 1:num_covar)
  shape.train <- DataSim$Y*(1 + DataSim$Y)
  rate.train <- 1 + DataSim$Y
  T.train <- rgamma(n, shape=shape.train, rate=rate.train)
  mu.train <- DataSim$Y*pgamma(tau, shape = rate.train+1, rate = rate.train) + tau*pgamma(tau, shape = rate.train, 
                                                      rate = rate.train, lower.tail = FALSE)
  #C.train <-  runif(n, min=10, max=2000) ## max = 50 or 2000
  C.train <- rgamma(n, shape=2.2, rate=cens_rate)
  Y.train <- pmin(T.train, C.train)
  delta.train <- ifelse(T.train <= C.train, 1, 0) ## mean delta train 50-60 % or 80-90 %

  ## test set
  DataSim.test <- sim.reg(n.test, coef = coef, mu = mu, Rho = Rho)
  X.test <- DataSim.test$Z
  colnames(X.test) <- paste0('X', 1:num_covar)
  shape.test <- DataSim.test$Y*(1 + DataSim.test$Y)
  rate.test <- 1 + DataSim.test$Y
  T.test <- rgamma(n.test, shape=shape.test, rate=rate.test)
  mu.test <- DataSim.test$Y*pgamma(tau, shape = rate.test+1, rate = rate.test) + 
    tau*pgamma(tau, shape = rate.test, rate = rate.test, lower.tail = FALSE)
  #C.test <-  runif(n, min=10, max=2000) ## max = 50 or 2000
  C.test <- rgamma(n.test, shape=2.2, rate=cens_rate)
  Y.test <- pmin(T.test, C.test)
  delta.test <- ifelse(T.test <= C.test, 1, 0) ## mean delta train 50-60 % or 80-90 %

  tryCatch({
  ### 1. AFT linear model
  AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
  XX <- model.matrix(Y.test ~ X.test)
  aft_linpred <- as.numeric(XX%*%AFT$coefficients)
  aft_sigma <- AFT$scale
  aft_sigsq <- aft_sigma*aft_sigma
  ## exp(aft_linpred) is an approximate fitted value
  ## For RMST, a more precise definition of fitted value is
  
  gt_prob <- pnorm((log(tau) - aft_linpred)/aft_sigma, lower.tail=FALSE)
  lt_prob <- pnorm((log(tau) - aft_sigsq - aft_linpred)/aft_sigma)
  AFT_fitted <- exp(aft_sigsq/2 + aft_linpred)*lt_prob + tau*gt_prob
  ## do either bootstrap or derive confidence interval directly
  
  ##### 2. AFT intercept-only model
  AFT_null <- survreg(Surv(Y.train, delta.train) ~ 1)
  aft_linpred <- rep(AFT_null$coefficients[1], nrow(X.test))
  aft_sigma <- AFT_null$scale
  aft_sigsq <- aft_sigma*aft_sigma
  
  gt_prob <- pnorm((log(tau) - aft_linpred)/aft_sigma, lower.tail=FALSE)
  lt_prob <- pnorm((log(tau) - aft_sigsq - aft_linpred)/aft_sigma)
  AFT_null_fitted <- exp(aft_sigsq/2 + aft_linpred)*lt_prob + tau*gt_prob
  
  #### 3. AFT_BART model
  AFT_BART <- abart(X.train, Y.train, delta.train, x.test=X.test)
  ndraw_abart <- nrow(AFT_BART$yhat.test)
  AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(X.test))
  for(k in 1:ndraw_abart) {
    aft_bart_mu <- AFT_BART$yhat.test[k,]
    aft_bart_sig <- AFT_BART$sigma[k]
    aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
    ## exp(aft_linpred) is an approximate fitted value
    ## For RMST, a more precise definition of fitted value is
    
    gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
    lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)
    
    AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
  }
  AFT_BART_fitted <- colMeans(AFT_fit_reps)
  AFT_BART_CI <- t(apply(AFT_fit_reps, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
  
  
  ## 4. Boosting with IPCW weights
  ipcw_weights <- IPCweights(x=Surv(Y.train, delta.train))
  IPW_boost <- glmboost(x=X.train[delta.train==1,], y=pmin(Y.train[delta.train==1], tau),
                        weights = ipcw_weights[delta.train==1])
  IPW_fitted <- as.numeric(predict(IPW_boost, newdata=X.test))
  
  ## 5. Cox-PH model without penalization of regression coefficients
  COXPH.mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
  #coxhaz <- basehaz(COXPH.mod)
  coxhaz <- survfit(COXPH.mod, x=X.train, y=Surv(Y.train, delta.train))
  H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$cumhaz),
                    yright=max(coxhaz$cumhaz))
  
  COXPH_fitted <- CoxExpectedSurv(X=X.test, beta_val=COXPH.mod$coefficients,
                                  time = coxhaz$time, H0.vals=coxhaz$cumhaz,
                                  tau=tau)
  
  ### 6. Cox-PH model with lasso penalty
  rcox_tmp <- cv.glmnet(x=X.train, y=Surv(Y.train, delta.train), family = "cox",
                        type.measure = "deviance")
  RCOXPH <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox",
                    lambda = rcox_tmp$lambda.min)
  Rcoxhaz <- survfit(RCOXPH, x=X.train, y=Surv(Y.train, delta.train))
  RH0fn <- approxfun(c(0, Rcoxhaz$time), c(0, Rcoxhaz$cumhaz),
                     yright=max(Rcoxhaz$cumhaz))
  lasso_betahat <- as.numeric(coef(RCOXPH))
  RCOXPH_fitted <- CoxExpectedSurv(X=X.test, beta_val=lasso_betahat,
                                   time=Rcoxhaz$time,
                                   H0.vals=Rcoxhaz$cumhaz, tau=tau)
  
  ## 8. RMST BART
  U_tau <- pmin(Y.train[delta.train==1], tau)
  sgrid <- seq(0, tau, length.out=100)
  
  delta_alpha <- 1
  kappa0 <- 1
  Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
  for(k in 1:(ndraws + burnIn + 1)) {
    Gmat[k,] <- DrawIPCW(U=Y.train, delta=delta.train, Utau=U_tau, sgrid=sgrid,
                         kappa0=kappa0, delta_alpha=delta_alpha)
  }
  Gmat <- 1/sqrt(Gmat)
  
  bart_mod <- RMST_BART(Y.train, delta.train, X.train, Gweights=Gmat,
                        x.test=X.test, tau=tau, k = 2,
                        ndpost=ndraws, nskip=burnIn)
  bart_fitted <- bart_mod$yhat.test.mean
  BART_CI <- t(apply(bart_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
  
  ## RMST BCART
  bcart_mod <- RMST_BART(Y.train, delta.train, X.train, Gweights=Gmat,
                         x.test=X.test, tau=tau, k = 2, ntree=1L,
                         ndpost=ndraws, nskip=burnIn)
  bcart_fitted <- bcart_mod$yhat.test.mean
  BCART_CI <- t(apply(bcart_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
  
  ## Recording RMSE
  rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.test)*(bcart_fitted - mu.test)))
  rmse_bart[j] <- sqrt(mean((bart_fitted - mu.test)*(bart_fitted - mu.test)))
  rmse_coxph[j] <- sqrt(mean((COXPH_fitted - mu.test)*(COXPH_fitted - mu.test)))
  rmse_rcoxph[j] <- sqrt(mean((RCOXPH_fitted - mu.test)*(RCOXPH_fitted - mu.test)))
  #rmse_sboost[j] <- sqrt(mean((SBOOST_fitted - mu.test)*(SBOOST_fitted - mu.test)))
  rmse_ipcw[j] <- sqrt(mean((IPW_fitted - mu.test)*(IPW_fitted - mu.test)))
  rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.test)*(AFT_fitted - mu.test)))
  rmse_aft_bart[j] <- sqrt(mean((AFT_BART_fitted - mu.test)*(AFT_BART_fitted - mu.test)))
  rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.test)*(AFT_null_fitted - mu.test)))
  cens_prop[j] <- mean(1 - delta.train) # also record censoring proportion
  CorCT[j] <- cor(C.train, T.train)
  
  
  ## Recording coverage
  coverage_aft_bart[j] <- mean((mu.test >= AFT_BART_CI[,1]) & (mu.test <= AFT_BART_CI[,2]))
  coverage_bcart[j] <- mean((mu.test >= BCART_CI[,1]) & (mu.test <= BCART_CI[,2]))
  coverage_bart[j] <- mean((mu.test >= BART_CI[,1]) & (mu.test <= BART_CI[,2]))
  }, error = function(e) {
    # Handle the error gracefully (e.g., print a message)
    cat("Error occurred in iteration", j, ":", conditionMessage(e), "\n")
    # You can also assign default values or do nothing to skip the error
    rmse_bcart[j] <- NA  # Assigning a default value (e.g., NA)
    rmse_bart[j] <- NA
    rmse_coxph[j] <- NA
    rmse_rcoxph[j] <- NA
    # rmse_sboost[j] <- NA
    rmse_ipcw[j] <- NA
    rmse_aft[j] <- NA
    rmse_aft_bart[j] <- NA
    rmse_aft_null[j] <- NA
    cens_prop[j] <- NA
    CorCT[j] <- NA
    coverage_aft_bart[j] <- NA
    coverage_bcart[j] <- NA
    coverage_bart[j] <- NA
  })
  ## report the results of coverage as a table.
}


nmethods <- 8
Results <- matrix(NA, nrow=nmethods, ncol=2)
rownames(Results) <- c("AFT Null", "CoxPH", "Cox glmnet", "AFT linear",
                       "ipcw", "AFT BART", "BCART", "BART")
colnames(Results) <- c("Mean RMSE", "Median RMSE")

Results[1,1:2] <- c(mean(na.omit(rmse_aft_null)), median(na.omit(rmse_aft_null)))
Results[2,1:2] <- c(mean(na.omit(rmse_coxph)), mean(na.omit(rmse_coxph)))
Results[3,1:2] <- c(mean(na.omit(rmse_rcoxph)), median(na.omit(rmse_rcoxph)))
Results[4,1:2] <- c(mean(na.omit(rmse_aft)), median(na.omit(rmse_aft)))
Results[5,1:2] <- c(na.omit(mean(rmse_ipcw)), median(na.omit(rmse_ipcw)))
Results[6,1:2] <- c(mean(na.omit(rmse_aft_bart)), median(na.omit(rmse_aft_bart)))
Results[7,1:2] <- c(mean(na.omit(rmse_bcart)), median(na.omit(rmse_bcart)))
Results[8,1:2] <- c(mean(na.omit(rmse_bart)), median(na.omit(rmse_bart)))

round(Results, 4)
