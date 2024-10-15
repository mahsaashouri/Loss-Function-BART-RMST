library(survival)
library(mboost)
library(rmstbart)
library(AFTrees)
library(glmnet)
source("DrawIPCW.R")

set.seed(1234)
f.test <- function(x) {10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-.5)^2+10*x[ , 4]+5*x[ , 5]}

ndraws <- 1000
burnIn <- 500
n <- 5000   # 1000
n.test <- 5000
num_covar <- 10 # Choose to be 10 or 50
nreps <- 100 # number of simulation replications

theta <- 0.01
## Choose censoring rate to be either "heavy" or "medium"
cens_rate <- "heavy"
#cens_rate <- "medium"

tau <- 25

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

rmse_bcart <- rmse_bart <- rmse_bart_dep <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rep(NA, nreps)
rmse_aft <- rmse_aft_bart <- rmse_aft_null <- rmse_ipcw <- rep(NA, nreps)
rmse_bcart_default <- rmse_bart_default <- rmse_bart_dep_default <- rmse_aft_bart_default <- rep(NA, nreps)

bias_bart_dep <- bias_aft_bart <- bias_bcart <- bias_bart <- bias_aft_bart_default <- bias_bcart_default <- bias_bart_default <- rep(NA, nreps)
bias_coxph <- bias_rcoxph <- bias_ipcw <- bias_bart_dep_default <- bias_aft <- bias_aft_null <- rep(NA, nreps)

coverage_bcart <- coverage_bart <- coverage_dep_bart <- coverage_aft_bart <- rep(NA, nreps)
coverage_bcart_default <- coverage_bart_default <- coverage_dep_bart_default <- coverage_aft_bart_default <- rep(NA, nreps)
cens_prop <- rep(NA, nreps)
CorCT <- rep(NA, nreps)

for(j in 1:nreps) {
  tryCatch({
    X.train <- matrix(runif(n*num_covar), n, num_covar)
    colnames(X.train) <- paste0('X', 1:num_covar)
    ET.train <- f.test(X.train)
    
    shape.train <- ET.train*(1 + ET.train)
    rate.train <- 1 + ET.train
    T.train <- rgamma(n, shape=shape.train, rate=rate.train)
    
    mu.train <- ET.train*pgamma(tau, shape = rate.train+1, rate = rate.train) + tau*pgamma(tau, shape = rate.train, rate = rate.train, lower.tail = FALSE)
    
    if(cens_rate=="heavy") {
      C.train <- rexp(n, rate=theta*ET.train)
    } else {
      C.train <- rgamma(n, shape=3, rate=theta*ET.train)
    }
    delta.train <- ifelse(T.train < C.train, 1, 0)
    Y.train <- pmin(T.train, C.train)
    U_tau <- pmin(T.train[T.train < C.train], tau)
    
    X.test <- matrix(runif(n.test*num_covar), n.test, num_covar)
    colnames(X.test) <- paste0('X', 1:num_covar)
    ET.test <- f.test(X.test)
    
    shape.test <- ET.test*(1 + ET.test)
    rate.test <- 1 + ET.test
    mu.test <- ET.test*pgamma(tau, shape = rate.test+1, rate = rate.test) + tau*pgamma(tau, shape = rate.test, rate = rate.test, lower.tail = FALSE)
    
    T.test <- rgamma(n.test, shape=ET.test*(1 + ET.test), scale = 1/(1 + ET.train))
    
    AFT_try <- survreg(Surv(exp(Y.train), delta.train) ~ X.train)
    eta_hat <- AFT_try$scale*AFT_try$scale
    
    X.test <- matrix(runif(n.test*num_covar), n.test, num_covar)
    colnames(X.test) <- paste0('X', 1:num_covar)
    ET.test <- f.test(X.test)
    
    shape.test <- ET.test*(1 + ET.test)
    rate.test <- 1 + ET.test
    mu.test <- ET.test*pgamma(tau, shape = rate.test+1, rate = rate.test) + tau*pgamma(tau, shape = rate.test, rate = rate.test, lower.tail = FALSE)
    
    ### 1. AFT linear model
    AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
    XX <- model.matrix(rep(1, nrow(X.test)) ~ X.test)
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
    
    ################################################################
    ### 7. BART-RMST
    delta_alpha <- 1
    kappa0 <- 1
    sgrid <- seq(0, tau, by=.1)
    Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
    for(k in 1:(ndraws + burnIn + 1)) {
      Gmat[k,] <- DrawIPCW(U=Y.train, delta=delta.train, Utau=U_tau, sgrid=sgrid,
                           kappa0=kappa0, delta_alpha=delta_alpha)
    }
    Gmat_orig <- 1/sqrt(Gmat)
    Gmeans <- colMeans(1/Gmat)
    
    cens_bart <- AFTrees(x.train=X.train, y.train=Y.train, status=1-delta.train,
                         ndpost=ndraws + burnIn, verbose=FALSE, k=2)
    Mucens_draws <- cens_bart$m.train[,delta.train==1]
    GmatDep <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
    for(k in 1:(ndraws + burnIn)) {
      for(h in 1:length(U_tau)) {
        log.time.points <- log(U_tau[h])
        AA <- (log.time.points - cens_bart$locations[k,] - Mucens_draws[k,h])/cens_bart$sigma[k]
        Cprob <- sum(pnorm(AA, lower.tail=FALSE)*cens_bart$mix.prop[k,])
        GmatDep[k,h] <- 1/Cprob
      }
    }
    GmatDeporig <- 1/sqrt(GmatDep)
    
    
    ##########################################################################
    ## Doing Cross-Validation with RMST_BART and RMST_BCART to find best eta
    #########################################################################
    nfolds <- 5
    eta_hat_vals <- c(0.1*eta_hat, 0.25*eta_hat, 0.5*eta_hat, 0.75*eta_hat, eta_hat, 1.5*eta_hat)
    CVscore <- matrix(NA, nrow=nfolds, ncol=length(eta_hat_vals))
    X.train.obs <- X.train[delta.train==1,]
    Y.train.obs <- Y.train[delta.train==1]
    ntrain <- nrow(X.train.obs)
    fold_memb <- sample(1:nfolds, size=ntrain, replace=TRUE)
    CVscore <- CVscoreDep <- CVscore_BCART <- matrix(NA, nrow=nfolds, ncol=length(eta_hat_vals))
    cens_dist <- survfit(Surv(Y.train, 1-delta.train) ~ 1)
    GKM <- stepfun(c(0, cens_dist$time), c(1, cens_dist$surv, min(cens_dist$surv)))
    for(u in 1:length(eta_hat_vals)) {
      Gmat <- sqrt(2*eta_hat_vals[u])*Gmat_orig
      GmatDep <- sqrt(2*eta_hat_vals[u])*GmatDeporig
      for(k in 1:nfolds) {
        Y.train_tmp <- Y.train.obs[fold_memb!=k]
        delta.train_tmp <- rep(1, sum(fold_memb!=k))
        X.train_tmp <- X.train.obs[fold_memb!=k,]
        
        Y.test_tmp <- Y.train.obs[fold_memb==k]
        delta.test_tmp <- rep(1, sum(fold_memb==k))
        X.test_tmp <- X.train.obs[fold_memb==k,]
        
        Gmat_train_tmp <- Gmat[,fold_memb!=k]
        Gmat_test_tmp <- Gmat[,fold_memb==k]
        
        GmatDep_train_tmp <- GmatDep[,fold_memb!=k]
        GmatDep_test_tmp <- GmatDep[,fold_memb==k]
        
        ww_dep <- colMeans(1/(GmatDeporig[,fold_memb==k]^2))
        ww <- 1/GKM(Y.test_tmp)
        
        bartmod_tmp <- RMST_BART(Y.train_tmp, delta.train_tmp, X.train_tmp, Gweights=Gmat_train_tmp,
                                 x.test=X.test_tmp, tau=tau, k = 2, ndpost=ndraws, nskip=burnIn)
        
        bartmod_dep_tmp <- RMST_BART(Y.train_tmp, delta.train_tmp, X.train_tmp, Gweights=GmatDep_train_tmp,
                                     x.test=X.test_tmp, tau=tau, k = 2,
                                     ndpost=ndraws, nskip=burnIn)
        
        yhat <- bartmod_tmp$yhat.test.mean
        yhat_dep <- bartmod_dep_tmp$yhat.test.mean
        
        CVscore[k, u] <- mean(ww*((Y.test_tmp - yhat)*(Y.test_tmp - yhat)))
        CVscoreDep[k,u] <- mean(ww_dep*((Y.test_tmp - yhat_dep)*(Y.test_tmp - yhat_dep)))
        #CVscoreDep[k,u] <- mean(ww*((Y.test_tmp - yhat_dep)*(Y.test_tmp - yhat_dep)))
        
      }
    }
    CVfinal <- colMeans(CVscore)
    CVfinalDep <- colMeans(CVscoreDep)
    eta_hat_star <- eta_hat_vals[which.min(CVfinal)]
    eta_hat_star_dep <- eta_hat_vals[which.min(CVfinalDep)]
    
    #############################
    ## RMST-BART with best value of eta
    Gmat <- sqrt(2*eta_hat_star)*Gmat_orig
    bart_mod <- RMST_BART(Y.train, delta.train, X.train, Gweights=Gmat,
                          x.test=X.test, tau=tau, k = 2,
                          ndpost=ndraws, nskip=burnIn)
    bart_fitted <- bart_mod$yhat.test.mean
    BART_CI <- t(apply(bart_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
    
    
    #############################
    ## RMST-DEP-BART with best value of eta
    GmatDep <- sqrt(2*eta_hat_star_dep)*GmatDeporig
    bart_dep_mod <- RMST_BART(Y.train, delta.train, X.train, Gweights=GmatDep,
                              x.test=X.test, tau=tau, k = 2.0,
                              ndpost=ndraws, nskip=burnIn)
    bart_dep_fitted <- bart_dep_mod$yhat.test.mean
    BART_dep_CI <- t(apply(bart_dep_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
    
    
    ## RMST-BART default
    Gmat <- sqrt(2*eta_hat*0.5)*Gmat_orig
    bart_mod_default <- RMST_BART(Y.train, delta.train, X.train, Gweights=Gmat,
                                  x.test=X.test, tau=tau, k = 2,
                                  ndpost=ndraws, nskip=burnIn)
    bart_fitted_default <- bart_mod_default$yhat.test.mean
    
    
    ## RMST-Dep-BART default
    GmatDep <- sqrt(2*eta_hat*0.5)*GmatDeporig
    bart_dep_mod_default <- RMST_BART(Y.train, delta.train, X.train, Gweights=GmatDep,
                                      x.test=X.test, tau=tau, k = 2.0,
                                      ndpost=ndraws, nskip=burnIn)
    bart_dep_fitted_default <- bart_dep_mod_default$yhat.test.mean
    
    #### 3. AFT_BART model
    ## For AFT-BART cross validation, vary the sigquant and sigdf parameters
    ## as in the original BART paper
    sq_cand <- c(0.90, 0.99, 0.75)
    sdf_cand <- c(3, 3, 10)
    nfolds <- 5
    X.train.obs <- X.train[delta.train==1,]
    Y.train.obs <- Y.train[delta.train==1]
    ntrain <- nrow(X.train.obs)
    fold_memb <- sample(1:nfolds, size=ntrain, replace=TRUE)
    CVscore_AFT <- matrix(NA, nrow=nfolds, ncol=length(sq_cand))
    cens_dist <- survfit(Surv(Y.train, 1-delta.train) ~ 1)
    GKM <- stepfun(c(0, cens_dist$time), c(1, cens_dist$surv, min(cens_dist$surv)))
    for(u in 1:length(sq_cand)) {
      for(k in 1:nfolds) {
        Y.train_tmp <- Y.train.obs[fold_memb!=k]
        delta.train_tmp <- rep(1, sum(fold_memb!=k))
        X.train_tmp <- X.train.obs[fold_memb!=k,]
        
        Y.test_tmp <- Y.train.obs[fold_memb==k]
        delta.test_tmp <- rep(1, sum(fold_memb==k))
        X.test_tmp <- X.train.obs[fold_memb==k,]
        
        ww <- 1/GKM(Y.test_tmp)
        
        AFT_BART_tmp <- abart(X.train_tmp, Y.train_tmp, delta.train_tmp, x.test=X.test_tmp,
                              ndpost=ndraws, nskip=burnIn, sigdf=sdf_cand[u], sigquant=sq_cand[u])
        
        yhat <- AFT_BART_tmp$yhat.test.mean
        
        CVscore_AFT[k, u] <- mean(ww*((Y.test_tmp - yhat)*(Y.test_tmp - yhat)))
      }
    }
    CVfinal <- colMeans(CVscore_AFT)
    sq_star <- sq_cand[which.min(CVfinal)]
    sdf_star <- sdf_cand[which.min(CVfinal)]
    
    ## This is AFT BART with best hyperparameters:
    AFT_BART <- abart(X.train, Y.train, delta.train, x.test=X.test,
                      ndpost=ndraws, nskip=burnIn, sigdf=sdf_star, sigquant=sq_star)
    ndraw_abart <- nrow(AFT_BART$yhat.test)
    AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(X.test))
    for(k in 1:ndraw_abart) {
      aft_bart_mu <- AFT_BART$yhat.test[k,]
      aft_bart_sig <- AFT_BART$sigma[k]
      aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
      ## exp(aft_linpred) is an approximate fitted value
      ## For RMST, a more precise definition of fitted value is:
      
      gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
      lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)
      
      AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
    }
    AFT_BART_fitted <- colMeans(AFT_fit_reps)
    AFT_BART_CI <- t(apply(AFT_fit_reps, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
    
    ###################################################################
    ## Now, do AFT BART with default values of sigdf and sigquant
    AFT_BART <- abart(X.train, Y.train, delta.train, x.test=X.test,
                      ndpost=ndraws, nskip=burnIn)
    ndraw_abart <- nrow(AFT_BART$yhat.test)
    AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(X.test))
    for(k in 1:ndraw_abart) {
      aft_bart_mu <- AFT_BART$yhat.test[k,]
      aft_bart_sig <- AFT_BART$sigma[k]
      aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
      ## exp(aft_linpred) is an approximate fitted value
      ## For RMST, a more precise definition of fitted value is:
      
      gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
      lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)
      
      AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
    }
    AFT_BART_fitted_default <- colMeans(AFT_fit_reps)
    
    rmse_bart_dep[j] <- sqrt(mean((bart_dep_fitted - mu.test)*(bart_dep_fitted - mu.test)))
    rmse_bart[j] <- sqrt(mean((bart_fitted - mu.test)*(bart_fitted - mu.test)))
    #rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.test)*(bcart_fitted - mu.test)))
    rmse_bart_default[j] <- sqrt(mean((bart_fitted_default - mu.test)*(bart_fitted_default - mu.test)))
    rmse_bart_dep_default[j] <- sqrt(mean((bart_dep_fitted_default - mu.test)*(bart_dep_fitted_default - mu.test)))
    #rmse_bcart_default[j] <- sqrt(mean((bcart_fitted_default - mu.test)*(bcart_fitted_default - mu.test)))
    rmse_coxph[j] <- sqrt(mean((COXPH_fitted - mu.test)*(COXPH_fitted - mu.test)))
    rmse_rcoxph[j] <- sqrt(mean((RCOXPH_fitted - mu.test)*(RCOXPH_fitted - mu.test)))
    rmse_ipcw[j] <- sqrt(mean((IPW_fitted - mu.test)*(IPW_fitted - mu.test)))
    rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.test)*(AFT_fitted - mu.test)))
    rmse_aft_bart[j] <- sqrt(mean((AFT_BART_fitted - mu.test)*(AFT_BART_fitted - mu.test)))
    rmse_aft_bart_default[j] <- sqrt(mean((AFT_BART_fitted_default - mu.test)*(AFT_BART_fitted_default - mu.test)))
    rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.test)*(AFT_null_fitted - mu.test)))
    
    
    ## Recording mean of fitted values
    bias_bart_dep <- mean(bart_dep_fitted - mu.test)
    bias_bart_dep_default <- mean(bart_dep_fitted_default - mu.test)
    bias_aft_bart[j] <- mean(AFT_BART_fitted - mu.test)
    bias_bcart[j] <- mean(bcart_fitted - mu.test)
    bias_bart[j] <- mean(bart_fitted - mu.test)
    bias_aft_bart_default[j] <- mean(AFT_BART_fitted_default - mu.test)
    bias_bcart_default[j] <- mean(bcart_fitted_default - mu.test)
    bias_bart_default[j] <- mean(bart_fitted_default - mu.test)
    bias_coxph[j] <- mean(COXPH_fitted - mu.test)
    bias_rcoxph[j] <- mean(RCOXPH_fitted - mu.test)
    bias_ipcw[j] <- mean(IPW_fitted - mu.test)
    bias_aft[j] <- mean(AFT_fitted - mu.test)
    bias_aft_null[j] <- mean(AFT_null_fitted - mu.test)
    
    cens_prop[j] <- mean(delta.train)
    CorCT[j] <- cor(C.train, T.train)
  }, error=function(e){})
}
nmethods <- 13
Results <- matrix(NA, nrow=nmethods, ncol=2)
rownames(Results) <- c("AFT Null", "CoxPH", "Cox glmnet", "AFT linear",
                       "ipcw", "AFT BART", "AFT BART Default",
                       "BCART", "BCART Default", "BART", "BART Default",
                       "BART-DEP", "BART-DEP Default")
colnames(Results) <- c("Mean RMSE", "Median RMSE")

Results[1,1:2] <- c(mean(na.omit(rmse_aft_null)), median(na.omit(rmse_aft_null)))
Results[2,1:2] <- c(mean(na.omit(rmse_coxph)), mean(na.omit(rmse_coxph)))
Results[3,1:2] <- c(mean(na.omit(rmse_rcoxph)), median(na.omit(rmse_rcoxph)))
Results[4,1:2] <- c(mean(na.omit(rmse_aft)), median(na.omit(rmse_aft)))
Results[5,1:2] <- c(mean(na.omit(rmse_ipcw)), median(na.omit(rmse_ipcw)))
Results[6,1:2] <- c(mean(na.omit(rmse_aft_bart)), median(na.omit(rmse_aft_bart)))
Results[7,1:2] <- c(mean(na.omit(rmse_aft_bart_default)), median(na.omit(rmse_aft_bart_default)))
#Results[8,1:2] <- c(mean(rmse_bcart), median(rmse_bcart))
#Results[9,1:2] <- c(mean(rmse_bcart_default), median(rmse_bcart_default))
Results[10,1:2] <- c(mean(na.omit(rmse_bart)), median(na.omit(rmse_bart)))
Results[11,1:2] <- c(mean(na.omit(rmse_bart_default)), median(na.omit(rmse_bart_default)))
Results[12,1:2] <- c(mean(na.omit(rmse_bart_dep)), median(na.omit(rmse_bart_dep)))
Results[13,1:2] <- c(mean(na.omit(rmse_bart_dep_default)), median(na.omit(rmse_bart_dep_default)))

round(Results, 4)

## n=250-k=50-censoring'Medium': two rounds NA 

Bias <- matrix(NA, nrow = 1, ncol = 11)
colnames(Bias) <- c('BART-DEP', 'BART-DEP-default','AFT-BART', 'AFT-BART-default', 'BCART', 'BCART-default', 'BART', 'BART-default', 'coxph', 'rcoxph', 'ipcw', 'aft', 'aft-null')
Bias[,1] <- mean(bias_bart_dep)
Bias[,2] <- mean(bias_bart_dep_default)
Bias[,3] <- mean(bias_aft_bart)
Bias[,4] <- mean(bias_aft_bart_default)
Bias[,5] <- mean(bias_bcart)
Bias[,6] <- mean(bias_bcart_default)
Bias[,7] <- mean(bias_bart)
Bias[,8] <- mean(bias_bart_default)
Bias[,9] <- mean(bias_coxph)
Bias[,10] <- mean(bias_rcoxph)
Bias[,11] <- mean(bias_ipcw)
Bias[,12] <- mean(bias_aft, na.rm = T)
Bias[,13] <- mean(bias_aft_null)
