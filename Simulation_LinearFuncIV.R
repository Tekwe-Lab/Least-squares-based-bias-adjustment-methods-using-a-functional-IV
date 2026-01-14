library(stats)
library(splines)
library(Matrix)
library(fda)
library(ald)
library(mvtnorm)
#library(boot)
library(SparseM)
library(ggplot2)
library(quantreg)
library(lqmm)
library(refund)
library(corpcor)
library(MASS)
library(popbio)
library(lme4)
library(plyr)
library(dplyr)
library(tidyr)
library(openxlsx)
library(gam)
library(mgcv)
library(reshape2)
library(mvnfast)
library(VGAM)
library(LaplacesDemon)
library(stringr)
library(gridExtra)




###############################################################
################# Data Generating Functions ###################
###############################################################
met = c("f1","f2","f3","f4","f5","f6","f7","f8")
Betafunc = function(t, met){
  switch(met,
         f1 = sin(2*pi*t),
         f2 = 1/(1+exp(4*2*(t-.5))), 
         f3 = sin(pi*(8*(t-.5))/2)/ (1 + (2*(8*(t-.5))^2)*(sign(8*(t-.5))+1)),
         f4 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(.5*sin(8*(t-.5))+1)),
         f5 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(sin(8*(t-.5))+1)),
         f6 = 1/(1+exp((t))),
         f7 = 1/(1+exp(4*2*(t-.5)))+1,
         f8 = 6/(1+exp(4*2*(t-.5)))+1
  )
}


###############################################################
################# BIC Calculation Function ####################
###############################################################
BIC.bs = function(k,Y,W,time_interval,Z.c,Z.b){
  ### Y is response, W is observed matrix of W(t);
  ### k is the number of basis and should be at least 4;
  ### t is the number of time point data were observed
  cubic_bs = bs(time_interval, df=k, degree=3, intercept = TRUE)
  pred = t(t(W)-colMeans(W))%*%cubic_bs
  
  mod = lm(Y ~ pred + Z.c + Z.b)
  
  re = BIC(mod)
  return(re)
}


####################################################################
################# Covariance Structure Function ####################
####################################################################
## This is a squared exponential function where the covariance depends on th distance between the points
## details see Jehav Functional paper (2021)
se_cov = function(sigma, t2, t1, l){
  ### sigma is the variance of ME
  ### s and t are different time points 
  ### l controls the covariance between different points
  ### larger value of l means larger range of dependence
  sigma*exp(-(t2-t1)^2/(2*l^2))
}

## This is a AR(1) correlation structure function
ar1_cov <- function(t, rho, sig){
  ### rho is the base correlation between time points
  ### t is the total number of time points 
  ### sig is the sigma of X(t)/W(t)
  exponent = abs(matrix(1:t-1, nrow = t, ncol = t, byrow = TRUE) - (1:t-1))
  (sig^2)*(rho^exponent)
}

## This is an unstructured correlation structure function
un_cov <- function(t, rho, sig){
  ### rho is the base correlation between time points
  ### t is the total number of time points 
  ### sig is the sigma of X(t)/W(t)
  cor.matrix = matrix(runif(t*t, max(rho-0.25,0), min(rho+0.25,1)), nrow = t, ncol = t)
  diag(cor.matrix) = 1 ##make the diagonal value to be 1
  make.positive.definite(as.matrix(forceSymmetric(sig^2*cor.matrix)))
}


######################################################################
################# Extrapolation function for SIMEX ###################
######################################################################
extra = function(k,beta_hat,lambda_seq,method){
  ########## linear ############
  if(method=="linear"){
    gamma = rep(0,k)
    for(m in 1:k)
    {
      testr<-beta_hat[m,]
      lr<-lm(testr~lambda_seq)
      coeff<-lr$coefficients
      gamma[m]<- coeff[1]+(-1)*coeff[2]
    }
    return(gamma)
  }
  
  ########### quadratic ##########
  else if(method=="quadratic"){
    gamma = rep(0,k)
    lambda2 = lambda_seq*lambda_seq
    for(m in 1:k)
    {
      testr<-beta_hat[m,]
      lr<-lm(testr~lambda_seq+lambda2)
      coeff<-lr$coefficients
      gamma[m]<- coeff[1]+(-1)*coeff[2]+(-1)*(-1)*coeff[3]
    }
    return(gamma)
  }
}


########################################################################################
################# Mean Squared Percentage Error Calculation Function ###################
########################################################################################
rmspe = function(beta, beta_hat, domain){
  f = (beta-beta_hat)^2
  len.f = length(f)
  result = sum((f[-1] + f[-len.f]) * diff(domain))/2
  f2 = beta^2
  len.f2 = length(f2)
  result2 = sum((f2[-1] + f2[-len.f2]) * diff(domain))/2
  return(100*sqrt(result/result2))
}


##########################################################
################# Simulation Function ####################
##########################################################
# sim_iter = 10       ## Number of simulations/iterations
#        n = 500      ## Sample size
#        t = 100      ## Number of timepoints for the functional data
#  me_dist = "normal" ## Distribution of measurement errors
#     sd_X = 3        ## Standard deviation of error_X(t)
#     sd_U = 2.5      ## Standard deviation of U(t)
#     sd_M = 1.5      ## Standard deviation of the instrumental variable 
#    rho_X = 0.5      ## Correlation coefficient
#    rho_U = 0.5      ## Correlation coefficient
#    rho_M = 0.5      ## Correlation coefficient
#   vcov_X = "CS"     ## Variance-covariance Matrix
#   vcov_U = "CS"     ## Variance-covariance Matrix
#   vcov_M = "CS"     ## Variance-covariance Matrix
#    delta = 0.5      ## Parameter for generating unknown delta(t) that quantifies the relationship between M(t) and X(t)
#   gamma1 = 2        ## Beta coefficient for continuous error-free covariate (scalar-valued)
#   gamma2 = 0.6      ## Beta coefficient for binary error-free covariate (scalar-valued)
#   sd_Z.c = 0.5      ## Standard deviation of the continuous EF covariate
# prob_Z.b = 0.6      ## Probability of the binary EF covariate
#     sd_Y = 0.1      ## Standard deviation of error for Y
#    seeds = 123
#     iter = 1

simulation_func = function(sim_iter, n, t, me_dist, 
                           sd_X, sd_U, sd_M, rho_X, rho_U, rho_M, vcov_X, vcov_U, vcov_M, delta,
                           gamma1, gamma2, sd_Z.c, prob_Z.b, sd_Y, 
                           seeds){
  
  set.seed(seeds)
  
  #----------------------------------------#
  ############# Parameter setup ############
  #----------------------------------------#
  a = seq(0, 1, length.out=t) ##create a sequence between 0,1 of length = t
  
  #### Covariance matrix for error terms
  ##variance-covariance matrix of X
  if (vcov_X == "CS"){
    Sigma_X = sd_X^2*(matrix(rho_X, nrow = t, ncol = t)+diag(rep(1-rho_X,t))) 
  } else if (vcov_X == "SE"){
    Sigma_X = sapply(a, function(i) se_cov(sd_X, a, i, 0.2))+diag(rep(sd_X^2-sd_X,t))
  } else if (vcov_X == "AR1"){
    Sigma_X = ar1_cov(t, rho_X, sd_X)
  } else if (vcov_X == "IND"){
    Sigma_X = diag(rep(sd_X^2, t))
  } else if (vcov_X == "UN"){
    Sigma_X = un_cov(t, rho_X, sd_X)
  }
  
  ##variance-covariance matrix of W measurement error
  if (vcov_U == "CS"){
    Sigma_U = sd_U^2*(matrix(rho_U, nrow = t, ncol = t)+diag(rep(1-rho_U,t))) 
  } else if (vcov_U == "SE"){
    Sigma_U = sapply(a, function(i) se_cov(sd_U, a, i, 0.2))+diag(rep(sd_U^2-sd_U,t))
  } else if (vcov_U == "AR1"){
    Sigma_U = ar1_cov(t, rho_U, sd_U)
  } else if (vcov_U == "IND"){
    Sigma_U = diag(rep(sd_U^2, t))
  } else if (vcov_U == "UN"){
    Sigma_U = un_cov(t, rho_U, sd_U)
  }
  
  ##variance-covariance matrix of M
  if (vcov_M == "CS"){
    Sigma_M = sd_M^2*(matrix(rho_M, nrow = t, ncol = t)+diag(rep(1-rho_M,t))) 
  } else if (vcov_M == "SE"){
    Sigma_M = sapply(a, function(i) se_cov(sd_M, a, i, 0.2))+diag(rep(sd_M^2-sd_M,t))
  } else if (vcov_M == "AR1"){
    Sigma_M = ar1_cov(t, rho_M, sd_M)
  } else if (vcov_M == "IND"){
    Sigma_M = diag(rep(sd_M^2, t))
  } else if (vcov_M == "UN"){
    Sigma_M = un_cov(t, rho_M, sd_M)
  }
  
  zero = rep(0,times=t)
  
  
  #--------------------------------------#
  ############# Result matrix ############
  #--------------------------------------#
  ##Oracle estimator
  beta_X = matrix(nrow=t,ncol=sim_iter)
  beta_EF.c = matrix(nrow=1,ncol=sim_iter)
  beta_EF.b = matrix(nrow=1,ncol=sim_iter)
  
  ##Naive estimator
  beta_naive_W = matrix(nrow=t,ncol=sim_iter) 
  beta_naive_EF.c = matrix(nrow=1,ncol=sim_iter)
  beta_naive_EF.b = matrix(nrow=1,ncol=sim_iter)
  
  ##SIMEX estimator
  beta_simex_W = matrix(nrow=t,ncol=sim_iter)
  beta_simex_EF.c = matrix(nrow=1,ncol=sim_iter)  
  beta_simex_EF.b = matrix(nrow=1,ncol=sim_iter) 
  
  ##PW-2SLS estimator
  beta_lm_W = matrix(nrow=t,ncol=sim_iter)
  beta_lm_EF.c = matrix(nrow=1,ncol=sim_iter)  
  beta_lm_EF.b = matrix(nrow=1,ncol=sim_iter) 
  
  ##MULTI-2SLS estimator
  beta_mlm_W = matrix(nrow=t,ncol=sim_iter)
  beta_mlm_EF.c = matrix(nrow=1,ncol=sim_iter)  
  beta_mlm_EF.b = matrix(nrow=1,ncol=sim_iter)
  
  ##number of basis selected by basis selection function
  selected_k_n=c()
  
  
  for (iter in 1:sim_iter){
    
    #------------------------------------------------------------------------------------------------------------------#
    #################################################### Simulation ####################################################
    #------------------------------------------------------------------------------------------------------------------#
    #### Simulate X(t) and W(t), n is the sample size
    err_X = mvrnorm(n, zero, Sigma_X) ##error term of X(t), assume the errors following normal distribution
    fx = "f7" ##function for generating X(t)
    X_t = matrix(rep(Betafunc(a,fx), n), nrow=n, byrow = TRUE) + err_X ##dim n*t
    
    if (me_dist == "normal"){
      U_t = mvrnorm(n, zero, Sigma_U) ##measurement error of X(t)
    } else if (me_dist == "t_dist"){
      U_t = mvnfast::rmvt(n, mu=zero, sigma=Sigma_U, df=n-1) ##measurement error of X(t)
    } else if (me_dist == "laplace"){
      U_t = rmvl(n, mu=zero, Sigma=Sigma_U) ##measurement error of X(t)
    }
    W_t = X_t + U_t
    
    #### Simulate M(t), n is the sample size
    delta_t = delta*Betafunc(a,"f1")+1
    err_M = mvrnorm(n, zero, Sigma_M) ##error term of M(t), following normal distribution
    M_t = t(apply(X_t, 1, function(s) delta_t*s)) + err_M ##dim n*t
    
    #### Simulate error free covariates
    Z.c = rnorm(n, 0, sd_Z.c)
    Z.b = rbinom(n, 1, prob_Z.b)
    
    #### Simulate outcome Y (with a normal distribution)
    fi = "f1"
    Y_norm = crossprod(t(X_t), Betafunc(a,fi))/t +            ##true association of X(t) and Y
      crossprod(t(Z.c), gamma1) + crossprod(t(Z.b), gamma2) + ##true association of error free covariates and Y
      rnorm(n, 0, sd=sd_Y)                                    ##error term of Y
    
    Y_norm_center = Y_norm - mean(Y_norm) ##center Y at mean
    
    #### Basis expansion of functional data
    # nbasis = ceiling(n^(1/4))+2 ##estimate based on sample size
    k_min = 5; k_max = 15 ##estimate based on BIC
    bic = sapply(k_min:k_max, function(s) {BIC.bs(s, Y = Y_norm_center, W = W_t, time_interval = a, Z.c = Z.c, Z.b = Z.b)})
    nbasis = (k_min:k_max)[which.min(bic)]
    selected_k_n[iter] = nbasis
    bs2 = bs(a, df = nbasis, degree = 3, intercept = TRUE) ##cubic basis
    
    
    #------------------------------------------------------------------------------------------------------------------------------------#
    ################################################## ME correction and Data analysis ###################################################
    #------------------------------------------------------------------------------------------------------------------------------------#
    #------------------------------------------#
    ################## Oracle ##################
    #------------------------------------------#
    X_i = t(t(X_t)-colMeans(X_t))%*%bs2 ##dim n*nbasis
    
    mod = lm(Y_norm_center ~ X_i + Z.c + Z.b) ##linear modeling
    
    beta_X[,iter] = crossprod(t(bs2),mod$coefficients[2:(nbasis+1)])*t
    beta_EF.c[,iter] = mod$coefficients[(nbasis+2)] 
    beta_EF.b[,iter] = mod$coefficients[(nbasis+3)] 
    
    #-----------------------------------------#
    ################## Naive ##################
    #-----------------------------------------#
    W_i = t(t(W_t)-colMeans(W_t))%*%bs2 ##dim n*nbasis
    
    mod = lm(Y_norm_center ~ W_i + Z.c + Z.b) ##linear modeling
    
    beta_naive_W[,iter] = crossprod(t(bs2),mod$coefficients[2:(nbasis+1)])*t
    beta_naive_EF.c[,iter] = mod$coefficients[(nbasis+2)] 
    beta_naive_EF.b[,iter] = mod$coefficients[(nbasis+3)]
    
    #-----------------------------------------#
    ################## SIMEX ##################
    #-----------------------------------------# 
    delta_t.est = colMeans(M_t)/colMeans(W_t) ##estimated delta 
    M_star = M_t/matrix(rep(delta_t.est, n), nrow = n, ncol = t, byrow = T)
    M_star_i = t(t(M_star)-colMeans(M_star))%*%bs2 
    
    Sigma_xx = cov(W_i, M_star_i, use="complete.obs") ##by assumption cov(W(t),M(t))/delta = Sigma_xx
    Sigma_ww = var(W_i, use="complete.obs") ##covariance matrix of observed surrogate var(W_ic) (centered)
    Sigma_uu = Sigma_ww - Sigma_xx ##covariance matrix of W|x from W=X+U
    Sigma_uu = make.positive.definite(as.matrix(forceSymmetric(Sigma_uu)))
    
    B = 100 ##number of replicates
    lambda = seq(0.0001,2.0001,.05) ##get a set of monotonically increasing small numbers
    gamma.simex = lapply(seq(1:B), function(b){
      sapply(lambda, function(s) {
        U_b = mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu, empirical = TRUE)
        W_lambda = W_i+(sqrt(s)*U_b)
        
        mod = lm(Y_norm_center ~ W_lambda + Z.c + Z.b) ##linear modeling
        return(mod$coefficients)
      })
    })
    gamma_simex.ave = Reduce("+", gamma.simex)/B ##average across B
    
    # mod.coefficients = as.vector(extra(p,gamma_simex.ave,lambda,"linear")) ##linear extrapolation
    mod.coefficients = as.vector(extra(nrow(gamma_simex.ave),gamma_simex.ave,lambda,"quadratic")) ##quadratic extrapolation
    # Note:
    # We tried to fit the nonlinear extrapolation methods and 
    # had a lot of computational issues with it and decided to just stick with the quadratic one.
    # Remove the nonlinear extrapolation for now.
    
    beta_simex_W[,iter] = crossprod(t(bs2),mod.coefficients[2:(nbasis+1)])*t
    beta_simex_EF.c[,iter] = mod.coefficients[(nbasis+2)] 
    beta_simex_EF.b[,iter] = mod.coefficients[(nbasis+3)]  
    
    #-------------------------------------------#
    ################## PW-2SLS ##################
    #-------------------------------------------#
    W.lm = matrix(ncol = t, nrow = n)
    for (i in 1:t){
      fit = lm(W_t[,i] ~ M_t[,i])
      W.lm[,i] = predict(fit)
    }
    W.lm_i = t(t(W.lm)-colMeans(W.lm))%*%bs2 ##dim n*nbasis
    
    mod = lm(Y_norm_center ~ W.lm_i + Z.c + Z.b) ##linear modeling
    
    beta_lm_W[,iter] = crossprod(t(bs2),mod$coefficients[2:(nbasis+1)])*t
    beta_lm_EF.c[,iter] = mod$coefficients[(nbasis+2)] 
    beta_lm_EF.b[,iter] = mod$coefficients[(nbasis+3)] 
    
    #----------------------------------------------#
    ################## MULTI-2SLS ##################
    #----------------------------------------------#
    M_i = t(t(M_t)-colMeans(M_t))%*%bs2 ##dim n*nbasis
    
    W.mlm = matrix(ncol = nbasis, nrow = n)
    for (i in 1:nbasis){
      fit = lm(W_i[,i] ~ M_i)
      W.mlm[,i] = predict(fit)
    }
    
    mod = lm(Y_norm_center ~ W.mlm + Z.c + Z.b) ##linear modeling
    
    beta_mlm_W[,iter] = crossprod(t(bs2),mod$coefficients[2:(nbasis+1)])*t
    beta_mlm_EF.c[,iter] = mod$coefficients[(nbasis+2)] 
    beta_mlm_EF.b[,iter] = mod$coefficients[(nbasis+3)]  
    
    
    # plot(a, sin(2*pi*a))
    # lines(a, beta_X[,iter], col="black")
    # lines(a, beta_naive_W[,iter], col="green")
    # lines(a, beta_simex_W[,iter], col="purple")
    # lines(a, beta_lm_W[,iter], col="orange")
    # lines(a, beta_mlm_W[,iter], col="red")
    
    
    # print(iter)
    
  } ### end of simulation loop
  
  re = list(beta_X    = beta_X,
            beta_EF.c = beta_EF.c,
            beta_EF.b = beta_EF.b,
          
            beta_naive_W    = beta_naive_W,
            beta_naive_EF.c = beta_naive_EF.c,
            beta_naive_EF.b = beta_naive_EF.b,
            
            beta_simex_W    = beta_simex_W,
            beta_simex_EF.c = beta_simex_EF.c,
            beta_simex_EF.b = beta_simex_EF.b,
            
            beta_lm_W    = beta_lm_W,
            beta_lm_EF.c = beta_lm_EF.c,
            beta_lm_EF.b = beta_lm_EF.b,
            
            beta_mlm_W    = beta_mlm_W,
            beta_mlm_EF.c = beta_mlm_EF.c,
            beta_mlm_EF.b = beta_mlm_EF.b,
            
            selected_k_n = selected_k_n)
  
  
  #--------------------------------------------------#
  ################## Result Summary ##################
  #--------------------------------------------------#
  m = t-1
  true_beta_X = Betafunc(a,fi)
  true_gamma1 = gamma1
  true_gamma2 = gamma2
  
  #### Get estimated coefficient
  ##Oracle
  hat_beta_X    = rowMeans(re$beta_X)
  hat_beta_EF.c = mean(re$beta_EF.c)
  hat_beta_EF.b = mean(re$beta_EF.b)
  mspe_beta_X    = rmspe(true_beta_X, hat_beta_X, a)
  bias_beta_X    = mean((hat_beta_X-true_beta_X)^2)
  bias_beta_EF.c = hat_beta_EF.c-true_gamma1
  bias_beta_EF.b = hat_beta_EF.b-true_gamma2
  var_beta_X    = sum((re$beta_X-hat_beta_X)^2)/(m*sim_iter)
  var_beta_EF.c = var(re$beta_EF.c[1,]) 
  var_beta_EF.b = var(re$beta_EF.b[1,]) 
  AIMSE_beta_X    = bias_beta_X + var_beta_X
  AIMSE_beta_EF.c = bias_beta_EF.c^2 + var_beta_EF.c
  AIMSE_beta_EF.b = bias_beta_EF.b^2 + var_beta_EF.b
  
  ##Naive
  hat_beta_naive_W    = rowMeans(re$beta_naive_W)
  hat_beta_naive_EF.c = mean(re$beta_naive_EF.c)
  hat_beta_naive_EF.b = mean(re$beta_naive_EF.b)
  mspe_beta_naive_W    = rmspe(true_beta_X, hat_beta_naive_W, a)
  bias_beta_naive_W    = mean((hat_beta_naive_W-true_beta_X)^2)
  bias_beta_naive_EF.c = hat_beta_naive_EF.c-true_gamma1
  bias_beta_naive_EF.b = hat_beta_naive_EF.b-true_gamma2
  var_beta_naive_W    = sum((re$beta_naive_W-hat_beta_naive_W)^2)/(m*sim_iter)
  var_beta_naive_EF.c = var(re$beta_naive_EF.c[1,])
  var_beta_naive_EF.b = var(re$beta_naive_EF.b[1,])
  AIMSE_beta_naive_W    = bias_beta_naive_W + var_beta_naive_W
  AIMSE_beta_naive_EF.c = bias_beta_naive_EF.c^2 + var_beta_naive_EF.c
  AIMSE_beta_naive_EF.b = bias_beta_naive_EF.b^2 + var_beta_naive_EF.b
  
  ##SIMEX
  hat_beta_simex_W    = rowMeans(re$beta_simex_W)
  hat_beta_simex_EF.c = mean(re$beta_simex_EF.c)
  hat_beta_simex_EF.b = mean(re$beta_simex_EF.b)
  mspe_beta_simex_W    = rmspe(true_beta_X, hat_beta_simex_W, a)
  bias_beta_simex_W    = mean((hat_beta_simex_W-true_beta_X)^2)
  bias_beta_simex_EF.c = hat_beta_simex_EF.c-true_gamma1
  bias_beta_simex_EF.b = hat_beta_simex_EF.b-true_gamma2
  var_beta_simex_W    = sum((re$beta_simex_W-hat_beta_simex_W)^2)/(m*sim_iter)
  var_beta_simex_EF.c = var(re$beta_simex_EF.c[1,])
  var_beta_simex_EF.b = var(re$beta_simex_EF.b[1,])
  AIMSE_beta_simex_W    = bias_beta_simex_W + var_beta_simex_W
  AIMSE_beta_simex_EF.c = bias_beta_simex_EF.c^2 + var_beta_simex_EF.c
  AIMSE_beta_simex_EF.b = bias_beta_simex_EF.b^2 + var_beta_simex_EF.b
  
  ##PW-2SLS
  hat_beta_lm_W    = rowMeans(re$beta_lm_W)
  hat_beta_lm_EF.c = mean(re$beta_lm_EF.c)
  hat_beta_lm_EF.b = mean(re$beta_lm_EF.b)
  mspe_beta_lm_W    = rmspe(true_beta_X, hat_beta_lm_W, a)
  bias_beta_lm_W    = mean((hat_beta_lm_W-true_beta_X)^2)
  bias_beta_lm_EF.c = hat_beta_lm_EF.c-true_gamma1
  bias_beta_lm_EF.b = hat_beta_lm_EF.b-true_gamma2
  var_beta_lm_W    = sum((re$beta_lm_W-hat_beta_lm_W)^2)/(m*sim_iter)
  var_beta_lm_EF.c = var(re$beta_lm_EF.c[1,])
  var_beta_lm_EF.b = var(re$beta_lm_EF.b[1,])
  AIMSE_beta_lm_W    = bias_beta_lm_W + var_beta_lm_W
  AIMSE_beta_lm_EF.c = bias_beta_lm_EF.c^2 + var_beta_lm_EF.c
  AIMSE_beta_lm_EF.b = bias_beta_lm_EF.b^2 + var_beta_lm_EF.b
  
  ##MULTI-2SLS
  hat_beta_mlm_W    = rowMeans(re$beta_mlm_W)
  hat_beta_mlm_EF.c = mean(re$beta_mlm_EF.c)
  hat_beta_mlm_EF.b = mean(re$beta_mlm_EF.b)
  mspe_beta_mlm_W    = rmspe(true_beta_X, hat_beta_mlm_W, a)
  bias_beta_mlm_W    = mean((hat_beta_mlm_W-true_beta_X)^2)
  bias_beta_mlm_EF.c = hat_beta_mlm_EF.c-true_gamma1
  bias_beta_mlm_EF.b = hat_beta_mlm_EF.b-true_gamma2
  var_beta_mlm_W    = sum((re$beta_mlm_W-hat_beta_mlm_W)^2)/(m*sim_iter)
  var_beta_mlm_EF.c = var(re$beta_mlm_EF.c[1,])
  var_beta_mlm_EF.b = var(re$beta_mlm_EF.b[1,])
  AIMSE_beta_mlm_W    = bias_beta_mlm_W + var_beta_mlm_W
  AIMSE_beta_mlm_EF.c = bias_beta_mlm_EF.c^2 + var_beta_mlm_EF.c
  AIMSE_beta_mlm_EF.b = bias_beta_mlm_EF.b^2 + var_beta_mlm_EF.b
  
  #### Combine results
  res_hat_betaX = data.frame(sample_size = n, timepoints = t, me_dist = me_dist,
                             sd_X = sd_X, sd_U = sd_U, sd_ratio_XW = sd_X/sd_U, sd_M = sd_M,
                             rho_X = rho_X, rho_U = rho_U, rho_M = rho_M, vcov_X = vcov_X, vcov_U = vcov_U, vcov_M = vcov_M, delta = delta,
                             gamma1 = gamma1, gamma2 = gamma2, sd_Z.c = sd_Z.c, prob_Z.b = prob_Z.b, sd_Y = sd_Y,
                             nbasis = mean(re$selected_k_n),
                             
                             true_beta_X = true_beta_X, 
                             hat_beta_X = hat_beta_X, 
                             hat_beta_mlm_W = hat_beta_mlm_W,
                             hat_beta_lm_W = hat_beta_lm_W, 
                             hat_beta_simex_W = hat_beta_simex_W,
                             hat_beta_naive_W = hat_beta_naive_W)
  
  res_hat_betaEF = data.frame(sample_size = n, timepoints = t, me_dist = me_dist,
                              sd_X = sd_X, sd_U = sd_U, sd_ratio_XW = sd_X/sd_U, sd_M = sd_M,
                              rho_X = rho_X, rho_U = rho_U, rho_M = rho_M, vcov_X = vcov_X, vcov_U = vcov_U, vcov_M = vcov_M, delta = delta,
                              gamma1 = gamma1, gamma2 = gamma2, sd_Z.c = sd_Z.c, prob_Z.b = prob_Z.b, sd_Y = sd_Y,
                              nbasis = mean(re$selected_k_n),
                              
                              true_beta_EF.c = true_gamma1, 
                              hat_beta_EF.c = hat_beta_EF.c,  
                              hat_beta_mlm_EF.c = hat_beta_mlm_EF.c, 
                              hat_beta_lm_EF.c = hat_beta_lm_EF.c, 
                              hat_beta_simex_EF.c = hat_beta_simex_EF.c, 
                              hat_beta_naive_EF.c = hat_beta_naive_EF.c, 
                              
                              true_beta_EF.b = true_gamma1, 
                              hat_beta_EF.b = hat_beta_EF.b, 
                              hat_beta_mlm_EF.b = hat_beta_mlm_EF.b, 
                              hat_beta_lm_EF.b = hat_beta_lm_EF.b, 
                              hat_beta_simex_EF.b = hat_beta_simex_EF.b, 
                              hat_beta_naive_EF.b = hat_beta_naive_EF.b)
  
  res_bias = data.frame(sample_size = n, timepoints = t, me_dist = me_dist, 
                        sd_X = sd_X, sd_U = sd_U, sd_ratio_XW = sd_X/sd_U, sd_M = sd_M,
                        rho_X = rho_X, rho_U = rho_U, rho_M = rho_M, vcov_X = vcov_X, vcov_U = vcov_U, vcov_M = vcov_M, delta = delta, 
                        gamma1 = gamma1, gamma2 = gamma2, sd_Z.c = sd_Z.c, prob_Z.b = prob_Z.b, sd_Y = sd_Y,
                        nbasis = mean(re$selected_k_n),
                        
                        mspe_beta_X = mspe_beta_X, 
                        mspe_beta_mlm_W = mspe_beta_mlm_W, 
                        mspe_beta_lm_W = mspe_beta_lm_W,
                        mspe_beta_simex_W = mspe_beta_simex_W, 
                        mspe_beta_naive_W = mspe_beta_naive_W, 
                        
                        bias_beta_X = bias_beta_X, 
                        bias_beta_mlm_W = bias_beta_mlm_W, 
                        bias_beta_lm_W = bias_beta_lm_W,
                        bias_beta_simex_W = bias_beta_simex_W, 
                        bias_beta_naive_W = bias_beta_naive_W, 
                        
                        bias_beta_EF.c = bias_beta_EF.c, 
                        bias_beta_mlm_EF.c = bias_beta_mlm_EF.c, 
                        bias_beta_lm_EF.c = bias_beta_lm_EF.c,
                        bias_beta_simex_EF.c = bias_beta_simex_EF.c, 
                        bias_beta_naive_EF.c = bias_beta_naive_EF.c, 
                        
                        bias_beta_EF.b = bias_beta_EF.b, 
                        bias_beta_mlm_EF.b = bias_beta_mlm_EF.b, 
                        bias_beta_lm_EF.b = bias_beta_lm_EF.b,
                        bias_beta_simex_EF.b = bias_beta_simex_EF.b, 
                        bias_beta_naive_EF.b = bias_beta_naive_EF.b)
  
  res_var = data.frame(sample_size = n, timepoints = t, me_dist = me_dist, 
                       sd_X = sd_X, sd_U = sd_U, sd_ratio_XW = sd_X/sd_U, sd_M = sd_M,
                       rho_X = rho_X, rho_U = rho_U, rho_M = rho_M, vcov_X = vcov_X, vcov_U = vcov_U, vcov_M = vcov_M, delta = delta, 
                       gamma1 = gamma1, gamma2 = gamma2, sd_Z.c = sd_Z.c, prob_Z.b = prob_Z.b, sd_Y = sd_Y,
                       nbasis = mean(re$selected_k_n),
                       
                       var_beta_X = var_beta_X, 
                       var_beta_mlm_W = var_beta_mlm_W, 
                       var_beta_lm_W = var_beta_lm_W,
                       var_beta_simex_W = var_beta_simex_W, 
                       var_beta_naive_W = var_beta_naive_W, 
                       
                       var_beta_EF.c = var_beta_EF.c, 
                       var_beta_mlm_EF.c = var_beta_mlm_EF.c, 
                       var_beta_lm_EF.c = var_beta_lm_EF.c,
                       var_beta_simex_EF.c = var_beta_simex_EF.c, 
                       var_beta_naive_EF.c = var_beta_naive_EF.c, 
                       
                       var_beta_EF.b = var_beta_EF.b, 
                       var_beta_mlm_EF.b = var_beta_mlm_EF.b, 
                       var_beta_lm_EF.b = var_beta_lm_EF.b,
                       var_beta_simex_EF.b = var_beta_simex_EF.b, 
                       var_beta_naive_EF.b = var_beta_naive_EF.b)
  
  res_AIMSE = data.frame(sample_size = n, timepoints = t, me_dist = me_dist, 
                         sd_X = sd_X, sd_U = sd_U, sd_ratio_XW = sd_X/sd_U, sd_M = sd_M,
                         rho_X = rho_X, rho_U = rho_U, rho_M = rho_M, vcov_X = vcov_X, vcov_U = vcov_U, vcov_M = vcov_M, delta = delta, 
                         gamma1 = gamma1, gamma2 = gamma2, sd_Z.c = sd_Z.c, prob_Z.b = prob_Z.b, sd_Y = sd_Y,
                         nbasis = mean(re$selected_k_n),
                         
                         AIMSE_beta_X = AIMSE_beta_X, 
                         AIMSE_beta_mlm_W = AIMSE_beta_mlm_W, 
                         AIMSE_beta_lm_W = AIMSE_beta_lm_W,
                         AIMSE_beta_simex_W = AIMSE_beta_simex_W, 
                         AIMSE_beta_naive_W = AIMSE_beta_naive_W, 
                         
                         AIMSE_beta_EF.c = AIMSE_beta_EF.c, 
                         AIMSE_beta_mlm_EF.c = AIMSE_beta_mlm_EF.c, 
                         AIMSE_beta_lm_EF.c = AIMSE_beta_lm_EF.c, 
                         AIMSE_beta_simex_EF.c = AIMSE_beta_simex_EF.c, 
                         AIMSE_beta_naive_EF.c = AIMSE_beta_naive_EF.c, 
                         
                         AIMSE_beta_EF.b = AIMSE_beta_EF.b,  
                         AIMSE_beta_mlm_EF.b = AIMSE_beta_mlm_EF.b, 
                         AIMSE_beta_lm_EF.b = AIMSE_beta_lm_EF.b,
                         AIMSE_beta_simex_EF.b = AIMSE_beta_simex_EF.b, 
                         AIMSE_beta_naive_EF.b = AIMSE_beta_naive_EF.b)
  
  
  res=list(re             = re,
           res_hat_betaX  = res_hat_betaX,
           res_hat_betaEF = res_hat_betaEF,
           res_bias       = res_bias,
           res_var        = res_var,
           res_AIMSE      = res_AIMSE)
  
  return(res)
  
} ### end of simulation function



