
library(splines)
library(Matrix)
library(lqmm)
library(MASS)
library(dplyr); library(tidyr)




###############################################################
################# BIC calculation Function ####################
###############################################################
BIC.bs = function(k,Y,W,time_interval,EF,wt){
  ### Y is response, W is observed matrix of W(t);
  ### k is the number of basis and should be at least 4;
  ### t is the number of time point data were observed
  cubic_bs = bs(time_interval, df=k, degree=3, intercept = TRUE)
  pred = t(t(W)-colMeans(W, na.rm = TRUE))%*%cubic_bs
  
  mod = lm(Y ~ pred + EF, weights = wt)
  
  re = BIC(mod)
  return(re)
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


################################################################
################# Data application Function ####################
################################################################

# df=nhanes_2005 %>% filter(paxhour >= 7 & paxhour <= 22)
# ID="seqn"
# Y="bmxbmi"
# Time="paxhour"
# FV="Weekday"
# IV="Weekend"
# EF=c("ridageyr","riagendr.Female","ridreth1_regroup.Black","ridreth1_regroup.Hispanic","ridreth1_regroup.Other","db.1")
# wt="wtint2yr"

linear_analysis = function(df,ID,Y,Time,FV,IV,EF,wt){
  
  ##########################################################################################
  ## Data Preparation 
  ##########################################################################################
  df = df %>% arrange(!!rlang::sym(ID), !!rlang::sym(Time))
  
  # Parameters
  t = length(unique(df[[Time]]))
  n = length(unique(df[[ID]]))
  a = seq(0, 1, length.out = t)
  
  # Outcome
  Y = df[!duplicated(df[[ID]]),] %>% pull(Y)
  
  # Functional variable with measurement error
  df.FV = df %>% dplyr::select(ID,Time,FV)
  W_t = df.FV %>% pivot_wider(id_cols = ID, names_from = Time, values_from = FV, names_prefix = "W.")
  W_t = as.matrix(W_t[,-1])
  
  # Instrumental variable
  df.IV = df %>% dplyr::select(ID,Time,IV)
  M_t = df.IV %>% pivot_wider(id_cols = ID, names_from = Time, values_from = IV, names_prefix = "M.")
  M_t = as.matrix(M_t[,-1])
  
  # Error-free covariates
  EF = as.matrix(df[!duplicated(df[[ID]]), EF])
  
  # Sample weight
  wt = df[!duplicated(df[[ID]]), wt]
  
  # Basis expansion
  # nbasis = 4+ceiling(n**0.1) ##based on sample size
  k_min = 5; k_max = 15 ##based on BIC
  bic = sapply(k_min:k_max, function(s) {BIC.bs(s, Y = Y, W = W_t, time_interval = a, EF = EF, wt = wt)})
  nbasis = (k_min:k_max)[which.min(bic)]
  bs2 = bs(a, df = nbasis, degree = 3, intercept = TRUE) ##cubic basis
  
  
  ##########################################################################################
  ## Data Analysis 
  ##########################################################################################
  #~~~~~~~~~~~~~~~#
  ##### Naive #####
  #~~~~~~~~~~~~~~~#
  W_i = t(t(W_t)-colMeans(W_t, na.rm = TRUE))%*%bs2
  
  mod = lm(Y ~ W_i + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_naive_intercept = coef[1]
  res_naive_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_naive_EF = coef[(nbasis+2):length(coef)]
  
  #~~~~~~~~~~~~~~~#
  ##### SIMEX #####
  #~~~~~~~~~~~~~~~#
  delta_t.est = colMeans(M_t, na.rm = TRUE)/colMeans(W_t, na.rm = TRUE) ##non-smoothed estimated delta(t) based on the assumption that M(t) = delta(t)*W(t)
  M_star = data.frame(M_t/matrix(rep(delta_t.est, n), nrow = n, ncol = t, byrow = TRUE))
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
      
      mod = lm(Y ~ W_lambda + EF, weights = wt)
      
      coef = summary(mod)$coefficients[,1]
      return(coef)
    })
  })
  gamma_simex.ave = Reduce("+", gamma.simex)/B ##average across B
  
  # mod.coefficients = as.vector(extra(p,gamma_simex.ave,lambda,"linear")) ##linear extrapolation
  mod.coefficients = as.vector(extra(nrow(gamma_simex.ave),gamma_simex.ave,lambda,"quadratic")) ##quadratic extrapolation
  
  res_simex_intercept = mod.coefficients[1]
  res_simex_W = crossprod(t(bs2),mod.coefficients[2:(nbasis+1)])
  res_simex_EF = mod.coefficients[(nbasis+2):length(mod.coefficients)]
  
  #~~~~~~~~~~~~#
  ##### LM #####
  #~~~~~~~~~~~~#
  W.lm = matrix(ncol = t, nrow = n)
  for (i in 1:t){
    fit = lm(W_t[,i] ~ M_t[,i])
    W.lm[,i] = predict(fit)
  }
  W.lm_i = t(t(W.lm)-colMeans(W.lm))%*%bs2 
  
  mod = lm(Y ~ W.lm_i + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_lm_intercept = coef[1]
  res_lm_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_lm_EF = coef[(nbasis+2):length(coef)]
  
  #~~~~~~~~~~~~~#
  ##### MLM #####
  #~~~~~~~~~~~~~#
  M_i = t(t(M_t)-colMeans(M_t))%*%bs2 ##dim n*nbasis
  
  W.mlm = matrix(ncol = nbasis, nrow = n)
  for (i in 1:nbasis){
    fit = lm(W_i[,i] ~ M_i)
    W.mlm[,i] = predict(fit)
  }
  
  mod = lm(Y ~ W.mlm + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_mlm_intercept = coef[1]
  res_mlm_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_mlm_EF = coef[(nbasis+2):length(coef)]
  
  
  # Beta Coefficient
  beta_intercept = data.frame(Naive = res_naive_intercept, SIMEX = res_simex_intercept, LM = res_lm_intercept, MLM = res_mlm_intercept)
  beta_W = data.frame(Time = c(7:22),
                      Naive = res_naive_W, SIMEX = res_simex_W, LM = res_lm_W, MLM = res_mlm_W)
  beta_EF = data.frame(Naive = res_naive_EF, SIMEX = res_simex_EF, LM = res_lm_EF, MLM = res_mlm_EF)
  
  
  res = list(a = a,
             beta_intercept = beta_intercept, 
             beta_W = beta_W, 
             beta_EF = beta_EF)
  return(res)
}


#########################################################
################# Bootstrap Function ####################
#########################################################

# df=nhanes_2005 %>% filter(paxhour >= 7 & paxhour <= 22)
# ID="seqn"
# Y="bmxbmi"
# Time="paxhour"
# FV="Weekday"
# IV="Weekend"
# EF=c("ridageyr","riagendr.Female","ridreth1_regroup.Black","ridreth1_regroup.Hispanic","ridreth1_regroup.Other","db.1")
# wt="wtint2yr"
# id_boot=id_boot_list[[1]]

linear_analysis_boot = function(df,ID,Y,Time,FV,IV,EF,wt,id_boot){
  
  ##########################################################################################
  ## Data Preparation 
  ##########################################################################################
  df = df %>% arrange(!!rlang::sym(ID),!!rlang::sym(Time))
  
  # Parameters
  t = length(unique(df[[Time]]))
  n = length(unique(df[[ID]]))
  a = seq(0, 1, length.out = t)
  
  df_boot = merge(df, id_boot, by = ID) %>%
    dplyr::group_by(.data[[Time]]) %>%  ## for same hour, rename the repeated subject
    mutate(!!ID := ave(as.character(.data[[ID]]), .data[[ID]], FUN = function(x)
      if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])) %>%
    ungroup() %>%
    as.data.frame() %>%
    arrange(!!rlang::sym(ID),!!rlang::sym(Time))
  
  # Outcome
  Y = df_boot[!duplicated(df_boot[[ID]]),] %>% pull(Y)
  
  # Functional variable with measurement error
  df.FV = df_boot %>% dplyr::select(ID,Time,FV)
  W_t = df.FV %>% pivot_wider(id_cols = ID, names_from = Time, values_from = FV, names_prefix = "W.")
  W_t = as.matrix(W_t[,-1])
  
  # Instrumental variable
  df.IV = df_boot %>% dplyr::select(ID,Time,IV)
  M_t = df.IV %>% pivot_wider(id_cols = ID, names_from = Time, values_from = IV, names_prefix = "M.")
  M_t = as.matrix(M_t[,-1])
  
  # Error-free covariates
  EF = as.matrix(df_boot[!duplicated(df_boot[[ID]]), EF])
  
  # Sample weight
  wt = df_boot[!duplicated(df_boot[[ID]]), wt]
  
  # Basis expansion
  # nbasis = 4+ceiling(n**0.1) ##based on sample size
  k_min = 5; k_max = 15 ##based on BIC
  bic = sapply(k_min:k_max, function(s) {BIC.bs(s, Y = Y, W = W_t, time_interval = a, EF = EF, wt = wt)})
  nbasis = (k_min:k_max)[which.min(bic)]
  bs2 = bs(a, df = nbasis, degree = 3, intercept = TRUE) ##cubic basis
  
  
  ##########################################################################################
  ## Data Analysis 
  ##########################################################################################
  #~~~~~~~~~~~~~~~#
  ##### Naive #####
  #~~~~~~~~~~~~~~~#
  W_i = t(t(W_t)-colMeans(W_t, na.rm = TRUE))%*%bs2
  
  mod = lm(Y ~ W_i + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_naive_intercept = coef[1]
  res_naive_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_naive_EF = coef[(nbasis+2):length(coef)]
  
  #~~~~~~~~~~~~~~~#
  ##### SIMEX #####
  #~~~~~~~~~~~~~~~#
  delta_t.est = colMeans(M_t, na.rm = TRUE)/colMeans(W_t, na.rm = TRUE) ##non-smoothed estimated delta(t) based on the assumption that M(t) = delta(t)*W(t)
  M_star = data.frame(M_t/matrix(rep(delta_t.est, n), nrow = n, ncol = t, byrow = TRUE))
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
      
      mod = lm(Y ~ W_lambda + EF, weights = wt)
      
      coef = summary(mod)$coefficients[,1]
      return(coef)
    })
  })
  gamma_simex.ave = Reduce("+", gamma.simex)/B ##average across B
  
  # mod.coefficients = as.vector(extra(p,gamma_simex.ave,lambda,"linear")) ##linear extrapolation
  mod.coefficients = as.vector(extra(nrow(gamma_simex.ave),gamma_simex.ave,lambda,"quadratic")) ##quadratic extrapolation
  
  res_simex_intercept = mod.coefficients[1]
  res_simex_W = crossprod(t(bs2),mod.coefficients[2:(nbasis+1)])
  res_simex_EF = mod.coefficients[(nbasis+2):length(mod.coefficients)]
  
  #~~~~~~~~~~~~#
  ##### LM #####
  #~~~~~~~~~~~~#
  W.lm = matrix(ncol = t, nrow = n)
  for (i in 1:t){
    fit = lm(W_t[,i] ~ M_t[,i])
    W.lm[,i] = predict(fit)
  }
  W.lm_i = t(t(W.lm)-colMeans(W.lm))%*%bs2 
  
  mod = lm(Y ~ W.lm_i + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_lm_intercept = coef[1]
  res_lm_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_lm_EF = coef[(nbasis+2):length(coef)]
  
  #~~~~~~~~~~~~~#
  ##### MLM #####
  #~~~~~~~~~~~~~#
  M_i = t(t(M_t)-colMeans(M_t))%*%bs2 ##dim n*nbasis
  
  W.mlm = matrix(ncol = nbasis, nrow = n)
  for (i in 1:nbasis){
    fit = lm(W_i[,i] ~ M_i)
    W.mlm[,i] = predict(fit)
  }
  
  mod = lm(Y ~ W.mlm + EF, weights = wt)
  
  coef = summary(mod)$coefficients[,1]
  res_mlm_intercept = coef[1]
  res_mlm_W = crossprod(t(bs2),coef[2:(nbasis+1)])
  res_mlm_EF = coef[(nbasis+2):length(coef)]
  
  
  # Beta Coefficient
  beta_intercept = data.frame(Naive = res_naive_intercept, SIMEX = res_simex_intercept, LM = res_lm_intercept, MLM = res_mlm_intercept)
  beta_W = data.frame(Time = c(7:22),
                      Naive = res_naive_W, SIMEX = res_simex_W, LM = res_lm_W, MLM = res_mlm_W)
  beta_EF = data.frame(Naive = res_naive_EF, SIMEX = res_simex_EF, LM = res_lm_EF, MLM = res_mlm_EF)
  
  
  res = list(a = a,
             id_boot = id_boot,
             df_boot = df_boot,
             beta_intercept = beta_intercept, 
             beta_W = beta_W, 
             beta_EF = beta_EF)
  return(res)
}



