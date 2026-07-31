
library(foreach); library(doParallel)


# Setup parallel backend to use many processors
# cores = detectCores()
cl = makeCluster(8) #not to overload your computer cores[1]-1
registerDoParallel(cl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 1 - Sample Size Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result = 
  foreach(ss = c(100,500,1000,5000)) %dopar% {
    source("Simulation_LinearFuncIV.R")
  
    res_sum <- simulation_func(sim_iter=500, n=ss, t=100, me_dist="normal", 
                               sd_X=1.5, sd_U=1, sd_M=1, rho_X=0.5, rho_U=0.5, rho_M=0.5, vcov_X="AR1", vcov_U="AR1", vcov_M="AR1", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
  
    res_sum
  }

save(result, file = "Sample Size Effect.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 2 - Measurement Error Distribution Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result = 
  foreach(med = c("normal","t_dist","laplace")) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist=med, 
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=0.5, rho_W=0.5, rho_M=0.5, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "ME Distribution Effect.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 3 - Variance-Covariance Structure & Correlation for the Functional Error Terms Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Independent structure
source("Simulation_LinearFuncIV.R")

result = simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal",
                         sd_X=1.5, sd_W=1, sd_M=1, rho_X=0, rho_W=0, rho_M=0, vcov_X="IND", vcov_W="IND", vcov_M="IND", delta = 0.5,
                         gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                         seeds=123)

save(result, file = "VCOV Effect IND.RData")

# AR(1) structure
result = 
  foreach(r_X = c(0.25,0.5,0.75)) %:%
  foreach(r_W = c(0.25,0.5,0.75)) %:%
  foreach(r_M = c(0.25,0.5,0.75)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal",
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=r_X, rho_W=r_W, rho_M=r_M, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "VCOV Effect AR1.RData")

# Compound symmetry structure
result = 
  foreach(r_X = c(0.25,0.5,0.75)) %:%
  foreach(r_W = c(0.25,0.5,0.75)) %:%
  foreach(r_M = c(0.25,0.5,0.75)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal",
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=r_X, rho_W=r_W, rho_M=r_M, vcov_X="CS", vcov_W="CS", vcov_M="CS", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "VCOV Effect CS.RData")

# Unstructured
result = 
  foreach(r_X = c(0.25,0.5,0.75)) %:%
  foreach(r_W = c(0.25,0.5,0.75)) %:%
  foreach(r_M = c(0.25,0.5,0.75)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal",
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=r_X, rho_W=r_W, rho_M=r_M, vcov_X="UN", vcov_W="UN", vcov_M="UN", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "VCOV Effect UN.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 4 - Magnitudes of the Functional Predictor Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result = 
  foreach(sdx = c(1.0,1.5,2.0,4.0)) %:%
  foreach(sdw = c(0.5,1.0,2.0)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal",
                               sd_X=sdx, sd_W=sdw, sd_M=1, rho_X=0.5, rho_W=0.5, rho_M=0.5, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "FV SD Effect.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 5 - Magnitudes of the Instrumental Variable Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Standard deviation associated with the instrumental variable
result = 
  foreach(sdm = c(0.5,1.0,2.0,4.0)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal", 
                               sd_X=1.5, sd_W=1, sd_M=sdm, rho_X=0.5, rho_W=0.5, rho_M=0.5, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = 0.5, 
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "IV SD Effect.RData")

# δ(t) associated with the instrumental variable
result = 
  foreach(d = c(0,0.25,0.5,0.75)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal", 
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=0.5, rho_W=0.5, rho_M=0.5, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = d,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=0.1, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "IV Delta Effect.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation Study 6 - Magnitudes of the Random Error of Y Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result = 
  foreach(sdy = c(0.05,0.1,0.5,1,2)) %dopar% {
    source("Simulation_LinearFuncIV.R")
    
    res_sum <- simulation_func(sim_iter=500, n=1000, t=100, me_dist="normal", 
                               sd_X=1.5, sd_W=1, sd_M=1, rho_X=0.5, rho_W=0.5, rho_M=0.5, vcov_X="AR1", vcov_W="AR1", vcov_M="AR1", delta = 0.5,
                               gamma1=2, gamma2=0.6, sd_Z.c=0.5, prob_Z.b=0.6, sd_Y=sdy, 
                               seeds=123)
    
    res_sum
  }

save(result, file = "Y Effect.RData")


# Stop cluster
stopCluster(cl)



