
library(foreach); library(doParallel)


# Setup parallel backend to use many processors
# cores = detectCores()
cl = makeCluster(8) #not to overload your computer cores[1]-1
registerDoParallel(cl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample Size Effect
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


save(result, file = "Output/Sample Size Effect.RData")


# Stop cluster
stopCluster(cl)