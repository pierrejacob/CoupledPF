# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(doRNG)
ncores <- 2
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)

dimension <- 5
ar <- get_ar(dimension)

load(paste0("ar", dimension, "data.RData"))
datalength <- 1000
observations <- matrix(observations[1:datalength,], ncol = dimension)


### prior
ar$dprior <- function(theta, ...){
  return(dnorm(theta, mean = 0, sd = 1, log = TRUE))
}
###


dosave <- TRUE

mcmciterations <- 100000
theta_init <- 0.4
proposal_covariance <- matrix(0.0001, 1, 1)
pmmh_parameters <- list(mcmciterations = mcmciterations, proposal_covariance = proposal_covariance)

# KF
kalman_module <- Module( "kalman_mod", PACKAGE = "CoupledPF")
kf_res <- kalman_mh(pmmh_parameters, ar, theta_init, dimension, observations, kalman_module)
filename <- paste0("ar.kfmh.D", dimension, "T", datalength, ".RData")
if (dosave){
  save(kf_res, file = filename)
}


