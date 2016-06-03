# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(doRNG)
ncores <- 10
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

nrep <- 5

nparticles <- 2^6
mcmciterations <- 5000
# theta_init <- 0.4
theta_inits <- runif(nrep, 0.37, 0.41)
proposal_covariance <- matrix(0.001, 1, 1)

# theta_init <- 0.4
pmmh_parameters <- list(mcmciterations = mcmciterations, nparticles = nparticles, 
                        proposal_covariance = proposal_covariance, rho_perturb = 0.999)


# # standard pmmh
foreach(irep = 1:nrep, .combine = c) %dorng% {
  pmmh_res <- pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  filename <- paste0("pmmh.repeat", irep, "T", datalength, "N", nparticles, "M", mcmciterations, ".RData")
  save(pmmh_res, pmmh_parameters, file = filename)
  1.
}
# # 
# 
# # coupled pmmh - systematic resampling
pmmh_parameters$resampling_scheme <- CR_systematic_given
pmmh_parameters$resampling_parameters <- list()
foreach(irep = 1:nrep, .combine = c) %dorng% {
  pmmh_res <- coupled_pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  filename <- paste0("cpmmh.systematic.repeat", irep, "T", datalength, "N", nparticles, "M", pmmh_parameters$mcmciterations, ".RData")
  save(pmmh_res, pmmh_parameters, file = filename)
  1.
}
#
# # coupled pmmh - index-matching resampling
pmmh_parameters$resampling_scheme <- CR_indexmatching_given
foreach(irep = 1:nrep, .combine = c) %dorng% {
  pmmh_res <- coupled_pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  filename <- paste0("cpmmh.index.repeat", irep, "T", datalength, "N", nparticles, "M", pmmh_parameters$mcmciterations, ".RData")
  save(pmmh_res, pmmh_parameters, file = filename)
  1.
}

# # coupled pmmh - hilbert
pmmh_parameters$resampling_scheme <- CR_hilbert_given
foreach(irep = 1:nrep, .combine = c) %dorng% {
  # pmmh_res <- hilbert_pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  pmmh_res <- coupled_pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  filename <- paste0("cpmmh.hilbert.repeat", irep, "T", datalength, "N", nparticles, "M", pmmh_parameters$mcmciterations, ".RData")
  save(pmmh_res, pmmh_parameters, file = filename)
  1.
}

# coupled pmmh - transport resampling
pmmh_parameters$resampling_scheme <- CR_transport_given
pmmh_parameters$resampling_parameters <- list(epsilon = 0.1, desired_alpha = 0.95)

foreach(irep = 1:nrep, .combine = c) %dorng% {
  pmmh_res <- coupled_pmmh(pmmh_parameters, ar, theta_inits[irep], observations)
  filename <- paste0("cpmmh.transport.repeat", irep, "T", datalength, "N", nparticles, "M", pmmh_parameters$mcmciterations, ".RData")
  save(pmmh_res, pmmh_parameters, file = filename)
  1.
}

