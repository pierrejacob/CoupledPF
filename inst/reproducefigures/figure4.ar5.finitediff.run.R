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


nparticles <- 2^7
theta_star <- 0.3
nrep <- 1000
hs <- c(0.001, 0.025, 0.05)
nhs <- length(hs)

###
resampling_parameters <- list(epsilon = 0.05, desired_alpha = 0.99)
resampling_schemes <- c(CR_systematic, CR_hilbert, CR_indexmatching, CR_transport)
resampling_schemes_labels <- c("systematic", "Hilbert", "index-matching", "transport")
filename <- paste0("ar_correlation_resamplingcomparison.R", nrep, "D", dimension, "T", datalength,
                   "N", nparticles,  ".RData")
experiments <- data.frame()
for (ih in 1:nhs){
  cat("h : ", ih, "/", nhs, "\n")
  for (iresampling in 1:length(resampling_schemes)){
    cat("# resampling : ", iresampling, "/", length(resampling_schemes), "\n")
    theta1 <- theta_star
    theta2 <- theta1
    theta1 <- theta1 - hs[ih]
    theta2 <- theta2 + hs[ih]
    experiments_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
      randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
      res <- coupled_particle_filter(nparticles, ar, theta1, theta2, observations, randomness,
                                     resampling_schemes[[iresampling]], resampling_parameters)
      d <- data.frame(method = resampling_schemes_labels[iresampling], irep = irep, ll1 = res[1], ll2 = res[2])
      d
    }
    experiments_$h <- hs[ih]
    experiments <- rbind(experiments, experiments_)
  }
  save(experiments, hs, nhs, nparticles, nrep, datalength, observations, 
       resampling_schemes_labels,
       file = filename)
}
save(experiments, hs, nhs, nparticles, nrep, datalength, observations, 
     resampling_schemes_labels,
     file = filename)

