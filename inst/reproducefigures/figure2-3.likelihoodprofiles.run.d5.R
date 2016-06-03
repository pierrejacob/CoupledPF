# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(doRNG)
library(ggthemes)
ncores <- 8
registerDoMC(cores = ncores)
# load custom theme for the plots
setmytheme()
# fix the random seed
set.seed(17)


dimension <- 5
# dimension <- 5
ar <- get_ar(dimension)

load(paste0("ar", dimension, "data.RData"))
# datalength <- 1000
datalength <- 1000

observations <- matrix(observations[1:datalength,], ncol = dimension)

theta <- 0.4

# nparticles <- 2^10 # For d = 5, T = 1000
nparticles <- 2^7 # For d = 1, T = 100
nrep <- 5
nrhos <- 20
rhos <- seq(from = 0.3, to = 0.4, length.out = nrhos)


##
pfll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  pf <- rep(0, nrhos)
  for (irho in 1:nrhos){
    theta <- rhos[irho]
    randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
    ll <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)$ll
    pf[irho] <- ll
  }
  pf
}


kfll <- foreach(rho = rhos, .combine = rbind) %dopar% {
  kalman_module <- Module( "kalman_mod", PACKAGE = "CoupledPF")
  LGModel <- new(kalman_module$LinearGaussian)
  LGModel$set_multivariate_parameters(rho, dimension)
  LGModel$set_observations(observations)
  Kalman <- new(kalman_module$Kalman)
  Kalman$setLinearGaussian(LGModel)
  Kalman$filtering()
  data.frame(rhos = rho, ll = Kalman$getLL())
  # ll <- kalman_loglikelihood(list(rho = rho, sigma = 1, eta = 1, tau = 1), observations)
  # data.frame(rhos = rho, ll = ll)
} 

# pfll <- cbind(rhos, pfll)
# pfll <- data.frame(pfll)
# pfll.melt <- melt(pfll, "rhos")
# g <- ggplot(df, aes(x = rhos, y = value, group = variable)) + geom_line()
# g <- g + theme(legend.position = "none")
# g <- g + xlab(expression(rho)) + ylab("log-likelihood")
# g <- g + geom_line(aes(y = ll), colour = "red")
# print(g)

### Common random numbers
common_pfll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  pf <- rep(0, nrhos)
  randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
  for (irho in 1:nrhos){
    theta[1] <- rhos[irho]
    ll <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)$ll
    pf[irho] <- ll
  }
  pf
}

### Index-matching resampling
resampling_scheme <- CR_indexmatching_given
resampling_parameters <- list()

indexll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  pf <- rep(0, nrhos)
  randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
  particle_system <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)
  pf[1] <- particle_system$ll
  for (irho in 1:nrhos){
    theta <- rhos[irho]
    particle_system <- coupled_pf_given(nparticles, ar, theta, observations, randomness,
                                        resampling_scheme, resampling_parameters, particle_system)
    pf[irho] <- particle_system$ll
  }
  pf
}


### transport resampling
resampling_scheme <- CR_hilbert_given
resampling_parameters <- list(epsilon = 0.01, desired_alpha = 0.99)

hilbertll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  pf <- rep(0, nrhos)
  randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
  particle_system <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)
  pf[1] <- particle_system$ll
  for (irho in 1:nrhos){
    theta <- rhos[irho]
    particle_system <- coupled_pf_given(nparticles, ar, theta, observations, randomness,
                                        resampling_scheme, resampling_parameters, particle_system)
    pf[irho] <- particle_system$ll
  }
  pf
}

### transport resampling
resampling_scheme <- CR_transport_given
resampling_parameters <- list(epsilon = 0.01, desired_alpha = 0.99)

transportll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  pf <- rep(0, nrhos)
  randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
  particle_system <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)
  pf[1] <- particle_system$ll
  for (irho in 1:nrhos){
    theta <- rhos[irho]
    particle_system <- coupled_pf_given(nparticles, ar, theta, observations, randomness,
                                        resampling_scheme, resampling_parameters, particle_system)
    pf[irho] <- particle_system$ll
  }
  pf
}


pfll.melt <- melt(cbind(rhos, data.frame(pfll)), id.vars = "rhos")
names(pfll.melt) <- c("rhos", "rep", "value") 

common.melt <- melt(cbind(rhos, data.frame(common_pfll)), id.vars = "rhos")
names(common.melt) <- c("rhos", "rep", "value") 
common.melt$method <- "common"

index.melt <- melt(cbind(rhos, data.frame(indexll)), id.vars = "rhos")
names(index.melt) <- c("rhos", "rep", "value") 
index.melt$method <- "index"

hilbert.melt <- melt(cbind(rhos, data.frame(hilbertll)), id.vars = "rhos")
names(hilbert.melt) <- c("rhos", "rep", "value") 
hilbert.melt$method <- "hilbert"

transport.melt <- melt(cbind(rhos, data.frame(transportll)), id.vars = "rhos")
names(transport.melt) <- c("rhos", "rep", "value") 
transport.melt$method <- "transport"

pfll.melt$method <- "indep."

all.melt <- rbind(pfll.melt, common.melt, index.melt, hilbert.melt, transport.melt)
all.melt$method <- factor(all.melt$method, levels = c("indep.", "common", "transport", "index", "hilbert"),
                          labels = c("independent", "CR+systematic", "CR+transport", "CR+index-matching", "CR+sorted"))

filename <- paste0("likelihoodprofiles.D", dimension, "T", datalength, "N", nparticles, ".RData")
save(kfll, all.melt, rhos, datalength, nparticles, dimension, file = filename)

