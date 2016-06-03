#### This scripts briefly shows how this package works
# First, some initialization
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
#########

# This is where the action starts. First we create a synthetic dataset
# from an autoregressive model, in dimension 1, with 100 observations.

dimension <- 1
ar <- get_ar(dimension)

alpha_star <- 0.95
A_star <- create_A(alpha_star, dimension)
datalength <- 100
observations <- matrix(nrow = datalength, ncol = dimension)
x_t <- fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
for (time in 1:datalength){
  x_t <- t(A_star %*% t(x_t)) + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
  observations[time,] <- x_t + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
}

# now we estimate the log-likelihood with various methods.
# we use the following number of particles
nparticles <- 2^7 # For d = 1, T = 100
# number of experiments
nrep <- 5
# number of parameter values
nrhos <- 20
# grid of parameter values at which to estimate the log-likelihood
rhos <- seq(from = 0.3, to = 0.4, length.out = nrhos)
# estimate the log-likelihood using independent bootstrap particle filters
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

# evaluate the log-likelihood exactly using the Kalman filter
kfll <- foreach(rho = rhos, .combine = rbind) %dopar% {
  kalman_module <- Module( "kalman_mod", PACKAGE = "CoupledPF")
  LGModel <- new(kalman_module$LinearGaussian)
  LGModel$set_multivariate_parameters(rho, dimension)
  LGModel$set_observations(observations)
  Kalman <- new(kalman_module$Kalman)
  Kalman$setLinearGaussian(LGModel)
  Kalman$filtering()
  data.frame(rhos = rho, ll = Kalman$getLL())
} 

# estimate the log-likelihood with coupled particle filters
# with index-coupled resampling, here called index-matching.
resampling_scheme <- CR_indexmatching_given
resampling_parameters <- list()
indexll <- foreach (irep = 1:nrep, .combine = cbind) %dorng% {
  theta <- rhos[1]
  pf <- rep(0, nrhos)
  randomness <- ar$generate_randomness(nparticles = nparticles, datalength = datalength)
  particle_system <- particle_filter_storeall(nparticles, ar, theta, observations, randomness)
  pf[1] <- particle_system$ll
  for (irho in 2:nrhos){
    theta <- rhos[irho]
    particle_system <- coupled_pf_given(nparticles, ar, theta, observations, randomness,
                                        resampling_scheme, resampling_parameters, particle_system)
    pf[irho] <- particle_system$ll
  }
  pf
}

# now pool the results together for plotting
pfll.melt <- melt(cbind(rhos, data.frame(pfll)), id.vars = "rhos")
names(pfll.melt) <- c("rhos", "rep", "value") 
pfll.melt$method <- "indep."
index.melt <- melt(cbind(rhos, data.frame(indexll)), id.vars = "rhos")
names(index.melt) <- c("rhos", "rep", "value") 
index.melt$method <- "index"
all.melt <- rbind(pfll.melt, index.melt)
all.melt$method <- factor(all.melt$method, levels = c("indep.", "index"),
                          labels = c("independent", "CRN+index-coupled"))
# plotting
g <- ggplot(all.melt, aes(x = rhos, y = value, group = interaction(rep, method), colour = rep)) + geom_line() + geom_point()
g <- g + theme(legend.position = "none") + facet_wrap(~ method)
g <- g + scale_color_colorblind()
g <- g + xlab(expression(theta)) + ylab("log-likelihood")
g <- g + scale_x_continuous(breaks = c(0.3, 0.35, 0.4))
g <- g + geom_line(data = kfll, aes(x = rhos, y = ll, group = NULL, colour = NULL), colour = "red", size = 2)
g
# ta-daaaah!

