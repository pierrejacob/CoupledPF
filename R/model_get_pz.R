# PZ model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
# transformed so that all parameters are in R
# theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq)
#'@rdname get_pz
#'@title Phytoplankton-Zooplankton model
#'@description This function returns a list with objects such as
#'* rinit  to sample from the initial distribution
#'* rtransition to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* generate_randomness to evaluate the measurement density
#'* perturb_randomness to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list 
#'@export
get_pz <- function(){
  logit <- function(z) log(z / (1 - z))
  expit <- function(z) 1 / (1 + exp(-z))
  
  rinit <- function(nparticles, theta, rand, ...){
    return(exp(matrix(log(2) + rand[1:(2*nparticles)], nrow = 2)))
  }
  rinit_rand <- function(nparticles, theta){
    return(rnorm(2*nparticles))
  }
  #  
  rtransition <- function(xparticles, theta, time, rand, ...){
    ra <- rand[(2*nparticles + (time-1)*nparticles + 1):(2*nparticles+ (time)*nparticles)]
    alphas <- theta[1] + exp(theta[2]) * ra
    xparticles <- pz_transition(xparticles, alphas, time-1, expit(theta[3:6]))
    return(xparticles)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(dnorm(x = observation, mean = log(xparticles[1,]), sd = 0.2, log = TRUE))
  }
  #
  generate_randomness <- function(nparticles, datalength){
    return(pz_generate_randomness_cpp(nparticles, datalength))
  }
  #
  perturb_randomness <- function(randomness, rho){
    return(pz_perturb_randomness_cpp(randomness, rho))
  }
  #
  pz_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                   rtransition_rand = rtransition_rand,
                   dmeasurement = dmeasurement,
                   generate_randomness = generate_randomness, 
                   perturb_randomness = perturb_randomness,
                   dimension = 2)
  return(pz_model)
}
