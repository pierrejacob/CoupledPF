# lorenz model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
# transformed so that all parameters are in R
# theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq)
#'@rdname get_lorenz
#'@title Lorenz 96 model
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
get_lorenz <- function(){
  rinit <- function(nparticles, theta, rand, ...){
    unif <- -1 + 4 * pnorm(rand[1:(8*nparticles)])
    return(matrix(unif, ncol = nparticles))
  }
  #  
  rtransition <- function(xparticles, theta, time, rand, ...){
    nparticles <- ncol(xparticles)
    ra <- rand[(8*nparticles + (time-1)*8*nparticles + 1):(8*nparticles + (time)*8*nparticles)]
    xparticles <- lorenz_transition(xparticles, time_start = 0.05*(time - 1), time_end = 0.05*time, dt = 0.05, theta[1])
    # print(dim(matrix(theta[2] * ra, ncol = nparticles)))
    # print(dim(xparticles))
    xparticles <- xparticles + matrix(theta[2] * ra, ncol = nparticles)
    return(xparticles)
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(fast_dmvnorm_transpose_cholesky(xparticles, observation, diag(0.5, 8, 8)))
  }
  #
  generate_randomness <- function(nparticles, datalength){
    return(lorenz_generate_randomness_cpp(nparticles, datalength))
  }
  #
  perturb_randomness <- function(randomness, rho){
    return(lorenz_perturb_randomness_cpp(randomness, rho))
  }
  #
  lorenz_model <- list(rinit = rinit, rtransition = rtransition,
                   dmeasurement = dmeasurement,
                   generate_randomness = generate_randomness, 
                   perturb_randomness = perturb_randomness,
                   dimension = 8)
  return(lorenz_model)
}
