#'@rdname get_ar
#'@title Hidden auto-regressive model
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
get_ar <- function(dimension){
  #
  rinit <- function(nparticles, theta, rand, ...){
    return(ar_rinit_rcpp(nparticles, theta, rand, dimension))
  }
  rtransition <- function(xparticles, theta, time, rand, precomputed, ...){
    # return(precomputed$A %*% xparticles + rand[])
    return(ar_rtransition_rcpp(xparticles, theta, time, rand, dimension, precomputed$A))
    # return(precomputed$A %*% xparticles + 
    #          matrix(rand[(dimension * nparticles * time + 1):(dimension * nparticles * (time+1))],
    #                 ncol = ncol(xparticles)))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, precomputed, ...){
    return(fast_dmvnorm_transpose_cholesky(xparticles, observation, precomputed$di))
  }
  
  generate_randomness <- function(nparticles, datalength){
    return(ar_generate_randomness_cpp(nparticles, datalength, dimension))
  }
  #
  perturb_randomness <- function(randomness, rho){
    return(ar_perturb_randomness_cpp(randomness, rho, dimension))
  }
  #
  precompute <- function(theta){
    A <- create_A(theta, dimension)
    return(list(A = A, di = diag(1, dimension, dimension)))
  }
  
  #
  ar_model <- list(rinit = rinit, rtransition = rtransition, 
               dmeasurement = dmeasurement, generate_randomness = generate_randomness,
               perturb_randomness = perturb_randomness, precompute = precompute, dimension = dimension)
  
  return(ar_model)
}
