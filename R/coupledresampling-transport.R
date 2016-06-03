#'@rdname CR_transport
#'@title Coupled Resampling: transport
#'@description This function performs coupled resampling based on transport
#' The parameters are given as a list with keys
#' * epsilon, which determines the approximation in the solution of the transport problem
#' * and desired_alpha, which determines implicitly the number of iterations performed by the iterative
#' algorithm.
#'@return Two vectors of ancestors, column-binded in a matrix.
#'@export
CR_transport <- function(xparticles1, xparticles2, normweights1, normweights2, 
                         parameters = list(epsilon = 0.1, desired_alpha = 0.9)){
  nparticles <- ncol(xparticles1)
  ancestors <- transport_cpp(xparticles1, xparticles2, normweights1, normweights2, runif(nparticles + 2), parameters$epsilon, parameters$desired_alpha)
  return(ancestors + 1)
}