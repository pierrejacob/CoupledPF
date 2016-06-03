# pmmh_parameters should contain mcmciterations and proposalcovariance
#'@export
kalman_mh <- function(pmmh_parameters, model, theta_init, dimension, observations,
                      kalman_module){
  kalman_ll <- function(theta){
    LGModel <- new(kalman_module$LinearGaussian)
    LGModel$set_multivariate_parameters(theta, dimension)
    LGModel$set_observations(observations)
    Kalman <- new(kalman_module$Kalman)
    Kalman$setLinearGaussian(LGModel)
    Kalman$filtering()
    return(Kalman$getLL())
  }
  current_theta <- theta_init
  theta_dim <- length(theta_init)
  mcmciterations <- pmmh_parameters$mcmciterations
  proposal_covariance <- pmmh_parameters$proposal_covariance
  #
  current_ll <- kalman_ll(current_theta)
  current_posterior <- current_ll + model$dprior(current_theta)
  pmmh_naccepts <- 0
  pmmh_chain <- matrix(nrow = mcmciterations, ncol = theta_dim)
  pmmh_chain[1,] <- current_theta
  loglikelihoods <- rep(0, mcmciterations)
  logposteriors <- rep(0, mcmciterations)
  loglikelihoods[1] <- current_ll
  logposteriors[1] <- current_posterior
  for (iteration in 2:mcmciterations){
    if (iteration %% 100 == 1){
      cat("iteration: ", iteration, " / ", mcmciterations, "\n")
      cat("acceptance rate: ", pmmh_naccepts / iteration * 100, "%\n")
    }
    proposal <- current_theta + fast_rmvnorm(1, rep(0, theta_dim), proposal_covariance)[1,]
    proposal_prior <- model$dprior(proposal)
    # first test whether prior density is > 0 (otherwise reject)
    if (!is.infinite(proposal_prior)){
      proposal_ll <- kalman_ll(proposal)
      # if error in computing the log-likelihood, set it to -Infinity
      if (inherits(proposal_ll, "try-error") || is.na(proposal_ll)){
        proposal_ll <- -Inf
      }
      proposal_posterior <- proposal_ll + proposal_prior
      if (log(runif(1)) < (proposal_posterior - current_posterior)){
        current_theta <- proposal
        current_ll <- proposal_ll
        current_posterior <- proposal_posterior
        pmmh_naccepts <- pmmh_naccepts + 1
      }
    }
    pmmh_chain[iteration,] <- current_theta
    loglikelihoods[iteration] <- current_ll
    logposteriors[iteration] <- current_posterior
  }
  return(list(chain = pmmh_chain, acceptance_rate = pmmh_naccepts / mcmciterations,
              loglikelihoods = loglikelihoods, logposteriors = logposteriors))
}