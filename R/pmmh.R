# pmmh_parameters should contain nparticles, mcmciterations and proposalcovariance
#'@export
pmmh <- function(pmmh_parameters, model, theta_init, observations){
  current_theta <- theta_init
  theta_dim <- length(theta_init)
  nparticles <- pmmh_parameters$nparticles
  mcmciterations <- pmmh_parameters$mcmciterations
  proposal_covariance <- pmmh_parameters$proposal_covariance
  datalength <- nrow(observations)
  #
  randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
  current_pf <- particle_filter_storeall(nparticles, model, current_theta, observations, randomness)
  current_ll <- current_pf$ll
  current_posterior <- current_ll + model$dprior(current_theta)
  pmmh_naccepts <- 0
  pmmh_chain <- matrix(nrow = mcmciterations, ncol = theta_dim)
  pmmh_chain[1,] <- current_theta
  proposals <- matrix(NA, nrow = mcmciterations, ncol = theta_dim)
  loglikelihoods <- rep(0, mcmciterations)
  proposal_loglikelihoods <- rep(NA, mcmciterations)
  proposal_priors <- rep(NA, mcmciterations)
  logposteriors <- rep(0, mcmciterations)
  loglikelihoods[1] <- current_ll
  logposteriors[1] <- current_posterior
  for (iteration in 2:mcmciterations){
    if (iteration %% 100 == 1){
      cat("iteration: ", iteration, " / ", mcmciterations, "\n")
      cat("acceptance rate: ", pmmh_naccepts / iteration * 100, "%\n")
    }
    proposal <- current_theta + fast_rmvnorm(1, rep(0, theta_dim), proposal_covariance)[1,]
    proposals[iteration,] <- proposal
    proposal_prior <- model$dprior(proposal)
    # first test whether prior density is > 0 (otherwise reject)
    if (!is.infinite(proposal_prior)){
      randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
      proposal_pf <- try(particle_filter_storeall(nparticles, model, proposal, observations, randomness))
      # if error in computing the log-likelihood, set it to -Infinity
      if (inherits(proposal_pf, "try-error") || is.na(proposal_pf$ll)){
        proposal_ll <- -Inf
      } else {
        proposal_ll <- proposal_pf$ll
      }
      proposal_loglikelihoods[iteration] <- proposal_ll
      proposal_priors[iteration] <- proposal_prior
      proposal_posterior <- proposal_ll + proposal_prior
      if (log(runif(1)) < (proposal_posterior - current_posterior)){
        current_theta <- proposal
        current_ll <- proposal_ll
        current_posterior <- proposal_posterior
        current_pf <- proposal_pf
        pmmh_naccepts <- pmmh_naccepts + 1
      }
    }
    pmmh_chain[iteration,] <- current_theta
    loglikelihoods[iteration] <- current_ll
    logposteriors[iteration] <- current_posterior
  }
  return(list(chain = pmmh_chain, acceptance_rate = pmmh_naccepts / mcmciterations,
              loglikelihoods = loglikelihoods))
}