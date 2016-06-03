#'@rdname coupled_particle_filter
#'@title coupled_particle_filter
#'@description coupled_particle_filter
#'@export
## this couples two filters by feeding the same seed + doing transport during resampling
coupled_particle_filter <- function(nparticles, model, theta1, theta2, observations, randomness, 
                                    coupled_resampling, resampling_parameters = list()){
  datalength <- nrow(observations)
  
  precomputed1 <- try(model$precompute(theta1))
  if (inherits(precomputed1, "try-error")){
    precomputed1 <- NULL
  }
  precomputed2 <- try(model$precompute(theta2))
  if (inherits(precomputed2, "try-error")){
    precomputed2 <- NULL
  }
  
  # initialization
  xparticles1 <- model$rinit(nparticles, theta1, randomness)
  normweights1 <- rep(1/nparticles, nparticles)
  ll1 <- 0
  #
  xparticles2 <- model$rinit(nparticles, theta2, randomness)
  normweights2 <- rep(1/nparticles, nparticles)
  ll2 <- 0
  # step t > 1
  for (time in 1:datalength){
    ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2, resampling_parameters)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    xparticles1 <- xparticles1[,ancestors1]
    xparticles2 <- xparticles2[,ancestors2]
    #
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = model$dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = model$dimension)
    #
    xparticles1 <- model$rtransition(xparticles1, theta1, time, randomness, precomputed1)
    logw1 <- model$dmeasurement(xparticles1, theta1, observations[time,], precomputed1)
    #
    xparticles2 <- model$rtransition(xparticles2, theta2, time, randomness, precomputed2)
    logw2 <- model$dmeasurement(xparticles2, theta2, observations[time,], precomputed2)
    # book-keeping
    maxlw1 <- max(logw1)
    w1 <- exp(logw1 - maxlw1)
    ll1 <- ll1 + maxlw1 + log(mean(w1))
    normweights1 <- w1 / sum(w1)
    #
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    ll2 <- ll2 + maxlw2 + log(mean(w2))
    normweights2 <- w2 / sum(w2)
  }
  return(c(ll1, ll2))
}
