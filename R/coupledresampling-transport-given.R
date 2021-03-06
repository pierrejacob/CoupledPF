#'@rdname CR_transport_given
#'@title Coupled Resampling: transport, given vector of ancestor
#'@description This function performs tranasprt resampling on the second system
#'@return A vector of ancestors
#'@export
CR_transport_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
                               parameters = list(epsilon = 0.1, desired_alpha = 0.9),
                               ancestors_ref, uniforms){
  
  return(transport_given_cpp(xparticles1, xparticles2, normweights1, normweights2, uniforms, parameters$epsilon, parameters$desired_alpha,
                             ancestors_ref - 1) + 1)
}




# CR_transport_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
#                                            parameters = list(epsilon = 0.1, desired_alpha = 0.9),
#                                            ancestors_ref){
#   nparticles <- nrow(xparticles1)
#   M <- cost_matrix(xparticles1, xparticles2)
#   epsilon <- parameters$epsilon
#   desired_alpha <- parameters$desired_alpha
#   # wd <- wasserstein_auto(p = normweights1, q = normweights2, M, epsilon * median(M), desired_alpha)
#   wd <- wasserstein_auto(p = normweights1, q = normweights2, M, epsilon * median(M), desired_alpha)
#   couplingmatrix <- wd$transportmatrix
#   
#   u_1 <- rowSums(couplingmatrix)
#   u_2 <- colSums(couplingmatrix)
#   #  
#   nu <- pmin(normweights1 / u_1, normweights2 / u_2)
#   alpha <- min(nu)
#   # R1 <- (normweights1 - alpha * u_1) / (1 - alpha)
#   R2 <- (normweights2 - alpha * u_2) / (1 - alpha)
#   coupled <- (runif(nparticles) < alpha)
#   ncoupled <- sum(coupled)
#   #
#   ancestors2 <- rep(0, nparticles)
#   if (ncoupled > 0){
#     for (iparticle in which(coupled)){
#       coupl_row <- couplingmatrix[ancestors_ref[iparticle],]
#       coupl_row <- coupl_row / sum(coupl_row)
#       ancestors2[iparticle] <- systematic_resampling_n(coupl_row, 1, runif(1))
#     }
#     if (ncoupled < nparticles){
#       u <- runif(1)
#       ancestors2[!coupled] <- systematic_resampling_n(R2, nparticles - ncoupled, u)
#     }
#   } else {
#     u <- runif(1)
#     ancestors2 <- systematic_resampling_n(R2, nparticles, u)
#   }
#   return(ancestors2)
# }
# 
# 
# 
