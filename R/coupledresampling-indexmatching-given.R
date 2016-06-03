#' #' @export
#' CR_indexmatching_given <- function(xparticles1, xparticles2, normweights1, normweights2,
#'                                    parameters = list(),
#'                                    ancestors_ref, uniforms){
#'   nparticles <- length(normweights1)
#'   # common measure
#'   nu <- pmin(normweights1, normweights2)
#'   alpha <- sum(nu)
#'   # check if the weight vectors are equal, in which case we don't need to sweat too much
#'   if (alpha > 1-1e-20){
#'     ancestors2 <- ancestors_ref
#'     return(ancestors2)
#'   }
#'   mu <- nu / alpha
#'   # residuals
#'   R1 <- normweights1 -  nu
#'   R1 <- R1 / (1 - alpha)
#'   
#'   R2 <- normweights2 -  nu
#'   R2 <- R2 / (1 - alpha)
#'   
#'   P <- diag(nu, nparticles, nparticles) + (1-alpha)* matrix(R1, ncol = 1) %*% matrix(R2, nrow = 1) 
#'   # print(P)
#'   ancestors2 <- rep(0, nparticles)
#'   for (i in 1:nparticles){
#'     ancestors2[i] <- sample(1:nparticles, 1, replace = TRUE, prob = P[ancestors_ref[i],])
#'   }
#'   #
#'   return(ancestors2)
#' }

#'@rdname CR_indexmatching_given
#'@title Coupled Resampling: index-matching, given vector of ancestor
#'@description This function performs index-matching resampling on the second system
#'@return A vector of ancestors
#'@export
CR_indexmatching_given <- function(xparticles1, xparticles2, normweights1, normweights2,
                                   parameters = list(),
                                   ancestors_ref, uniforms){
  nparticles <- length(normweights1)
  # uniforms <- runif(nparticles + 1)
  return(indexmatching_given_cpp(nparticles, normweights1, normweights2, uniforms, ancestors_ref - 1) + 1)
}
#' #' @export
#' CR_indexmatching_given <- function(xparticles1, xparticles2, normweights1, normweights2,
#'                                              parameters = list(),
#'                                              ancestors_ref, uniforms){
#'   nparticles <- length(normweights1)
#'   # common measure
#'   nu <- pmin(normweights1, normweights2)
#'   alpha <- sum(nu)
#'   # check if the weight vectors are equal, in which case we don't need to sweat too much
#'   if (alpha > 1-1e-20){
#'     ancestors2 <- ancestors_ref
#'     return(ancestors2)
#'   }
#'   mu <- nu / alpha
#'   # residuals
#'   R2 <- normweights2 -  nu
#'   R2 <- R2 / (1 - alpha)
#'   #
#'   coupled <- (runif(nparticles) < alpha)
#'   ncoupled <- sum(coupled)
#'   ancestors2 <- rep(0, nparticles)
#'   if (ncoupled > 0){
#'     ancestors2[coupled] <- ancestors_ref[coupled]
#'     if (ncoupled < nparticles){
#'       # u <- runif(1)
#'       ancestors2[!coupled] <- sample(1:nparticles, nparticles - ncoupled, replace = TRUE, prob = R2)
#'       #systematic_resampling_n(R2, nparticles - ncoupled, u)
#'     }
#'   } else {
#'     # u <- runif(1)
#'     # ancestors2 <- systematic_resampling_n(R2, nparticles, u)
#'     ancestors2 <- sample(1:nparticles, nparticles, replace = TRUE, prob = R2)
#'   }
#'   return(ancestors2)
#' }
