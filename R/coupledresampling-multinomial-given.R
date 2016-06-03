#'@export
CR_multinomial_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                parameters = list(),
                                ancestors_ref, uniforms){
  nparticles <- ncol(xparticles1)
  ancestors2 <- multinomial_resampling_n(normweights2, nparticles)
  return(ancestors2)
}
