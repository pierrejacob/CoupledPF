#'@rdname CR_hilbert
#'@title Coupled Resampling: Hilbert
#'@description This function performs coupled resampling based on sorting the samples
#'according to their projection on the Hilbert space filling curve.
#'@return Two vectors of ancestors, column-binded in a matrix.
#'@export
CR_hilbert <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  nparticles <- ncol(xparticles1)
  horder1 <- hilbert_order(xparticles1)
  horder2 <- hilbert_order(xparticles2)
  nw_sorted1 <- normweights1[horder1]
  nw_sorted2 <- normweights2[horder2]
  u <- runif(1)
  ancestors1 <- systematic_resampling_n(nw_sorted1, nparticles, u)
  ancestors2 <- systematic_resampling_n(nw_sorted2, nparticles, u)
  return(cbind(horder1[ancestors1], horder2[ancestors2]))
}