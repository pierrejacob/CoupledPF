#'@export
CR_hilbert_given <- function(xparticles1, xparticles2, normweights1, normweights2, parameters = list(),
                             ancestors_ref, uniforms){
  nparticles <- ncol(xparticles1)
  horder2 <- hilbert_order(xparticles2)
  nw_sorted2 <- normweights2[horder2]
  ancestors2 <- systematic_resampling_n(nw_sorted2, nparticles, uniforms[1])
  return(horder2[ancestors2])
}