% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_get_stochvol.R
\name{get_stochvol}
\alias{get_stochvol}
\title{Stochastic Volatility model}
\usage{
get_stochvol(dimension)
}
\value{
A list
}
\description{
This function returns a list with objects such as
* rinit  to sample from the initial distribution
* rtransition to sample from the transition
* dtransition to evaluate the transition density
* dmeasurement to evaluate the measurement density
* generate_randomness to evaluate the measurement density
* perturb_randomness to evaluate the measurement density
* dimension, which represents the dimension of the latent process
}

