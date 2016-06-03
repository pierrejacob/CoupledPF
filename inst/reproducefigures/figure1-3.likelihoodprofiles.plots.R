# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(doRNG)
library(ggthemes)
ncores <- 8
registerDoMC(cores = ncores)
# load custom theme for the plots
setmytheme()
# fix the random seed
set.seed(17)


# dimension <- 5
dimension <- 1
load(paste0("ar", dimension, "data.RData"))
datalength <- 100
# datalength <- 1000

observations <- matrix(observations[1:datalength,], ncol = dimension)

theta <- 0.4

# nparticles <- 2^10 # For d = 5, T = 1000
nparticles <- 2^7 # For d = 1, T = 100
nrep <- 5
nrhos <- 20
rhos <- seq(from = 0.3, to = 0.4, length.out = nrhos)

filename <- paste0("likelihoodprofiles.D", dimension, "T", datalength, "N", nparticles, ".RData")
load(file = filename)
all.melt$method <- factor(all.melt$method, levels = c("independent", "CR+systematic", "CR+sorted", "CR+index-matching",  "CR+transport"),
labels = c("independent", "CRN+systematic", "CRN+sorted", "CRN+index-coupled", "CRN+transport"))

g <- ggplot(all.melt %>% filter(method %in% c("independent", "CRN+systematic")), aes(x = rhos, y = value, group = interaction(rep, method), colour = rep)) + geom_line() + geom_point()
g <- g + theme(legend.position = "none") + facet_wrap(~ method)
g <- g + scale_color_colorblind()
g <- g + xlab(expression(theta)) + ylab("log-likelihood")
g <- g + scale_x_continuous(breaks = c(0.3, 0.35, 0.4))
g <- g + geom_line(data = kfll, aes(x = rhos, y = ll, group = NULL, colour = NULL), colour = "red", size = 2)
g

filename <- paste0("likelihoodprofiles.independent.D", dimension, "T", datalength, "N", nparticles, ".pdf")
ggsave(plot = g, filename = filename, width = 10, height = 5)

g <- ggplot(all.melt %>% filter(method %in% c("CRN+transport", "CRN+index-coupled", "CRN+sorted")), aes(x = rhos, y = value, group = interaction(rep, method), colour = rep)) + geom_line() + geom_point()
g <- g + facet_wrap(~ method) + theme(legend.position = "none")
g <- g + scale_color_colorblind()
g <- g + xlab(expression(theta)) + ylab("log-likelihood")
g <- g + geom_line(data = kfll, aes(x = rhos, y = ll, group = NULL, colour = NULL), colour = "red", size = 2)
g <- g + scale_x_continuous(breaks = c(0.3, 0.35, 0.4))
g
filename <- paste0("likelihoodprofiles.crn.D", dimension, "T", datalength, "N", nparticles, ".pdf")
ggsave(plot = g, filename = filename, width = 15, height = 5)


dimension <- 5
datalength <- 1000
filename <- paste0("likelihoodprofiles.D", dimension, "T", datalength, "N", nparticles, ".RData")
load(file = filename)
all.melt$method <- factor(all.melt$method, levels = c("independent", "CR+systematic", "CR+sorted" , "CR+index-matching", "CR+transport"),
                          labels = c("independent", "CRN+systematic",  "CRN+sorted", "CRN+index-coupled", "CRN+transport"))


g <- ggplot(all.melt %>% filter(method %in% c("independent", "CRN+systematic")), aes(x = rhos, y = value, group = interaction(rep, method), colour = rep)) + geom_line() + geom_point()
g <- g + theme(legend.position = "none") + facet_wrap(~ method)
g <- g + scale_color_colorblind()
g <- g + xlab(expression(theta)) + ylab("log-likelihood")
g <- g + scale_x_continuous(breaks = c(0.3, 0.35, 0.4))
g <- g + geom_line(data = kfll, aes(x = rhos, y = ll, group = NULL, colour = NULL), colour = "red", size = 2)
# g

g <- ggplot(all.melt %>% filter(method %in% c("CRN+sorted", "CRN+index-coupled", "CRN+transport")), aes(x = rhos, y = value, group = interaction(rep, method), colour = rep)) + geom_line() + geom_point()
g <- g + facet_wrap(~ method) + theme(legend.position = "none")
g <- g + scale_color_colorblind()
g <- g + xlab(expression(theta)) + ylab("log-likelihood")
g <- g + geom_line(data = kfll, aes(x = rhos, y = ll, group = NULL, colour = NULL), colour = "red", size = 2)
g <- g + scale_x_continuous(breaks = c(0.3, 0.35, 0.4))
g
filename <- paste0("likelihoodprofiles.crn.D", dimension, "T", datalength, "N", nparticles, ".pdf")
ggsave(plot = g, filename = filename, width = 15, height = 5)

