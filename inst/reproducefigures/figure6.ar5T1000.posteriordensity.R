# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(doRNG)
library(ggthemes)
library(mcmcse)
library(coda)

ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
nrep <- 5
Ns <- 256
datalength <- 1000
niterations <- 5000

methods <- c("pmmh", "cpmmh.hilbert", "cpmmh.index", "cpmmh.transport")

posterior.df <- data.frame()
for (irep in 1:nrep){
  for (nparticles in Ns){
    for (method in methods){
      load(paste0(method, ".repeat", irep, "T", datalength, "N", nparticles, "M", niterations, ".RData"))
      pmmh.chain.df <- melt(pmmh_res$chain)
      names(pmmh.chain.df) <- c("iteration", "component", "chain")
      pmmh.chain.df$rep <- irep
      pmmh.chain.df$N <- nparticles
      pmmh.chain.df$method <- method
      posterior.df <- rbind(posterior.df, pmmh.chain.df)
    }
  }
}

posterior.df %>% tail

posterior.df$method <- factor(posterior.df$method, levels = methods,
                                   labels = c("independent  ", "CRN+sorted  ", "CRN+index-coupled  ", "CRN+transport  "))

load("ar.kfmh.D5T1000.RData")
kfchain.df <- melt(kf_res$chain)
kfchain.df %>% tail
names(kfchain.df) <- c("iteration", "component", "chain")
posterior.df %>% tail
unique(posterior.df$method)
g <- ggplot(posterior.df %>% filter(method %in% c("independent  ", "CRN+index-coupled  ")), 
            aes(x = chain, group = interaction(method, rep))) +
  geom_density(adjust = 1) + facet_wrap(~ method)
g <- g + geom_density(data=kfchain.df, aes(x = chain, group = NULL), fill = "red", colour = "red", alpha = 0.25, 
                      adjust=1)
g <- g + ylim(0,300)
g <- g +  xlab(expression(theta))
g
filename <- paste0("pmmh.ar5.density.pdf")
ggsave(plot = g, filename = filename, width = 15, height = 5)


# ?geom_density
