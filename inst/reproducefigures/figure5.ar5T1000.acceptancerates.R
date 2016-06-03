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


# nrep <- 10
# Ns <- c(64, 128, 256, 512)
nrep <- 5
Ns <- c(64,128,256)
datalength <- 1000
niterations <- 5000
# methods <- c("pmmh", "cpmmh.systematic",  "cpmmh.hilbert", "cpmmh.index")
methods <- c("pmmh", "cpmmh.systematic", "cpmmh.hilbert", "cpmmh.index", "cpmmh.transport")
# methods <- c("cpmmh.hilbert", "cpmmh.index", "cpmmh.transport")

accept.df <- data.frame()
for (irep in 1:nrep){
  for (nparticles in Ns){
    for (method in methods){
      a <- try(load(paste0(method, ".repeat", irep, "T", datalength, "N", nparticles, "M", niterations, ".RData")))
      if (inherits(a, "try-error")){
        pmmh_res <- list(acceptance_rate = NA)
      }
      accept.df_ <- data.frame(i = irep, method = method, N = nparticles, acceptrate = pmmh_res$acceptance_rate)
      accept.df <- rbind(accept.df, accept.df_)
    }
  }
}
# accept.df
# colors for independent, systematic, sorted, index-matching, transport
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#0072B2", "#F0E442")
# colors for independent, systematic, sorted, index-matching, transport
cbbPalette <- cbbPalette[c(1,2,3,4,5)]
# cbbPalette <- cbbPalette[c(3,4,5)]

accept.summary.df <- accept.df %>% group_by(method, N) %>% summarise(m = mean(100*acceptrate), s = sd(100*acceptrate))
accept.summary.df
accept.summary.df$N <- factor(accept.summary.df$N)
accept.summary.df$method <- factor(accept.summary.df$method, levels = methods,
                                   labels = c("independent  ", "CRN+systematic  ", "CRN+sorted  ", "CRN+index-coupled  ", "CRN+transport  "))
                                    # labels = c("CRN+sorted  ", "CRN+index-matching  ", "CRN+transport  "))
g <- ggplot(accept.summary.df, 
            aes(x = N, y = m, group = method, fill = method))
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g + scale_fill_manual(name = "", 
                           guide = guide_legend(title.position = "top", nrow = 1),
                           values = cbbPalette) + ylab("acceptance rate (%)")
g <- g + ylim(0, 8)
g <- g + geom_text(aes(label = round(m, 2)), position = position_dodge(width = 0.9), vjust = -0.1, size = 6)
g
filename <- paste0("pmmh.ar5.acceptancerates.pdf")
ggsave(plot = g, filename = filename, width = 13, height = 6)


