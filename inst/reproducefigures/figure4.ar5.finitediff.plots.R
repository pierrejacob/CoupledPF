# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPF)
library(ggthemes)
library(doRNG)
kalman_module <- Module( "kalman_mod", PACKAGE = "CoupledPF")
ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)

load("ar_correlation_resamplingcomparison.R1000D5T1000N128.RData")
dimension <- 5
experiments$dimension <- dimension

kf.finitediff <- function(theta, h){
  kalman_ll <- function(alpha){
    LGModel <- new(kalman_module$LinearGaussian)
    LGModel$set_multivariate_parameters(alpha, dimension)
    LGModel$set_observations(observations)
    Kalman <- new(kalman_module$Kalman)
    Kalman$setLinearGaussian(LGModel)
    Kalman$filtering()
    return(Kalman$getLL())
  }
  ll_plus <- kalman_ll(theta + h)
  ll_minus <- kalman_ll(theta - h)
  kfgrad <- (ll_plus - ll_minus) / (2*h)
  return(kfgrad)
}
hs <- unique(experiments$h)
kf.df <- foreach (h = hs, .combine = rbind) %do% {
  data.frame(dimension = dimension, h = h, kfgrad = kf.finitediff(0.3, h))
}

# allkf.df <- rbind(allkf.df, kf.df)
# allkf.df
# summary(all.df %>% filter(dimension == 1) %>% select(h))

experiments$method <- factor(experiments$method, levels = c("systematic", "Hilbert", "index-matching", "transport"),
                             labels = c("CRN+systematic  ", "CRN+sorted  ", "CRN+index-coupled  ", "CRN+transport  "))
# 
kf.df
experiments %>% head
experiments %>% group_by(dimension, h, method) %>% summarise(grad = mean((ll2 - ll1)/(2*h)))

experiments %>% group_by(dimension, h, method) %>% summarise(m1 = mean(ll1), m2 = mean(ll2))
experiments %>% group_by(dimension, h, method) %>% summarise(v1 = var(ll1), v2 = var(ll2)) %>% ungroup() %>%
  mutate(varfd = (v1 + v2) / (4*h^2))

# colors for independent, systematic, sorted, index-matching, transport
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#0072B2", "#F0E442")
cbbPalette <- cbbPalette[2:5]

# true gradient
correlation.df <- experiments %>% group_by(h, dimension, method) %>% summarise(correlation = cor(ll1, ll2))
head(correlation.df)
correlation.df$h <- factor(correlation.df$h)
correlation.df <- correlation.df %>% mutate(gain = 1/(1 - correlation))
g <- ggplot(correlation.df, 
            aes(x = h, y = gain, group = method, fill = method))
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g + scale_fill_manual(name = "", values = cbbPalette, guide = guide_legend(title.position = "top", nrow = 1)) + ylab("variance reduction factor")
g <- g + geom_text(aes(label = round(gain, 0)), position = position_dodge(width = 0.9), vjust = -0.1, size = 6)
g <- g + ylim(0, 530)
g

filename <- paste0("finitediff.ar5.R", nrep, "D", dimension, "T", datalength, "N", nparticles, ".pdf")
ggsave(plot = g, filename = filename, width = 12, height = 6)

