#install.packages("vqtl")
#install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
library("tidyverse")
library("ggplot2")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/vqtl/StressVQTL")

### Inbred Plot Work ###

stressprodInbred = readRDS(file = "data/inbredplotprep.RDS")

# #cross plots of p-values for stress
# xlims <- c(-0.3, 0.5)
# ylims <- c(0.6, 1.6)
# plot(p1 <- mean_var_plot_model_based(cross = stressprod,
#                                      phenotype.name = stressprod$pheno$stress,
#                                      focal.groups = NULL,
#                                      point_size = 3,
#                                      xlim = xlims,
#                                      ylim = ylims))
# 
# ggplot2::ggsave(plot = p1, filename = 'plots/mean_var_plot_inbred.pdf', height = 3, width = 4)

plot(stressprodInbred, main = "Inbred vQTL Plots")

so_p1 <- scanone(cross = stressprodInbred,
                 pheno.col = 'stress')
sov_p1 <- scanonevar(cross = stressprodInbred,
                     mean.formula = stress ~ mean.QTL.add + mean.QTL.dom,
                     var.formula = ~ var.QTL.add + var.QTL.dom)

ymax <- 6

plot(p1_LOD_scan <- plot(x = sov_p1, y = so_p1, ymax = ymax))


### Hybrid Plot Work ###

stressprodHybrid = readRDS(file = "data/hybridplotprep.RDS")

plot(stressprodHybrid)




### running through all of Corty's code to get graphs ###


so_p1 <- scanone(cross = stressprodInbred, pheno.col = 'stress')
sov_p1 <- scanonevar(cross = stressprodInbred,
                     mean.formula = stress ~ mean.QTL.add,
                     var.formula = ~ var.QTL.add)
