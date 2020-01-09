#install.packages("vqtl")
#install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
library("tidyverse")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/vqtl/StressVQTL")

stressprod = readRDS(file = "data/inbredplotprep.RDS")
#cross plots of p-values for stress

xlims <- c(-0.3, 0.5)
ylims <- c(0.6, 1.6)
plot(p1 <- mean_var_plot_model_based(cross = stressprod,
                                     phenotype.name = "stress",
                                     focal.groups = NULL,
                                     point_size = 3,
                                     xlim = xlims,
                                     ylim = ylims))

ggplot2::ggsave(plot = p1, filename = 'plots/mean_var_plot_inbred.pdf', height = 3, width = 4)
