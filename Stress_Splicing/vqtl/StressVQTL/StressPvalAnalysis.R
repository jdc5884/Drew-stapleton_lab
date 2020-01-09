### Analyzing vqtl Stress P-vals ###
library("tidyverse")

setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/vqtl")


datH = read.csv(file = "AdditiveModelHybridStress_Output.csv")
datHvqtl = datH[c(2:5,8:9,16:19)] #dataset containing only vqtl variables
small5H = head(sort(datHvqtl$vQTL.asymp.p), 5)
boxplot(datHvqtl$vQTL.asymp.p, main = "Boxplot of vQTL asymptotic P-values",
        col = "light blue", ylab = "Asymptotic P-value", xlab = "Hybrid Stress")
#boxplot of the 5 smallest p-vals
boxplot(small5H, main = "Boxplot of 5 smallest asymptotic P-values", 
        col = "steelblue4", ylab = "Asymptotic P-value", xlab = "Hybrid Stress")

sig_pH = datHvqtl %>% filter(datHvqtl$vQTL.asymp.p<0.0001) 

datI = read.csv(file = "AdditiveModelInbredStress_Output.csv")
datIvqtl = datI[c(2:5,8:9,16:19)] #dataset containing only vqtl variables
small5I = head(sort(datIvqtl$vQTL.asymp.p), 5)
boxplot(datIvqtl$vQTL.asymp.p, main = "Boxplot of vQTL asymptotic P-values",
        col = "light pink", ylab = "Asymptotic P-value", xlab = "Inbred Stress")
#boxplot of the 5 smallest p-vals
boxplot(small5I, main = "Boxplot of 5 smallest asymptotic P-values", 
        col = "pink3", ylab = "Asymptotic P-value", xlab = "Inbres Stress")

#boxplot of all p-values
boxplot(datHvqtl$vQTL.asymp.p, datIvqtl$vQTL.asymp.p, main = "Boxplot of vQTL asymptotic P-values",
        col = c("light blue", "light pink"), ylab = "Asymptotic P-value",
        names = c("Hybrid", "Inbred"))
#boxplot of both of the smallest 5
boxplot(small5H, small5I, main = "Boxplot of vQTL asymptotic P-values",
        col = c("steelblue4", "pink3"), ylab = "Asymptotic P-value",
        names = c("Hybrid", "Inbred"))
