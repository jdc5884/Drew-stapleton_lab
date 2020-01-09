#install.packages("vqtl")
#install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
library("tidyverse")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/vqtl/StressVQTL")

stressprod = read.cross(file = "../../Heirarchical/FullInbvqtlinput.csv")
stressprod <- drop.nullmarkers(stressprod)

##### CORTY code #####
stressprod$pheno$month = factor(stressprod$pheno$month)
# ggplot of stress type and height
# library(tidyverse)
# stressprod$pheno %>%
#   ggplot(mapping = aes(x = BreedType, y = stress)) +
#   geom_jitter(width = 0.2)

stressprod <- calc.genoprob(stressprod)
saveRDS(object = stressprod,
        file = 'data/inbredplotprep.RDS')


outv <- scanonevar(cross = stressprod,
                   mean.formula = stress ~ mean.QTL.add,
                   var.formula = ~ var.QTL.add ,
                   return.covar.effects = TRUE)


#outv$result %>% glimpse()

write.csv(outv$result, file = "AdditiveModelInbredStress_Output.csv")

### don't need .dom ###