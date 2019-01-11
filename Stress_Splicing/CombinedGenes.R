### Combining Genomes and Mo### 
library(stringr)

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
dat = read.csv(file = "Plant_Height.csv", header = TRUE)

#Take out unneeded IBMB###, NA, B73 loci
dat = dat[-(907:938),-(4:5)]

#Create Categorical Variables for PH207*Mo### and Mo### by gene breed
BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Outbred")
dat = cbind(dat, BreedType)

#Add in SNP info from Marker data CSV, beginning with column six
snp = read.csv(file = "IBM94markerset08seq.csv")
snp = snp[,-(1:5)]

#Create zero matrix to which data will input
relevant = data.frame(matrix(rep(0,length(dat$Genotype)*dim(snp)[1]), ncol = dim(snp)[1]))

#From "Genotype" values, match Plant Height data and Marker data by detecting last three digits of Mo###'s
dat2 = sapply(str_sub(dat$Genotype,-3,-1), function(x){
  column = which(str_sub(colnames(snp),-3,-1) == x)
  vect = data.frame(as.character(snp[,column]))
  return(vect)
})

#Add matched Mo### values to data frame 
dat2 = as.data.frame(matrix(unlist(dat2), nrow = dim(dat)[1], byrow = TRUE))

#Running into errors beginning with colnames(dat2)#
colnames(dat2) = colnames()
library(beepr)
beep()
dim(dat3);dim(snp)

#####Adding back in the Trait info#####
dat3 = cbind(dat2$Height,dat2[,1:3],dat3)
colnames(dat3) = c(colnames(dat2[1:3]),"Height",as.character(snp$markername))
dat3[1:10,1:10]
#####Adding marker location and chromosome#####
aux = matrix(snp$incre_new, nrow= 1)
aux = rbind(aux,snp$Chromosome)
other = as.data.frame(matrix(rep(0,8), nrow = 2))
aux = cbind(other,aux)
colnames(aux) = rep("",3239)
colnames(dat3) = rep("",3239)
dat4 = rbind(aux,dat3)
colnames(dat4) = c("Height", colnames(dat2[1:3]),as.character(snp$markername))
dat4[1:10,1:10]
write.csv(dat4, file = ,
          row.names = FALSE)
beep()
#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####