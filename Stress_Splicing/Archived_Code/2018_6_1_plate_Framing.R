#Separating qPCR 6-1-2018 plate data into Calibrated and Experimental data frames

library(tidyverse)
setwd("/Users/andrewnorris/stapleton_lab/Stress_Splicing")
#Import qPCR plate raw data
dat=read.csv(file = "2018_6_1_plate.csv", header=FALSE)
#Remove plate ID from raw data and blank row
dat=dat[-c(1,5),]
#Transpose raw data so headers at top
datT=as.data.frame(t(dat), header=TRUE)
#Give column names
colnames(datT)=c("reaction_type", "sampleID", "starting_quantity", "Cycle_number")
#Remove extra labels row
datT=datT[-1,]
write.csv(datT, file="2018_6_1_plate_T.csv")
#Create new transposed data set
dat2=read.csv(file="2018_6_1_plate_T.csv", header=TRUE)
#Remove labels column
dat2=dat2[,-1]

#Create data frame for Calibrated values
calib_df = dat2 %>% filter(str_detect(sampleID, "g"))
write.csv(calib_df, file="Calibrated_Data_Frame.csv")

#Create data frame for Experimental values
exp_df = dat2 %>% filter(str_detect(sampleID, "g")==FALSE)
write.csv(exp_df, file="Experimental_Data_Frame.csv")



