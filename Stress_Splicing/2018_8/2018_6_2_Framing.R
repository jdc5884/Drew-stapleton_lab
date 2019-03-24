#Separating qPCR output 2018_6_2 plate data into Calibrated and Experimental data frames

library(tidyverse)
setwd("/Users/andrewnorris/stapleton_lab/Stress_Splicing/2018_8")
#Import qPCR plate raw data -- CSV file including derivatives is in MB's repository, copied to AN's
dat=read.csv(file = "2018_6_1_qPCR_output_withHeaders.csv", header=FALSE)
#Remove plate ID from raw data and blank row
dat=dat[-c(1,5),]
#Transpose raw data so headers at top
datT=as.data.frame(t(dat), header=TRUE)
#Give column names
colnames(datT)=c("reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
#Remove extra labels row
datT=datT[-1,]
write.csv(datT, file="2018_6_1_qPCR_output_withHeaders_T.csv")
#Create new transposed data set
dat2=read.csv(file="2018_6_1_qPCR_output_withHeaders_T.csv", header=TRUE) 
#Format starting quantity values 
format(dat2$starting_quantity, scientific=FALSE)
#Remove labels column
dat2=dat2[,-1]

#Create data frame for Calibrated values
calib_df = dat2 %>% filter(str_detect(sampleID, "g"))
write.csv(calib_df, file="2018_6_1_Calibrated_Data_Frame.csv")

#Create data frame for Experimental values
exp_df = dat2 %>% filter(str_detect(sampleID, "g")==FALSE)
write.csv(exp_df, file="2018_6_1_Experimental_Data_Frame.csv")

