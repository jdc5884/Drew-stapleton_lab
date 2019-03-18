#Separating qPCR output 2018_(MONTH)_(#) plate data into Calibrated and Experimental data frames

library(tidyverse)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")

# PC Directory
#setwd()

#In the case of having two data sets per individual month separated by test1 or all_products,
#use this code to combine, prior to the following transpositions:

dat.allP<-read.csv(file = "2018_11_1_plate",head=TRUE)
dat.test1<-read.csv(file = "2018_11_2_plate",head=TRUE)
dat=merge(dat.allP, dat.test1, head=TRUE)

#In the case of having one data set per individual month containing both test1 and all_products, use this code: 
#dat=read.csv(file = "2018_6_1_qPCR_output_withHeaders.csv", header=FALSE)


#### INCLUDE CODE TO TAKE MB's DERIVATIVE VALUES AND INPUT/REPLACE RAW CYCLE DATA ---> work w/ MB on this ####


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

