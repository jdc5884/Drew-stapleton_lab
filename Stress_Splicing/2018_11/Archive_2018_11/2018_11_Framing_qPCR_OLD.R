#Separating qPCR output 2018_11 plate data into Calibrated and Experimental data frames

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")

# PC Directory
#setwd()

#Import the separate CSV file containing calculated cycle derivative values (cpD1, cpD2) for (MONTH)

#In the case of having two separate CSV files of calculated derivatives,
#use this code to combine, prior to the following transpositions:

deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv",header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv",header=FALSE)
deriv=cbind(deriv.1, deriv.2)

#In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)

#Transpose derivatives to be in equivalent format as raw plate data
deriv=as.data.frame(t(deriv), header=TRUE)
#Remove extra labels row, 
deriv=deriv[-1,]
#Reorder columns to match plate data frame
deriv=deriv[c(2,3,4,5)]
#Use first row data as column names
names(deriv) <- lapply(deriv[1, ], as.character)
deriv <- deriv[-1,] 
#Sort data by sample ID
deriv=deriv[order(deriv$unique_sampleID),]
#Drop extra labels row, automatically ordered to the last row
deriv=deriv[-nrow(deriv),]
#Give column names
colnames(deriv)=c("reaction_type", "sampleID", "starting_quantity", "cpD1")
#Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] <- grepl('NTC', deriv$sampleID)
#Remove NTC samples and indicator (T/F) column
ntc <- which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv=deriv[,-c(5)]

#CALIBRATED DATA FRAME
#Create data frame for Calibrated values -- by detecting if a "g" is in sample ID, i.e. there is a known starting quantity
calib_df = deriv %>% filter(str_detect(sampleID, "g"))
# Indicate if sample is 'Plus' or 'Minus'
calib_df['sampleID_Plus'] <- grepl('plus', calib_df$sampleID)
# Remove 'Minus' values (include only gblock+ values)
minus <- which(!calib_df$sampleID_Plus)
calib_df = calib_df[-minus,]
# Sort data by starting quantity, remove 'Plus' indicator column
calib_df = calib_df[with(calib_df,order(calib_df$starting_quantity, calib_df$reaction_type)),] 
calib_df=calib_df[,-c(5)]

###Change starting_quantity to have consistent number of digits inc. trailing zeros
#format(calib_df$starting_quantity, trim=FALSE, nsmall=5)
#signif(calib_df$starting_quantity, digits=7)
#formatC(calib_df$starting_quantity, big.interval=2, small.interval=3)
#format.data.frame(calib_df$starting_quantity, )
#format(calib_df$starting_quantity, big.mark=2, small.mark=3)
#str_pad(calib_df$starting_quantity, 7, pad="0")
#sprintf

# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_df$starting_quantity)){
  sq <- calib_df$starting_quantity[i]
  if(i %% 6 == 1){
    startq <- c(startq,sq,sq,sq)
  }
  val <- toString(calib_df$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 <- c(test1, calib_df$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP <- c(allP, calib_df$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib_data = cbind(startq, test1, allP)

#EXPERIMENTAL DATA FRAME
#Create data frame for Experimental values
exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)

#Prev code...
# Read in Expiremental Data
#exp_data <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv", head=TRUE)    

# Format starting quantity as numeric
#exp_data <- exp_data[c(4,2,3,5,6)]
# Remove first extra labeling row --->> ADD THIS TO FRAMING CODE SO WONT MANUALLY
#exp_data <- exp_data[-1,]
# Sort data by sample ID
#exp_data <- exp_data[order(exp_data$sampleID),]
# Remove NTC (negative control) values, not to be included in adjustments --> already completed in prev framing
#exp_data <- exp_data[-(187:196),]
# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_df$sampleID)){
  id.exp <- toString(exp_df$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp <- c(sampleID.exp, id.exp)
  }
  val <- toString(exp_df$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp <- c(test1.exp, exp_df$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp <- c(allP.exp, exp_df$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID
exp_data <- cbind(sampleID.exp, test1.exp, allP.exp)




#write.csv(exp_df, file="2018_6_1_Experimental_Data_Frame.csv")


###________________________
#In the case of having two data sets per individual month separated by test1 or all_products,
#use this code to combine, prior to the following transpositions:

dat.allP<-read.csv(file = "2018_11_1_plate.csv",header=FALSE)
dat.test1<-read.csv(file = "2018_11_2_plate.csv",header=FALSE)
dat=cbind(dat.allP, dat.test1)

#In the case of having one data set per individual month containing both test1 and all_products, use this code: 
#dat=read.csv(file = "2018_6_1_qPCR_output_withHeaders.csv", header=FALSE)


#Remove plate ID from raw data and blank row
dat=dat[-c(1,5),]
#Transpose raw data so headers at top
dat=as.data.frame(t(dat), header=TRUE)
#Use first row data as column names
names(dat) <- lapply(dat[1, ], as.character)
dat <- dat[-1,] 
#Sort data by sample ID
dat=dat[order(dat$unique_sampleID_number),]
#Remove raw cycle data values to replace with corresponding derivative values
dat=dat[,-c(4:44)]

#### INCLUDE CODE TO TAKE MB's DERIVATIVE VALUES AND INPUT/REPLACE RAW CYCLE DATA ---> work w/ MB on this #### 

#Extract columns containing derivative values
CPvals=data.frame(deriv$cpD1)



#Create new transposed data set
#dat2=read.csv(file="2018_6_1_qPCR_output_withHeaders_T.csv", header=TRUE) 
#Format starting quantity values 
format(deriv$starting_quantity, scientific=FALSE)

