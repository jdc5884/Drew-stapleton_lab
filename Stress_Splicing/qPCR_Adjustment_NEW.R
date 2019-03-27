### QPCR PLATE FRAMING & ADJUSTMENT MODEL ###

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")

# PC Directory
#setwd()

### READ IN DERIVATIVE DATA###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)

### INITIAL DATA FRAMING ###
# Remove extra labels row, 
deriv = deriv[-1,]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Remove extra labels row
deriv=deriv[-1,]
# Rename columns
colnames(deriv)=c("reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# Remove extra labels row
deriv=deriv[-1,]
write.csv(deriv, file="2018_11_qPCR_output_withHeaders.csv")
# Create new transposed data set
deriv2=read.csv(file="2018_11_qPCR_output_withHeaders.csv", header=TRUE) 
# Remove unneeded labeling column
deriv2=deriv2[,-1]
# Indicate if sample is NTC (negative control)
deriv2['sampleID_NTC'] = grepl('NTC', deriv2$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv2$sampleID_NTC)
deriv2 = deriv2[-ntc,]
deriv2 = deriv2[,-c(5,6)]
# Indicate if sample is 'Plus' or 'Minus'
deriv2['sampleID_Minus'] = grepl('minus', deriv2$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv2$sampleID_Minus)
deriv2 = deriv2[-minus,]
deriv2 = deriv2[,-c(5)]
# Write CSV file 
write.csv(deriv2, file="2018_11_qPCR_output_New.csv")
# Re-read in CSV values --> this allows for-loop to work properly
deriv2=read.csv(file="2018_11_qPCR_output_New.csv", header=TRUE) 
# Remove labels column
deriv2=deriv2[,-1]
### COMPLETED INITIAL DATA FRAMING ###

### CALIBRATED DATA FRAME ###
# Create/Write data frame for Calibrated values
calib_df = deriv2 %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_df = calib_df[order(calib_df$starting_quantity),]



# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_data$starting_quantity)){
  sq <- calib_data$starting_quantity[i]
  if(i %% 6 == 1){
    startq = c(startq,sq,sq,sq)
  }
  val <- toString(calib_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 = c(test1, calib_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP = c(allP, calib_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib_data = cbind(startq, test1, allP)
# Format starting quantity values as decimals, not scientific notation
calib_data[,1]=format(calib_data[,1], scientific=FALSE)
write.csv(calib_data, file="2018_11_Calibrated_Data_Frame.csv")
### COMPLETED CALIBRATED DATA FRAME ###

### ADJUSTMENT MODEL - CALIBRATED DATA ### 
data=read.csv(file="2018_11_Calibrated_Data_Frame.csv", header=TRUE) 
# Create empty vectors for for-loop input
adj_val = c()
allP = c()
startq = c()
ratio = data$allP/data$test1
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allP[3*i - 2], data$allP[3*i - 1], data$allP[3*i])
  t_y <- c(data$test1[3*i - 2], data$test1[3*i - 1], data$test1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}
adjusted_test1 <- test1 + adj_val
# Creating the adjustment model lm(y-axis~x-axis)
# Changed adj_val^2 to adj_val to try to make the model better --> VERIFY THIS WITH DR. WANG
adj_model <- lm(adj_val~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
summary(adj_model)  
par(mfrow = c(2,2))
plot(adj_model)
### COMPLETED ADJUSTMENT MODEL - CALIBRATED DATA ###

### EXPERIMENTAL DATA FRAME ###
# Create/Write data frame for Calibrated values
exp_df = deriv2 %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_df = exp_df[order(exp_df$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
exp_df = exp_df[-1,]
exp_df = exp_df[-nrow(exp_df),]
write.csv(exp_df, file="2018_11_Experimental_Data_Frame.csv")
# Create new transposed data set
exp_data = read.csv(file="2018_11_Experimental_Data_Frame.csv", header=TRUE) 
# Remove extra labeling column
exp_data = exp_data[,-1]
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]


# Remove 
# Include sampleID in all output to match with vQTL analysis

# expdata2 will be used to compare to the for-loop output, exp_data, to verify processed correctly
expdata2=exp_data

#REMOVE SAMPLES WITH MISSING ALLP/TEST1 VALUES; 
#CONFIRM WITH DR.S IF SAMPLES WITH UNUSUAL CP VALUES SHOULD BE OMITTED
exp_data=exp_data ##Finish this code to remove certain sampleID's -- PRIOR to for-loop

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID
exp_data = cbind(sampleID.exp, test1.exp, allP.exp)
# Write CSV file
write.csv(exp_data, file="2018_11_Experimental_Data_Frame.csv")
### COMPLETED EXPERIMENTAL DATA FRAME ###

### ADJUSTMENT MODEL - EXPERIMENTAL DATA ### 
# Using the adjustment model on the expiremental data
new = data.frame(ratio = exp_data$all_productsPrimers_Cp1/exp_data$test1_Cp1)
predict(adj_model, new , interval = "confidence")
#---> Fill in any remaining parts of experimental adjusmtent model

### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###


#___________________________________#
## REMAINING QUESTIONS ABOUT MODEL:

# DO WE INCLUDE THE 50-50 SHIFT IN OUR MODEL/PROJECT? --> CONFIRM WITH DR.S, WITH CALIBRATED OR EXPERIMENTAL?

# 50/50 Adjustment
test1_50_50 <- data$test1Cp150_50
all_p_50_50 <- data$allproductsCp1_50_50

# The ratio of the cp values
ratio_50 <- all_p_50_50/test1_50_50
# Appending to the Calibration data
cbind(data, ratio_50)
# Create empty vectors for for-loop to input cpD1 values
adj_val_50_50 = vector()
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allproductsCp1_50_50[3*i - 2], data$allproductsCp1_50_50[3*i - 1], data$allproductsCp1_50_50[3*i])
  t_y <- c(data$test1Cp150_50[3*i - 2], data$test1Cp150_50[3*i - 1], data$test1Cp150_50[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val_50_50 <- c(adj_val_50_50, adj, adj, adj)
}

# Creating the adjustment model lm(y-axis~x-axis)
adj_model_50 <- lm(adj_val_50_50^2~ratio_50) #50:50 shift - another type of calibration method   
#Run adjustment model on experimental data to get 
summary(adj_model_50)
plot(adj_model_50)

##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##
##__________________________________##

#####______OLD/UNUSED CODE______######

#INITIAL FRAMING#
#Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
#Remove extra labels row, 
deriv = deriv[-1,]
#Reorder columns to match plate data frame
deriv = deriv[c(2,3,4,5)]
#Use first row data as column names
names(deriv) = lapply(deriv[1, ], as.character)
deriv <- deriv[-1,] 
#Sort data by sample ID
deriv = deriv[order(deriv$unique_sampleID),]
#Drop extra labels row, automatically ordered to the last row
deriv = deriv[-nrow(deriv),]
#Give column names
colnames(deriv) = c("reaction_type", "sampleID", "starting_quantity", "cpD1")
#Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
#Remove NTC samples and indicator (T/F) column
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(5,6)]
#Remove 'Minus' values (include only gblock+ values)
deri = which(!calib_df$sampleID_Plus)
calib_df = calib_df[-minus,]

#Create data frame for Calibrated values
calib_df = dat2 %>% filter(str_detect(sampleID, "g"))
write.csv(calib_df, file="2018_6_1_Calibrated_Data_Frame.csv")

#Create data frame for Experimental values
exp_df = dat2 %>% filter(str_detect(sampleID, "g")==FALSE)
write.csv(exp_df, file="2018_6_1_Experimental_Data_Frame.csv")

#______________NEW FRAMING --> Replace with working previous____________________#
#Import the separate CSV file containing calculated cycle derivative values (cpD1, cpD2) for (MONTH)

#In the case of having two separate CSV files of calculated derivatives,
#use this code to combine, prior to the following transpositions:

deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
deriv=cbind(deriv.1, deriv.2)

#In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)

#INITIAL FRAMING#
#Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
#Remove extra labels row, 
deriv = deriv[-1,]
#Reorder columns to match plate data frame
deriv = deriv[c(2,3,4,5)]
#Use first row data as column names
names(deriv) = lapply(deriv[1, ], as.character)
deriv <- deriv[-1,] 
#Sort data by sample ID
deriv = deriv[order(deriv$unique_sampleID),]
#Drop extra labels row, automatically ordered to the last row
deriv = deriv[-nrow(deriv),]
#Give column names
colnames(deriv) = c("reaction_type", "sampleID", "starting_quantity", "cpD1")
#Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
#Remove NTC samples and indicator (T/F) column
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(5)]

# Read in Expiremental Data
calib_df <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv", head=TRUE)     
# Format starting quantity as numeric
calib_df <- calib_df[c(4,2,3,5,6)]
# Remove first extra labeling row --->> ADD THIS TO FRAMING CODE SO WONT MANUALLY
#exp_data <- exp_data[-1,]
# Sort data by sample ID
calib_df <- calib_df[order(calib_df$sampleID),]
write.csv(calib_df, file="2018_6_1_Calibrated_Data_Frame.csv")

#Create data frame for Calibrated values -- by detecting if a "g" is in sample ID, i.e. there is a known starting quantity
calib_df = deriv %>% filter(str_detect(sampleID, "g"))
# Indicate if sample is 'Plus' or 'Minus'
calib_df['sampleID_Plus'] = grepl('plus', calib_df$sampleID)
# Remove 'Minus' values (include only gblock+ values)
minus = which(!calib_df$sampleID_Plus)
calib_df = calib_df[-minus,]
# Sort data by starting quantity, remove 'Plus' indicator column
calib_df = calib_df[with(calib_df,order(calib_df$starting_quantity, calib_df$reaction_type)),] 
calib_df = calib_df[,-c(5)]
# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_df$starting_quantity)){
  sq <- calib_df$starting_quantity[i]
  if(i %% 6 == 1){
    startq = c(startq,sq,sq,sq)
  }
  val <- toString(calib_df$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 = c(test1, calib_df$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP = c(allP, calib_df$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib_data = cbind(startq, test1, allP)

#EXPERIMENTAL DATA FRAME
#Create data frame for Experimental values
exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_df$sampleID)){
  id.exp = toString(exp_df$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_df$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_df$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_df$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID
exp_data = cbind(sampleID.exp, test1.exp, allP.exp)

##ADJUSTMENT MODEL### #data = calibrated d.f. ==> calib_data#
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$all_productsCp1[3*i - 2], data$all_productsCp1[3*i - 1], data$all_productsCp1[3*i])
  t_y <- c(data$test1_Cp1[3*i - 2], data$test1_Cp1[3*i - 1], data$test1_Cp1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}

adjusted_test1 <- test1 + adj_val

# Creating the adjustment model lm(y-axis~x-axis)
adj_model <- lm(adj_val^2~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
summary(adj_model)  

par(mfrow = c(2,2))
plot(adj_model)

# Using the adjustment model on the expiremental data
new = data.frame(ratio = exp_data$all_productsPrimers_Cp1/exp_data$test1_Cp1)
predict(adj_model, new , interval = "confidence")

# 50/50

test1_50_50 <- data$test1Cp150_50
all_p_50_50 <- data$allproductsCp1_50_50

# The ratio of the cp values
ratio_50 <- all_p_50_50/test1_50_50
# Appending to the Calibration data
cbind(data, ratio_50)

adj_val_50_50 = vector()

# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allproductsCp1_50_50[3*i - 2], data$allproductsCp1_50_50[3*i - 1], data$allproductsCp1_50_50[3*i])
  t_y <- c(data$test1Cp150_50[3*i - 2], data$test1Cp150_50[3*i - 1], data$test1Cp150_50[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val_50_50 <- c(adj_val_50_50, adj, adj, adj)
}

# Creating the adjustment model lm(y-axis~x-axis)
adj_model_50 <- lm(adj_val_50_50^2~ratio_50) #50:50 shift - another type of calibration method   
                                              #Run adjustment model on experimental data to get 
summary(adj_model_50)
plot(adj_model_50)

# Convert starting quantity values into data frame to bind columns
#startq = as.data.frame(startq)

# Transpose test1.exp
#test1.exp <- t(test1.exp)
# Create allP.exp as a data frame; Transpose to single column
#allP.exp <- as.data.frame(allP.exp)
#allP.exp <- t(allP.exp)
# Make startq.exp as data frame
#sampleID.exp <- as.data.frame(sampleID.exp)
#e <- cbind(startq.exp)

#exp_data <- as.double(exp_data$cpD1)
### ERROR:  After converting values from integers to double (req'd to complete exp.df), 
# error returned "$ operator is invalid for atomic vectors"  

###ERROR in exp_data:  Will not read values as double, not integers, even though values are being pulled from 
# equivalent-types CSV files 

# Create two data frames from reaction type (test1 or all_products) 
#test1.exp_data = exp_data %>% filter(str_detect(exp_data$reaction_type, "test1") == TRUE)
#allP.exp_data = exp_data %>% filter(str_detect(exp_data$reaction_type, "all_products") == TRUE)
# Merge two data frames by sample ID
#exp.df <- merge(test1.exp_data, allP.exp_data, by="sampleID")

  # Itterating through each set of (3) observations performing U-Stats on each set of inputs
  for (i in 1:(nrow(exp_data)/3)){
    t_x <- c(exp_data$all_productsCp1[3*i - 2], exp_data$all_productsCp1[3*i - 1], exp_data$all_productsCp1[3*i])
    t_y <- c(exp_data$test1_Cp1[3*i - 2], exp_data$test1_Cp1[3*i - 1], exp_data$test1_Cp1[3*i])
    adj <- mean(outer(t_x, t_y, "-"))
    adj_val <- c(adj_val, adj, adj, adj)
  }

# To complete following for-loop:
## cp Values of Calibration data
## test1 <- data$cpD1
## all_prod <- data$all_productsCp1
## #The ratio of the cp values
## ratio <- all_prod/test1
###  Appending to the Calibration data
## cbind(data, ratio)
# Adjustment value Vector
## adj_val = vector()

#Old calibration data frame code:

# Create empty vectors for for-loop to input cpD1 values
#test1 = c()
#allP = c()
#startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
#for(i in 1:length(calib_data$starting_quantity)){
  #sq <- calib_data$starting_quantity[i]
  #if(i %% 6 == 1){
    startq <- c(startq,sq,sq,sq)
  }
  val <- toString(calib_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 <- c(test1, calib_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP <- c(allP, calib_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib.df = cbind(startq, test1, allP)

#Old experimental data frame code:

# Read in Expiremental Data
exp_data <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv", head=TRUE)     
# Format starting quantity as numeric
exp_data <- exp_data[c(4,2,3,5,6)]
# Remove first extra labeling row --->> ADD THIS TO FRAMING CODE SO WONT MANUALLY
#exp_data <- exp_data[-1,]
# Sort data by sample ID
exp_data <- exp_data[order(exp_data$sampleID),]

#Use first row data as column names
#names(deriv) = lapply(deriv[1, ], as.character)
#deriv <- deriv[-1,] 

setwd("/Users/andrewnorris/stapleton_lab/Stress_Splicing/2018_6")
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

#Remove extra labels row
#datT=datT[-1,]
write.csv(datT, file="2018_6_1_qPCR_output_withHeaders_T.csv")
#YEAR_MONTH_PLATE_qPCR_output.csv

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
startq.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$starting_quantity)){
  sq <- exp_data$starting_quantity[i]
  if(i %% 6 == 1){
    startq = c(startq,sq,sq,sq)
  }
  val <- toString(calib_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 = c(test1, calib_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP = c(allP, calib_data$cpD1[i])
  }
}
#Bind test1 and allProd cpD1 values by starting quantity
calib_data = cbind(startq, test1, allP)
#Format starting quantity values as decimals, not scientific notation
calib_data[,1]=format(calib_data[,1], scientific=FALSE)
write.csv(calib_data, file="2018_11_Calibrated_Data_Frame.csv")
###COMPLETED CALIBRATED DATA FRAME###

###ADJUSTMENT MODEL### 
data=read.csv(file="2018_11_Calibrated_Data_Frame.csv", header=TRUE) 

adj_val = c()
allP = c()
startq = c()

# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$all_productsCp1[3*i - 2], data$all_productsCp1[3*i - 1], data$all_productsCp1[3*i])
  t_y <- c(data$test1_Cp1[3*i - 2], data$test1_Cp1[3*i - 1], data$test1_Cp1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}

adjusted_test1 <- test1 + adj_val

# Creating the adjustment model lm(y-axis~x-axis)
adj_model <- lm(adj_val^2~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
summary(adj_model)  

par(mfrow = c(2,2))
plot(adj_model)

# Using the adjustment model on the expiremental data
new = data.frame(ratio = exp_data$all_productsPrimers_Cp1/exp_data$test1_Cp1)
predict(adj_model, new , interval = "confidence")

#Create/Write data frame for Experimental values
exp_df = deriv.2 %>% filter(str_detect(sampleID, "g")==FALSE)
write.csv(exp_df, file="2018_6_1_Experimental_Data_Frame.csv")

#Create new transposed data set
dat2=read.csv(file="2018_6_1_qPCR_output_withHeaders_T.csv", header=TRUE) 
#Format starting quantity values 
format(dat2$starting_quantity, scientific=FALSE)
#Remove labels column
dat2=dat2[,-1]
