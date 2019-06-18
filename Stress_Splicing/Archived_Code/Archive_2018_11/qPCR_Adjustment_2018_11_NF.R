########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")
# PC Directory
#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing/qPCR")

### READ IN DERIVATIVE DATA ###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)

########################################################## 
################### Initial Data Framing #################
########################################################## 

# Remove extra labels row and column 
deriv = deriv[-1,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Rename columns
colnames(deriv)=c("reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# Remove extra labels row
deriv=deriv[-1,]
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(5,6)]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
deriv = deriv[,-c(5)]

### COMPLETED INITIAL DATA FRAMING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 

# Create/Write data frame for Calibrated values
calib_df = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_df = calib_df[order(calib_df$starting_quantity),]
calib_df$starting_quantity = as.numeric(as.character(calib_df$starting_quantity))
calib_df$cpD1 = as.numeric(as.character(calib_df$cpD1))
calib_data = calib_df
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
calib_data = as.data.frame(cbind(startq, test1, allP))
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)
ratio = calib_data$allP/calib_data$test1
# Append ratios to data set
calib_data = cbind(calib_data, ratio)

### COMPLETED CALIBRATED DATA FRAME ###


########################################################## 
##### Ordinal Logicistic Regression Calibrated Data ######
##########################################################
calib_data$startq = ordered(calib_data$startq, levels = levels(calib_data$startq))
calib_data$ratio = allP/test1
# Ordinal Logistic
OLR = polr(startq~ratio,data = calib_data, Hess = TRUE)
summary(OLR)
(ctable <- coef(summary(OLR)))

### COMPLETED LOGISTIC REGRESSION ###


########################################################## 
########### ADJUSTMENT MODEL - Calibrated Data ###########
########################################################## 

# Create empty vectors for for-loop input
data = as.data.frame(calib_data)
data$test1 = as.numeric(as.character(data$test1))
data$allP = as.numeric(as.character(data$allP))
adj_val = c()
allP = c()
startq = c()
ratio =data$allP/data$test1
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allP[3*i - 2], data$allP[3*i - 1], data$allP[3*i])
  t_y <- c(data$test1[3*i - 2], data$test1[3*i - 1], data$test1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}
adjusted_test1 <- test1 + adj_val 
# Append adjusted test1 values and adjustment value to data set
calib_data=cbind(data,adjusted_test1,adj_val)
# Write Calibrated Data CSV --> Used in "qPCR_Plotting" code for visuals
#write.csv(file="YEAR_MONTH_Calibrated_DF", calib_data)


### LINEAR MODEL ### --> NO LONGER BEING USED
# Creating the adjustment model lm(y-axis~x-axis)
# Changed adj_val^2 to adj_val to try to make the model better --> VERIFY THIS WITH DR. WANG
#adj_model <- lm(adj_val~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
#summary(adj_model)
#par(mfrow = c(2,2))
#plot(adj_model)
#dev.off()
#scatter.smooth(ratio, adj_val)
#abline(ratio, adj_val)

### COMPLETED ADJUSTMENT MODEL - CALIBRATED DATA ###


########################################################## 
################ Experimental Data Framing ###############
########################################################## 
# Create/Write data frame for Experimental values
exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_df = exp_df[order(exp_df$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
exp_df = exp_df[-1,]
exp_df = exp_df[-nrow(exp_df),]
exp_df$sampleID = as.numeric(as.character(exp_df$sampleID))
exp_df$cpD1 = as.numeric(as.character(exp_df$cpD1))
exp_data = exp_df
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
countsne2$Var1 = as.numeric(as.character(countsne2$Var1)) #---> CHECK IF THIS IS NECESSARY
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]
### Report invalid observations ###
# Send CSV of removed sampleID's to Dr. S (invalid obs), with additional plots of raw cycle values for invalid obs
# Write CSV file to send Dr. S for investigation
### WORK ON --> add derivative values in to the CSV file
### WORK ON --> creating a separate CSV file with samples with unusual derivatives
#write.csv(file="2018_11_SamplesToInvestigate", countsne2)
#write.csv(file="YEAR_MONTH_SamplesToInvestigate", countsne2)

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
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$sampleID.exp = as.numeric(as.character(exp_data$sampleID.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))
# Calculate ratios for experimental data 
ratio.exp = exp_data$allP.exp/exp_data$test1.exp
# Append ratios to data set
exp_data = cbind(exp_data, ratio.exp)
# Write Experimental Data CSV --> Used in "qPCR_Plotting" code for visuals
#write.csv(file="YEAR_MONTH_Experimental_DF", exp_data)

### COMPLETED EXPERIMENTAL DATA FRAME ###

########################################################## 
########## ADJUSTMENT MODEL - Experimental Data ##########
########################################################## 

# Using the adjustment model on the expiremental data
new = data.frame(ratio = exp_data$allP/exp_data$test1)
predict(adj_model, new , interval = "confidence")
#---> Fill in any remaining parts of experimental adjusmtent model

### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###