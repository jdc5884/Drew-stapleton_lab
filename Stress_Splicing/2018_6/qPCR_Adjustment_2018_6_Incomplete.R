########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 


library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")
# PC Directory
#setwd(~/Desktop/GIThub/StapletonLab/StressSplicing/qPCR/)

### READ IN DERIVATIVE DATA###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
#deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
#deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
#deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
deriv=read.csv(file = "2018_6_1_qPCR_output_withHeaders.csv", header=FALSE)

########################################################## 
################### Initial Data Framing #################
########################################################## 

# Remove extra labels row and column 
deriv = deriv[-c(1,5),-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Rename columns
colnames(deriv)=c("reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# Remove extra labels row
deriv=deriv[-1,]
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
# IF "minus" RETURNS EMPTY VALUES, COMMENT OUT COMMAND BELOW
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

### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############ ADJUSTMENT MODEL Calibrated Data ############
########################################################## 

# Create empty vectors for for-loop input
data = as.data.frame(calib_data)
data$test1 = as.numeric(as.character(data$test1))
data$allP = as.numeric(as.character(data$allP))
adj_val = c()
allP = c()
startq = c()
# Calculate ratio values
ratio =data$allP/data$test1
# Append ratios to data set
data=cbind(data,ratio)
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allP[3*i - 2], data$allP[3*i - 1], data$allP[3*i])
  t_y <- c(data$test1[3*i - 2], data$test1[3*i - 1], data$test1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}
# Calculate adjusted test1 value
adjusted_test1 <- test1 + adj_val
# Append adjusted test1 values and adjustment value to data set
calib_data=cbind(data,adjusted_test1,adj_val)
# Creating the adjustment model lm(y-axis~x-axis)
# Changed adj_val^2 to adj_val to try to make the model better --> VERIFY THIS WITH DR. WANG
adj_model <- lm(adj_val~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
summary(adj_model) 
# Plot adjustment model 
par(mfrow = c(2,2))
plot(adj_model)
dev.off()
# Plot ratios vs. adjustment values
plot(ratio, adj_val, main="Ratios vs. Adjustment Values", col="blue", xlab="Ratio", ylab="Adjustment Value")
abline(lm(adj_val~ratio), col="blue")

### COMPLETED ADJUSTMENT MODEL - CALIBRATED DATA ###

########################################################## 
################ Experimental Data Framing ###############
########################################################## 
# Create/Write data frame for Calibrated values
exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_df = exp_df[order(exp_df$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
#exp_df = exp_df[-1,]
#exp_df = exp_df[-nrow(exp_df),]
#exp_df$sampleID = as.numeric(as.character(exp_df$sampleID))
#exp_df$cpD1 = as.numeric(as.character(exp_df$cpD1))
exp_data = exp_df
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
# Find counts of each unique sampleID; for sample with a count not equal to 2
# Send removed sampleID's to Dr. S, with additional plots of raw cycle values for each of these samples
# Write CSV file of samples with count not equal to 2 to send to Dr. S for investigation
### To work on --> add derivative values in to the CSV file
### To work on --> creating a separate CSV file with samples with unusual derivatives
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
#write.csv(file="2018_11_SamplesToInvestigate", countsne2)
# Remove samples without both allP and test1 from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

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
exp_data = as.data.frame(exp_data)
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))
# Create column of ratios
ratio.exp = exp_data$allP.exp/exp_data$test1.exp
# Append ratio column to data frame
exp_data = cbind(exp_data, ratio.exp)
exp_data$ratio.exp = as.numeric(as.character(exp_data$ratio.exp))

# Write CSV file
#write.csv(exp_data, file="2018_11_Experimental_Data_Frame.csv")

### COMPLETED EXPERIMENTAL DATA FRAME ###

### ADJUSTMENT MODEL - EXPERIMENTAL DATA ### 
# Using the adjustment model on the expiremental data
exp.data = data.frame(ratio = exp_data$allP.exp/exp_data$test1.exp)

predict(adj_model, exp.data, interval = "confidence")

#---> Fill in any remaining parts of experimental adjusmtent model

### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###

##PLOTS##
#AllP#
hist(data$allP, xlim=c(0,50), ylim=c(0,100), col=rgb(1,0,0,0.5), main='Histogram of All Products', xlab='All Products Derivative')
hist(exp_data$allP.exp, xlim=c(0,50), ylim=c(0,100), add=T, col=rgb(0,0,1,0.5))
legend("topleft",
        c("Calibration", "Experimental"),
          fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
#dev.off()
#Test1#
hist(data$test1, xlim=c(0,30), ylim=c(0,80), col=rgb(1,0,0,0.5), main='Histogram of Test 1', xlab='Test 1 Derivative')
hist(exp_data$test1.exp, xlim=c(0,30), ylim=c(0,80), add=T, col=rgb(0,0,1,0.5))
legend("topleft",
       c("Calibration", "Experimental"),
       fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
#dev.off()
#Ratios - Calibrated#
hist(data$ratio, xlim=c(0,3), ylim=c(0,70), col=rgb(1,0,0,0.5), main='Histogram of Ratios', xlab='Ratio')
#Ratios - Experimental#
# Values excluded from histogram that will be further investigated later (their index)
x=c(8,82,141,148,149,153,161,170,172,175,180,188)
exp_data2=exp_data[-x,]
hist(exp_data2$ratio.exp, xlim=c(0,3), ylim=c(0,70), col=rgb(0,0,1,0.5), add=T)
legend("topleft",
       c("Calibration", "Experimental"),
       fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
#Calib Plot - S.Q. vs. Ratios
plot(calib_data$startq, calib_data$ratio, xlab='Starting Quantity', ylab='Ratio', 
      main='Calibrated Data - Starting Quantities vs. Ratios')
#Calib Plot - Test1 vs. Ratio
plot(calib_data$test1, calib_data$ratio, xlab='Test 1 Derivative', ylab='Ratio', 
      main='Calibrated Data - Test 1 Derivative vs. Ratio')
