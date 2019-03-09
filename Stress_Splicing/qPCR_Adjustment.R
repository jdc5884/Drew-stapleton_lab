# qPCR Adjustment Model

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)

# Mac Directory
setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")

# PC Directory
#setwd()

# Read in Calibration Data
calib_data <- read.csv(file = "2018_6_1_Calibrated_Data_Frame_with_Derivatives.csv")
# Format starting quantity as numeric, not in scientific notation
options(scipen=5)
# Reorder columns with starting quantity first, remove unneeded ordering column (how can I better word this...)
calib_data <- calib_data[c(4,2,3,5,6)]
# Indicate if sample is 'Plus' or 'Minus'
calib_data['sampleID_Plus'] <- grepl('plus', calib_data$sampleID)
# Remove 'Minus' values 
bad <- which(!calib_data$sampleID_Plus)
calib_data = calib_data[-bad,]
# Sort data by starting quantity
calib_data = calib_data[order(calib_data$starting_quantity),]
# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_data$starting_quantity)){
  sq <- calib_data$starting_quantity[i]
  if(i %% 6 == 1){
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
c = as.data.frame(startq)
calib.df = cbind(c, test1, allP)

# Read in Expiremental Data
exp_data <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv", head=TRUE)     
# Format starting quantity as numeric
exp_data <- exp_data[c(4,2,3,5,6)]
# Remove first extra labeling row
#exp_data <- exp_data[-1,]
# Sort data by sample ID
exp_data <- exp_data[order(exp_data$sampleID),]
# Remove NTC values
exp_data <- exp_data[-(187:196),]

#exp_data <- as.double(exp_data$cpD1)
    ### ERROR:  After converting values from integers to double (req'd to complete exp.df), 
      # error returned "$ operator is invalid for atomic vectors"  

    ###ERROR in exp_data:  Will not read values as double, not integers, even though values are being pulled from 
      # equivalent-types CSV files 

# Create two data frames from reaction type (test1 or all_products) 
test1.exp_data = exp_data %>% filter(str_detect(exp_data$reaction_type, "test1") == TRUE)
allP.exp_data = exp_data %>% filter(str_detect(exp_data$reaction_type, "all_products") == TRUE)
# Merge two data frames by sample ID
#exp.df <- merge(test1.exp_data, allP.exp_data, by="sampleID")



# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp <- toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp <- c(sampleID.exp, id.exp)
  }
  val <- toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp <- c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
#    print(exp_data$cpD1[i])
    allP.exp <- c(allP.exp, exp_data$cpD1[i])
  }
}
# Transpose test1.exp
test1.exp <- t(test1.exp)
# Create allP.exp as a data frame; Transpose to single column
allP.exp <- as.data.frame(allP.exp)
allP.exp <- t(allP.exp)
# Make startq.exp as data frame
sampleID.exp <- as.data.frame(sampleID.exp)
#e <- cbind(startq.exp)
exp.df <- cbind(sampleID.exp, test1.exp, allP.exp)



####________________________________________###
### Old / Unused Code ###

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


# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$all_productsCp1[3*i - 2], data$all_productsCp1[3*i - 1], data$all_productsCp1[3*i])
  t_y <- c(data$test1_Cp1[3*i - 2], data$test1_Cp1[3*i - 1], data$test1_Cp1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}

adjusted_test1 <- test1 + adj_val

# Creating the adjustment model lm(y-axis~x-axis)
adj_model <- lm(adj_val^2~ratio)
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
adj_model_50 <- lm(adj_val_50_50^2~ratio_50)

summary(adj_model_50)
plot(adj_model_50)

# Unused old code:

## If zero matrix is the way to go... #Create zero matrix to which data will input
# relevant = data.frame(matrix(rep(0,length(dat$sampleID)*dim(exp_data)[1]), ncol = dim(exp_data)[1]))
##test1 <- data$reaction_type("test1")