# qPCR Adjustment

library(tidyr)

# Mac Directory
setwd("/Users/andrewnorris/stapleton_lab/Stress_Splicing/2018_6")

# PC Directory
#setwd()

# Read in Calibration Data
data <- read.csv(file = "2018_6_1_Calibrated_Data_Frame_with_Derivatives.csv")
# Format starting quantity as numeric, not in scientific notation
options(scipen=5)
# Reorder columns with starting quantity first, remove unneeded ordering column (how can I better word this...)
data <- data[c(4,2,3,5,6)]
# Indicate if sample is 'Plus' or 'Minus'
data['sampleID_Plus'] <- grepl('plus', data$sampleID)

# Read in Expiremental Data
exp_data <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv")     
# Format starting quantity as numeric
exp_data <- exp_data[c(4,2,3,5,6)]
# Remove first extra labeling row
exp_data <- exp_data[-1,]

## If zero matrix is the way to go... #Create zero matrix to which data will input
# relevant = data.frame(matrix(rep(0,length(dat$sampleID)*dim(exp_data)[1]), ncol = dim(exp_data)[1]))
##test1 <- data$reaction_type("test1")

# cp Values of Calibration data
test1 <- data$cpD1
all_prod <- data$all_productsCp1
# The ratio of the cp values
ratio <- all_prod/test1
# Appending to the Calibration data
cbind(data, ratio)
# Adjustment value Vector
adj_val = vector()


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
