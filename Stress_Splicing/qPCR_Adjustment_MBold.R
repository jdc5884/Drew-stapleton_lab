# qPCR Adjustment using U-Statistics
# Author: Michael Byrd

# Mac Directory
setwd("/Users/mbyrd/StapletonLab/Stapleton_Lab/qPCR/adjustment")

# Windows Directory
# setwd("C:/Users/mbyrd/Documents/Stapleton/Stapleton_Lab/qPCR/adjustment")


# Calibration Data
data <- read.csv(file = "Neill2018calibrations.csv")
# Expiremental Data
exp_data <- read.csv(file = "corrected17June_Neill Thesis RNA samples.csv")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
# cp Values of Calibration data
test1 <- data$test1_Cp1
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
