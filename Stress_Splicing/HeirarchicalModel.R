### HEIARCHICAL MODEL ###

library(stringr)
library(dplyr)
library(pracma)

# MONTH 1 (2018_6 / JUNE) CALIBRATED DATA FRAME
calib_data_6 = read.csv("./2018_6/calib_2018_6.csv")[,-1]
calib_data_6$month ='june'
calib_data_6

# MONTH 2 (2018_8 / AUGUST) CALIBRATED DATA FRAME
calib_data_8 = read.csv("./2018_8/calib_2018_8.csv")[,-1]
calib_data_8$month ='aug'
calib_data_8

# MONTH 3 (2018_11 / NOVEMBER) CALIBRATED DATA FRAME
calib_data_11 = read.csv("./2018_11/calib_2018_11.csv")[,-1]
calib_data_11$month ='nov'
calib_data_11


# Combined Calib d.f. for all months
calib_data = rbind(calib_data_6, calib_data_8, calib_data_11)

# Create dummy varible columns for each month
calib_data$june = ifelse(str_detect(calib_data[,4], "june"), 1, 0)
calib_data$aug = ifelse(str_detect(calib_data[,4], "aug"), 1, 0)
calib_data = calib_data[,-4]

# Ordinal Logistic Regression Model 
model = polr(as.factor(calib_data$startq) ~ ., data=calib_data, Hess = TRUE)

(ctable <- coef(summary(model)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
options(scipen=999)
## combined table
(ctable <- cbind(ctable, "p value" = p))

# Linear Model
lin_model = lm(as.factor(calib_data$startq) ~ ., data=calib_data)
lin_model





# #### Month 1 (2018_6) CT ####
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
# m11 = read.csv(file = "2018_6_1_qPCR_Output.csv", header=FALSE)
# #Initial framing
# m11 = m11[,-1]
# # Transpose derivatives to be in equivalent format as raw plate data
# m11 = as.data.frame(t(m11), header=FALSE)
# m11 = m11[,-6]
# # Rename columns
# colnames(m11)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1")
# # Remove NTC, gblock_minus values
# m11['sampleID_NTC'] = grepl('NTC', m11$sampleID)
# ntc = which(m11$sampleID_NTC)
# m11 = m11[-ntc,]
# m11 = m11[,-6]
# m11['sampleID_Minus'] = grepl('minus', m11$sampleID)
# minus = which(m11$sampleID_Minus)
# m11 = m11[-minus,]
# m11 = m11[,-6]
# m11$cpD1 = as.numeric(as.character(m11$cpD1))
# # Remove CT vals <10
# m11 = m11 %>% filter(m11$cpD1 >= 10)
# # Create calibrated data frame
# calib_data = m11 %>% filter(str_detect(sampleID, "g"))

# 
# #Month 2 (2018_8) CT
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_8")
# m21 = read.csv(file = "2018_8_1_qPCR_Output.csv", header=FALSE)
# m21 = m21[,-1]
# m22 = read.csv(file = "2018_8_2_qPCR_Output.csv", header=FALSE)
# m22 = m22[,-1]
# m23 = read.csv(file = "2018_8_3_qPCR_Output.csv", header=FALSE)
# m23 = m23[,-1]
# 
# #Month 3 (2018_11) CT
# setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
# m31 = read.csv(file = "2018_11_1_qPCR_Output.csv", header=FALSE)
# m31 = m31[,-1]
# m32 = read.csv(file = "2018_11_2_qPCR_Output.csv", header=FALSE)
# m32 = m32[,-1]
# 
# #Combine months
# month = data.frame(m11, m21, m22, m23, m31, m32)
# month = as.data.frame(t(month))
# month$june = ifelse(str_detect(month[,1], "2018_6"), 1, 0)
# month$V1 = str_replace(month$V1, "2018_9", "2018_8")
# month$aug = ifelse(str_detect(month[,1], "2018_8"), 1, 0)
# 
# #Remove CT<10, NTC, gblock_minus, V6 (CT2)
# month$V5 = as.numeric(as.character(month$V5))
# month$V6 = as.numeric(as.character(month$V6))
# month = month %>% filter(month$V5 >= 10)
# month = month %>% filter(month$V6 >= 10)
# month['NTC'] = grepl('NTC', month$V3)
# ntc = which(month$NTC)
# month = month[-ntc,]
# month = month[,-9]
# month['Minus'] = grepl('minus', month$V3)
# minus = which(month$Minus)
# month = month[-minus,]
# month = month[,-c(6,9)]

