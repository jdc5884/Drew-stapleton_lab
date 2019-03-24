# Visual Representations of Calibrated and Experimental Data Frames

setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_6")
library(stringr)
library(tidyverse)

## Calibrated Data Frame
par(mfrow = c(2,1))
# Read in Calibration Data
calib.data <- read.csv(file = "2018_6_1_Calibrated_Data_Frame_with_Derivatives.csv")
# Format starting quantity as numeric, not in scientific notation
options(scipen=5)
# Reorder columns with starting quantity first, remove unneeded ordering column
calib.data <- calib.data[c(4,2,3,5,6)]
calib.data <- calib.data[order(calib.data$starting_quantity),]
# Sort data by reaction type (all_products / test1), ascending starting quantities
calib.data <- calib.data[order(calib.data$reaction_type),]
# Create new column to indicate gblock_minus values
calib.data['sampleID_minus'] <- grepl('minus', calib.data$sampleID)
# Remove gblock_minus values
minus <- which(calib.data$sampleID_minus)
calib.data = calib.data[-minus,]
# Take natural log of starting quantities to linearly fit to calibration data
sq.ln <- log10(calib.data$starting_quantity)
# Bind natural log column to calib.data frame
calib.data <- cbind(sq.ln, calib.data)


# Separate by reaction type  
all_prod=calib.data[c(1:27),] 
test1=calib.data[c(28:54),]
# Boxplot
calib.plot.all_prod <- boxplot(cpD1~sq.ln, data=all_prod, xlab="Nat. Log of Starting Quantity", ylab="Derivative", 
                      main="Calibration Derivatives by All Products Starting Quantity", 
                      col=c("darkgreen"), ylim=c(0,34))
calib.plot.test1 <- boxplot(cpD1~sq.ln, data=test1, xlab="Nat. Log of Starting Quantity", ylab="Derivative", 
                               main="Calibration Derivatives by Test 1 Starting Quantity", 
                               col=c("blue"), ylim=c(0,34))

#calib.plot.all_prod <- boxplot(cpD1~starting_quantity, data=all_prod, xlab="Starting Quantity", ylab="Derivative", 
                               main="Calibration Derivatives by All Products Starting Quantity", 
                               col=c("darkgreen"), ylim=c(0,34))
#calib.plot.test1 <- boxplot(cpD1~starting_quantity, data=test1, xlab="Starting Quantity", ylab="Derivative", 
                            main="Calibration Derivatives by Test 1 Starting Quantity", 
                            col=c("blue"), ylim=c(0,34))

dev.off()

## Experimental Data Frame
# Read in Expiremental Data
exp.data <- read.csv(file = "2018_6_1_Experimental_Data_Frame_with_Derivatives.csv")     
# Format starting quantity as numeric
exp.data <- exp.data[c(4,2,3,5,6)]
# Remove first extra labeling row
exp.data <- exp.data[-1,]
# Order data by reaction type
exp.data <- exp.data[order(exp.data$reaction_type),]
# Separate by reaction type  
all_prod_exp=exp.data[c(1:99),] 
test1_exp=exp.data[c(100:195),]


# Change cpD1-all_products column into vector for use in dotplot
cpD1_all_products <- as.integer(all_prod_exp$cpD1)
# Change cpD1-test1 column into vector for use in dotplot
cpD1_test1 <- as.integer(test1_exp$cpD1)

# Histogram
par(mfrow = c(2,1))
hist(cpD1_all_products, xlim = c(0,120), ylim = c(0,80))
hist(cpD1_test1, xlim = c(0,120), ylim = c(0,80))
dev.off()

# Dotplot
dotchart(cpD1_all_products,labels=all_prod$sampleID,
         main="Experimental Derivatives by Sample ID and All Products Reaction Type",
         xlab="Derivative", xlim=c(0,120))
dotchart(cpD1_test1,labels=test1$sampleID,
         main="Experimental Derivatives by Sample ID and Test1 Reaction Type",
         xlab="Derivative", xlim=c(0,120))
###^^^ Need to adjust y-axes to clearly differentiate sampleID --- Histogram seems more useful than dotplot
