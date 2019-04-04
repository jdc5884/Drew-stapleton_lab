########################################################## 
############### QPCR PLATE & MODEL PLOTTING ##############
########################################################## 

### READ IN CALIBRATED & EXPERIMENTAL DATA FRAMES (OUTPUT FROM "qPCR_Adjustment_YEAR_MO" code) ###
# Read in Calibrated Data Frame
data=read.csv(file="YEAR_MONTH_Calibrated_DF", header=TRUE)
# Read in Experimental Data Frame
exp_data=read.csv(file="YEAR_MONTH_Experimental_DF", header=TRUE)

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