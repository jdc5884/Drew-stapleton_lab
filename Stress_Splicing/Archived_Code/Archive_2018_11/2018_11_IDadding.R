## adding unique sample id numbers to 2018_11 plate data

dat = read.csv(file = "2018_11_withStress.csv", header = TRUE)
sampPlan = read.csv(file = "../2016_Clayton/Field Book (2016) - Clayton - Sampling Plan_TIDIED.csv", header = TRUE)

colnames(dat)[2] = "sampleID"
colnames(sampPlan)[3] = "sampleID"


full = merge(sampPlan,dat, by = "sampleID")
full = full[, -c(2, 13, 14, 15)]
write.csv(full, file = "2018_11_ID&Stress.csv")
