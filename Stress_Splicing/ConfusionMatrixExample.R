library('caret')
numLlvs <- 4
confusionMatrix(
  factor(sample(rep(letters[1:numLlvs], 200), 50)),
  factor(sample(rep(letters[1:numLlvs], 200), 50)))