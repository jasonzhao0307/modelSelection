require(caret)
require(mlbench)
require(DESeq2)
require(ggplot2)
require(plotROC)
require(caretEnsemble)
require(pROC)
library("caTools")

############### functions define ##############

# get threshold of prediction probability for calculating best sensitivit + specificity
Get_Classification_Threshold <- function(predict, response) {
    perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
    df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
    df[which.max(df$sens + df$spec), "cut"]
}

# input a vector of predicted probability and true label, output sensitivity and specificity at best cut-off.
Get_Sensitivity_And_Specificity <- function(predict, response, provide.cutoff = FALSE,
                                            cutoff.value = 0.5,
                                            prevalence){
  if (provide.cutoff == FALSE){
    prob.cutoff <- Get_Classification_Threshold(predict, response)
  } else{
    prob.cutoff <- cutoff.value
  }

  pred <- ifelse(predict > prob.cutoff, 1, 0)
  # Must specify which class is positive !!!
  # also, in this case, we need to use character.
  # provide prevalence!!
  cm = confusionMatrix(as.character(pred), as.character(response), positive = "1", prevalence = prevalence)
  sensi <- cm$byClass[1]
  speci <- cm$byClass[2]
  ppv <- cm$byClass[3]
  npv <- cm$byClass[4]

  return(list(sens=sensi,spec=speci,ppv=ppv,npv=npv))
}


MaxFilter <- function(df, max.value = 10){
  df.filtered <- df[which(apply(df,1,max) >= max.value),]
  return(df.filtered)
}


deseq2_norm_rle <- function(dataFrame){
# RLE normalization: relative log expression
scaling.dataFrame <- estimateSizeFactorsForMatrix(dataFrame)
dataFrame.scaled <- dataFrame
for(i in 1:ncol(dataFrame)){
  dataFrame.scaled[,i] <- dataFrame[,i]/scaling.dataFrame[i]
}

return(dataFrame.scaled)
}



### centering and scaling data
normData <- function(training, testing){
  preProcValues <- preProcess(training, method = c("center", "scale"))
  trainTransformed <- predict(preProcValues, training)
  testTransformed <- predict(preProcValues, testing)
  return(list(trainTransformed = trainTransformed,
              testTransformed = testTransformed))
}



# prepare the dataset: user needs to provide which class in y is positive
prepareData <- function(df.training, df.testing, y.train, y.test, positive){
  y.train.tmp <- as.character(y.train)
  y.test.tmp <- as.character(y.test)
  y.train <- rep("N", length(y.train))
  y.test <- rep("N", length(y.test))
  y.train[y.train.tmp == as.character(positive)] <- "P"
  y.test[y.test.tmp == as.character(positive)] <- "P"
  training <- cbind(df.training, y.train)
  testing <- cbind(df.testing, y.test)
  training$y.train <- as.factor(training$y.train)
  testing$y.test <- as.factor(testing$y.test)
  colnames(training)[ncol(training)] <- "Class"
  colnames(testing)[ncol(testing)] <- "Class"
  return(list(training = training, testing = testing))
}


### Final model selection function
modelSelection <- function(df.training,
                           model.names,
                           num.cv.fold = 10,
                           num.cv.repeat = 10,
                           num.param.tune = 12,
                           seed.use = 100,
                           show.some.models = FALSE){
# const define
model.names.example <- c("lda", "naive_bayes", "gbm", "glmnet", "ranger", "svmLinear",
                 "svmRadial", "xgbLinear", "xgbTree")

if (show.some.models == TRUE){
  print(paste0("Here are some models you could try: ", paste(model.names.example,
                                                             collapse = ",")))
}
# Model training general parameters setting
fitControl <- trainControl(method = "repeatedcv",
                           number = num.cv.fold,
                           repeats = num.cv.repeat,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           savePredictions = TRUE,
                           ## parallel in windows:
                           allowParallel = TRUE,
                           ## Evaluate performance using
                           ## the following function
                           summaryFunction = twoClassSummary)

# train multiple models
list.model.fit <- list()
for (i in 1:length(model.names)){
  print(paste0("Now training using: ", model.names[i]))
 list.model.fit[[i]] <- caretModelTraining(df.training,
                                           model.names[i],
                                           fitControl,
                                           num.param.tune,
                                           seed.use)
}
names(list.model.fit) <- model.names

# add ensemble model



tmp.list <- list.model.fit
class(tmp.list) <- "caretList"

### correlation between models
roc.scatter.plot <- xyplot(resamples(tmp.list))
roc.cor <- modelCor(resamples(tmp.list))

set.seed(seed.use)
greedy_ensemble <- caretEnsemble(
  tmp.list, 
  metric="ROC",
  trControl=fitControl)


list.model.fit[[(length(list.model.fit) + 1)]] <- greedy_ensemble$ens_model
names(list.model.fit)[length(list.model.fit)] <- "Ensemble"



# compare models using their best params
resamps <- resamples(list.model.fit)
#print(summary(resamps))
model.summary <- (summary(resamps))$statistics
model.summary <- lapply(model.summary, function(x) {x[,3:4]})
df.model.summary <- as.data.frame(model.summary[[1]])
for (i in 2:length(model.summary)){
  df.model.summary <- cbind(df.model.summary, model.summary[[i]])
}
colnames(df.model.summary) <- c("AUC-median", "AUC-mean",
                                "Sens-median", "Sens-mean",
                                "Spec-median", "Spec-mean")


# Visualize the model comparison in AUC, sens and spec
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
plot.handle <- bwplot(resamps, layout = c(3, 1))


### find the best model by finding the largest median AUC
model.best.index <- which(df.model.summary[,1] ==
                            max(df.model.summary[,1]))
model.best <- list.model.fit[[model.best.index]]


print("Finished.")

return(list(list.model.fit = list.model.fit,
            df.model.summary = df.model.summary,
            comparison.plot = plot.handle,
            model.best = model.best,
            model.ensemble = greedy_ensemble,
            roc.scatter.plot = roc.scatter.plot,
            roc.cor = roc.cor))
}



### Individual model training
# train classification model, with the training matrix being:
# dim = m*n, m = #sample, n = (#features + 1)
# The last column is the label factor with colname: "Class"
caretModelTraining <- function(df.training, model.name, fitControl, num.param.tune, seed.use){
  set.seed(seed.use)
  if (model.name %in% c("glm", "glmnet")){
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "ROC")
  } else{
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "ROC",
                 verbose = FALSE)
  }

  return(modelFit)
}


### predict and evaluate functions
predictAndEvaluation <- function(model.best,
                                 test.data,
                                 prevalence,
                                 is.ensemble = FALSE){
    # create a numerical class var
    class.num <- as.character(test.data$Class)
    class.num[class.num == "P"] <- 1
    class.num[class.num == "N"] <- 0
    class.num <- as.numeric(class.num)

    ## predict
    test.prediction <- predict(model.best, newdata = test.data)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")

    ### ROC and AUC
    if (is.ensemble == FALSE){
        g <- ggplot(test.prediction.prob, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))
    
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob[,2]
    } else{
        prob.df <- data.frame(P = 1-test.prediction.prob, N = test.prediction.prob)
        g <- ggplot(prob.df, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))
    
        #### get probabilities for POSITIVE
        model.prob = 1-test.prediction.prob
    }


    ## get sens and spec, ppv and npv (require prevalence!)
    test.prediction.evaluation <- Get_Sensitivity_And_Specificity(predict = model.prob,
                                    response = class.num,
                                    prevalence = 0.1)

    ## combine all metrics
    test.prediction.evaluation$auc <- auc.value
    evaluation.names <- names(test.prediction.evaluation)
    test.prediction.evaluation <- sapply(test.prediction.evaluation, function(x) x)
    names(test.prediction.evaluation) <- evaluation.names

    ## return the results
    return(list(positive.prob = model.prob,
                test.prediction.evaluation = test.prediction.evaluation,
                roc = g))

}
