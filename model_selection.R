#########################
# Author: Yue (Jason) Zhao
# Github: https://github.com/jasonzhao0307




############### Load packages ##############

require(caret)
require(mice)
require(mlbench)
require(DESeq2)
require(ggplot2)
require(plotROC)
require(caretEnsemble)
require(pROC)
library(caTools)
library(gmodels)

############### functions define ##############

# get threshold of prediction probability for calculating best sensitivit + specificity
Get_Classification_Threshold <- function(predict, response) {
    perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
    df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@y.values[[1]], spec = perf@x.values[[1]])
    df[which.max(df$sens + df$spec), "cut"]
}



### Sensitivity Maximize threshold at 0.9
Get_Classification_Threshold_Sen <- function(predict, response, sen_min = 0.9, spe_min = 0.5) {
    perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
    df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@y.values[[1]], spec = perf@x.values[[1]])
    df <- df[-1,]
    df_tmp <- df[which(df$sens >= sen_min),]
    df_tmp <- df_tmp[which(df_tmp$spec >= spe_min),]
    if (nrow(df_tmp) == 0){
        print("Can't meet both sensitivity and specificity target! Will use max(sen+spe) instead")
        cut_off_report <- df[which.max(df$sens + df$spec), "cut"]
    }else{
        cut_off_report <- df_tmp[which.max(df_tmp$sens + df_tmp$spec), "cut"]
    }
    cut_off_report

}



# input a vector of predicted probability and true label, output sensitivity and specificity at best cut-off.
Get_Sensitivity_And_Specificity <- function(predict, response, provide.cutoff = FALSE,
                                            cutoff.value = 0.5,
                                            prevalence = 0.1){
  if (provide.cutoff == FALSE){
    prob.cutoff <- Get_Classification_Threshold(predict, response)
  } else{
    prob.cutoff <- cutoff.value
  }

  pred <- ifelse(predict >= prob.cutoff, 1, 0)
  # Must specify which class is positive !!!
  # also, in this case, we need to use character.
  # provide prevalence!!
  cm = confusionMatrix(as.factor(pred), as.factor(response), positive = "1", prevalence = prevalence)
  sensi <- cm$byClass[1]
  speci <- cm$byClass[2]
  ppv <- cm$byClass[3]
  npv <- cm$byClass[4]

  return(list(sens=sensi,spec=speci,ppv=ppv,npv=npv))
}


# input a vector of predicted probability and true label, output sensitivity and specificity at best cut-off.
Get_Sensitivity_And_Specificity_Sen <- function(predict, response, provide.cutoff = FALSE,
                                            cutoff.value = 0.5,
                                            prevalence = 0.1,
                                            sen_min = sen_min,
                                            spe_min = spe_min){
  if (provide.cutoff == FALSE){
    prob.cutoff <- Get_Classification_Threshold_Sen(predict, response, sen_min, spe_min)
  } else{
    prob.cutoff <- cutoff.value
  }

  pred <- as.character(ifelse(predict >= prob.cutoff, 1, 0))
  # Must specify which class is positive !!!
  # also, in this case, we need to use character.
  # provide prevalence!!
  cm = confusionMatrix(as.factor(pred), as.factor(response), positive = "1", prevalence = prevalence)
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
  if (length(unique(y.train)) == 2  & sum(unique(y.train) %in% c("P","N")) == 2){
    training <- cbind(df.training, y.train)
    testing <- cbind(df.testing, y.test)
  }else{
    y.train.tmp <- as.character(y.train)
    y.test.tmp <- as.character(y.test)
    y.train <- rep("N", length(y.train))
    y.test <- rep("N", length(y.test))
    y.train[y.train.tmp == as.character(positive)] <- "P"
    y.test[y.test.tmp == as.character(positive)] <- "P"
    training <- cbind(df.training, y.train)
    testing <- cbind(df.testing, y.test)
  }

  #training$y.train <- as.factor(training$y.train)
  #testing$y.test <- as.factor(testing$y.test)
  training$y.train <- as.character(training$y.train)
  testing$y.test <- as.character(testing$y.test)
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




print("Finished.")

return(list(list.model.fit = list.model.fit,
            df.model.summary = df.model.summary,
            comparison.plot = plot.handle,
            model.ensemble = greedy_ensemble,
            roc.cor = roc.cor))
}



### Final model selection function, maximizing sensitivity
modelSelectionSen <- function(df.training,
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
 list.model.fit[[i]] <- caretModelTrainingSen(df.training,
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
roc.cor <- modelCor(resamples(tmp.list))

set.seed(seed.use)
greedy_ensemble <- caretEnsemble(
  tmp.list,
  metric="Sens",
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




print("Finished.")

return(list(list.model.fit = list.model.fit,
            df.model.summary = df.model.summary,
            comparison.plot = plot.handle,
            model.ensemble = greedy_ensemble,
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




### Individual model training
# train classification model, with the training matrix being:
# dim = m*n, m = #sample, n = (#features + 1)
# The last column is the label factor with colname: "Class"
caretModelTrainingSen <- function(df.training, model.name, fitControl, num.param.tune, seed.use){
  set.seed(seed.use)
  if (model.name %in% c("glm", "glmnet")){
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "Sens",
                 maximize = TRUE)
  } else{
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "Sens",
                 maximize = TRUE,
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
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
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
        prob.df <- data.frame(P = test.prediction.prob, N = 1-test.prediction.prob)
        g <- ggplot(prob.df, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob
        if (auc.value < 0.5){
          prob.df <- data.frame(P = 1-test.prediction.prob, N = test.prediction.prob)
          g <- ggplot(prob.df, aes(m=P, d=class.num)) +
          geom_roc(n.cuts=0) +
          coord_equal() +
          style_roc()
          auc.value <- round((calc_auc(g))$AUC, 4)
          #### get probabilities for POSITIVE
          model.prob = 1-test.prediction.prob
        }
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))


    }

    #return(list(predict = model.prob, response = class.num))
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



### predict and evaluate functions
predictAndEvaluationSen <- function(model.best,
                                 test.data,
                                 prevalence,
                                 is.ensemble = FALSE,
                                 sen_min = 0.9,
                                 spe_min = 0.2){
    # create a numerical class var
    class.num <- as.character(test.data$Class)
    class.num[class.num == "P"] <- 1
    class.num[class.num == "N"] <- 0
    class.num <- as.numeric(class.num)

    ## predict
    test.prediction <- predict(model.best, newdata = test.data)
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
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
        prob.df <- data.frame(P = test.prediction.prob, N = 1-test.prediction.prob)
        g <- ggplot(prob.df, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob
        if (auc.value < 0.5){
          prob.df <- data.frame(P = 1-test.prediction.prob, N = test.prediction.prob)
          g <- ggplot(prob.df, aes(m=P, d=class.num)) +
          geom_roc(n.cuts=0) +
          coord_equal() +
          style_roc()
          auc.value <- round((calc_auc(g))$AUC, 4)
          #### get probabilities for POSITIVE
          model.prob = 1-test.prediction.prob
        }
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))


    }

    #return(list(predict = model.prob, response = class.num))
    ## get sens and spec, ppv and npv (require prevalence!)
    test.prediction.evaluation <- Get_Sensitivity_And_Specificity_Sen(predict = model.prob,
                                    response = class.num,
                                    prevalence = 0.1,
                                    sen_min = sen_min,
                                    spe_min = spe_min)

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


predictAndEvaluationSenBootstrap <- function(model_use,
                                             data_test,
                                             sen_min,
                                             spe_min,
                                             seed_use = 1,
                                             n_boot = 50,
                                             prevalence = 0.1
                                             ){
set.seed(seed_use)
sen.vec <- c()
spe.vec <- c()
ppv.vec <- c()
npv.vec <- c()
auc.vec <- c()

for (i in 1:n_boot){
    index <- sample(1:nrow(data_test),size = nrow(data_test),replace = TRUE)
    print(paste0("Round: ", i))
    result_tmp <- predictAndEvaluationSen(model.best = model_use,
                                    test.data = data_test[index,],
                                    prevalence = prevalence,
                                    is.ensemble = FALSE,sen_min = sen_min,spe_min = spe_min)

    sen.vec <- c(sen.vec, as.numeric(result_tmp$test.prediction.evaluation[1]))
    spe.vec <- c(spe.vec, as.numeric(result_tmp$test.prediction.evaluation[2]))
    ppv.vec <- c(ppv.vec, as.numeric(result_tmp$test.prediction.evaluation[3]))
    npv.vec <- c(npv.vec, as.numeric(result_tmp$test.prediction.evaluation[4]))
    auc.vec <- c(auc.vec, as.numeric(result_tmp$test.prediction.evaluation[5]))

}

print("Finished!")

test_set_prediction_evaluation <- NULL
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(sen.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(spe.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(ppv.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(npv.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(auc.vec),4)))

rownames(test_set_prediction_evaluation) <- c("Sensitivity",
                                              "Specificity",
                                              "PPV",
                                              "NPV",
                                              "AUC")


return(test_set_prediction_evaluation)

}

### predict and evaluate functions
predictAndReturnProb <- function(model.best,
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
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
    ### ROC and AUC
    if (is.ensemble == FALSE){
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob[,2]
    } else{
        model.prob = test.prediction.prob
    }

    ## return the results
    return(list(positive.prob = model.prob))
}
