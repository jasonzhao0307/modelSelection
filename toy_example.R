# load all functions
source("model_selection.R")



########   Example   ########

# load example dataset
# we have 4 objects here:
# - df.training, with rows being samples, columns being variables.
# - df.testing, with rows being samples, columns being variables.
# - y.train, a binary vector
# - y.test, a binary vector
data.example <- readRDS("model_selection_example_data.RDS")


# prepare the data to make the format useable by the pipeline
# [Important!] User needs to specify which class in y is the positive class
data.prepared <- prepareData(df.training = data.example$df.training,
                             df.testing = data.example$df.testing,
                             y.train = data.example$y.train,
                             y.test = data.example$y.test,
                             positive = "cancer")

## Data normalization
data.norm <- normData(training = data.prepared$training,
                      testing = data.prepared$testing)


### Model selection

model.selection.output <- modelSelection(df.training = data.norm$trainTransformed,
                      model.names = c("gbm", "svmLinear", "glmnet"),
                      num.cv.fold = 5,
                      num.cv.repeat = 5,
                      num.param.tune = 6,
                      seed.use = 800,
                      show.some.models = TRUE)

## Prediction on testing set and evaluation
## [Note]: user needs to provide the prevalence for the positive case

# The performance of the single best model
predicttion.single.best <- predictAndEvaluation(model.best = model.selection.output$model.best,
                                    test.data = data.norm$testTransformed,
                                    prevalence = 0.1,
                                    is.ensemble = FALSE)

# The performance of the ensemble model
predicttion.ensemble <- predictAndEvaluation(model.best = model.selection.output$model.ensemble,
                                    test.data = data.norm$testTransformed,
                                    prevalence = 0.1,
                                    is.ensemble = TRUE)



