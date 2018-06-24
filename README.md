# modelSelection

##### Yue (Jason) Zhao 


modelSelection is a R pipeline designed for finding best binary classification model, predicting the testing data, and evaluating the prediction.


#### The pipeline contains the following files:

* model_selection.R
	* defines all functions
* model_selection_example_data.RDS
	* contains a toy dataset
* toy_eample.R
	* provides a toy example on how to use this pipeline


#### Inputs
* 4 data objects:
	* df.training, with rows being samples, columns being variables.
	* df.testing, with rows being samples, columns being variables.
	* y.train, a binary vector
	* y.test, a binary vector