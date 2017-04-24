# Simple test.
#
#  create a model and make sure it works.
test.ModelConstruction <- function() {
	data(fs1)
	model <- flowFPModel(fs1)
	checkTrue(!is.null(model))
	checkTrue(is(model) == "flowFPModel")
}

# More advanced test.
#
# create a model with a name, and three parameters, check that the model got 
# created correctly.
test.ModelAccessors <- function() {
	data(fs1)
	TEST_NAME = "bluepudding"
	model <- flowFPModel(fs1, nRecursions=3, parameters=c(1,3,7), name=TEST_NAME)
	checkTrue(!is.null(model))
	checkTrue(is(model) == "flowFPModel")
	checkTrue(nFeatures(model) == 8)
	checkTrue(name(model) == TEST_NAME)
	checkTrue(all(parameters(model) == c("FS Lin","FL1 Log","FL5 Log")))
	
}

test.CheckAutoLevel <- function() {
	data(fs1)
	TEST_NAME = "bluepudding"
	model <- flowFPModel(fs1, nRecursions="auto", parameters=c(1,5,6), name=TEST_NAME)
	checkTrue(!is.null(model))
	checkTrue(is(model)[1] == "flowFPModel")
	checkTrue(nRecursions(model) == 11)
}

test.NonFlowData <- function() {
	mat <- matrix(1:100, nrow=10, ncol=10)
	
	checkException(flowFPModel(mat), silent=TRUE)
}

test.BadParams <-function() {

	data(fs1)
	checkException(flowFPModel(fs1, parameters=c("FS Lin","FL1 Log","FL5 typeO")), 
	               silent=FALSE)

	checkException(flowFPModel(fs1, parameters=c(1,25)), 
	               silent=FALSE)
}
