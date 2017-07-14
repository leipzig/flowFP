##
## Package: flowFP
## File: AllClasses.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##
##  This declarative file contains the class definitions for the three main
## classes used and exported by this package. It also contains a private class
## definition, used to hold bin boundaries.
##


## =========================================================================
## flowFPModel - Definition.
##
## parameters - a vector of either parameter names e.g. c("FSC-H", "SSC-H")
##              or the parameter number e.g c(1,5,7) used to construct the model
##              no other parameters from the file will be used.
##
## nRecursions- Number of levels of sub-division. The number of leaf bins is
##              equal to 2 ^ this value. nRecursions = 5 means 2^5 or 32 leaf
##              bins. The Constructor check to see if nRecursions is too big
##              based on the total number of events present in the training
##              set.
##
## trainingSet - List of sample names from the flowSet or flowFrame used to create the
##               model.
## trainingSetParams - All of the Parameters from the training flowSet or flowFrame.
##
## dequantize - If TRUE the FCS data in each flowFrame has a small ramp added
##              to each event. for each event in each parameter we add alpha * i.
##              where i goes from 1 to number of events in the particular 
##              flowFrame. We do this to eliminate duplicate values in the FCS
##              data section. Since we use the middle value to sub-divide, this
##              make it possible to keep the bins almost perfectly balanced, and
##              the results are reproducible.
##
## split_val & split_axis - These lists, are kept together and describe the sub-
##                          division process. They have 'nRecursions' items and
##                          each item is a vector that represents the splits for
##                          that level. So the first level has exactly one val
##                          and axis. The lengths go up as powers of 2.
##                          When the model is used by the fingerprint constructor,
##                          it has to split the events recursively, to arrive
##                          at the correct result.
##
## binBoundary - This slot contains lower left, and upper right coordinates for
##               each bin boundary in the model.
##
## .cRecursions - This private slot holds the value of the number of recursions
##                used to construct the model, this value should not be set directly
##                by the user. It defines the maximum resolution for the model.
##
## .tmp_tags - Are just used for bookkeeping, and are set to empty after the model
##            has been constructed.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("flowFPModel", representation(name="character", parameters="vector", 
                                       nRecursions="numeric",
                                       trainingSet= "vector", 
                                       trainingSetParams="vector",
                                       dequantize="logical", split_val="list",
                                       split_axis="list", binBoundary="list",
                                       .cRecursions="numeric",
                                       .tmp_tags="vector"),
         prototype = list(name="", parameters=vector(length=0), nRecursions=0,
                          trainingSet=vector(mode="character"),
                          trainingSetParams=vector(mode="character"),
                          dequantize=FALSE, split_val=list(), split_axis=list(),
                          binBoundary=list(),
                          .cRecursions=0,
                          .tmp_tags=vector(mode="integer", length=0))
)

setClass("binBoundary", representation(ll="vector", ur="vector"),
         prototype = list(ll=vector(length=0), ur=vector(length=0))
)

## =========================================================================
## flowFP - Definition.
##
## counts - an integer matrix, organized one instances per row, cols are the
##          number of events that populate that bin in the model.
##
## tags - A hairy array, (each row can be of different length), The rows are
##        instances, and the length is the number of events in the flowFrame.
##        The values are the integer bin number.
##
## sampleNames - Character vector of names for each instance.
##
## sampleClasses - factor, used to hold clinical class for an instance. e.g. 
##                  (Sick, Well), or (good, evil) ... not an OO class name.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("flowFP",
         representation(counts = "matrix",
                        tags = "list",
                        medians = "list",
                        sampleNames = "vector",
                        sampleClasses = "factor"), 
         contains="flowFPModel",
         prototype = list(counts = matrix(nrow = 0, ncol = 0),
                          tags = list(),
                          medians = list(),
                          sampleNames = vector(mode="character"),
                          sampleClasses = factor()
                        ),
         validity=function(object){
         	return (TRUE)
         }
)

## =========================================================================
## flowFPPlex - Definition.
##
## fingerprints - a list of flowFP's
##
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("flowFPPlex",
         representation(fingerprints = "list"), 
         prototype = list(fingerprints = list()),
         validity=function(object){
         	return (TRUE)
         }
)
