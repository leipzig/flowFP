##
## Package: flowFP
## File: flowFPModel.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

## =========================================================================
## Constructor for flowFPModel.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
flowFPModel <- function(fcs, name="Default Model", parameters=NULL, 
                        nRecursions='auto', dequantize=TRUE, sampleSize=NULL,
                        excludeTime=TRUE)
{

	checkType(fcs, c("flowSet", "flowFrame"), "flowFPModel")

	if (dequantize)
		fcs = dequantize(fcs)

	if(is(fcs,"flowSet"))
	{
		trainingSet = sampleNames(fcs)
		if (is.null(sampleSize)) {
			training_events = median(fsApply(fcs, nrow)) # for auto nRecursions
			fcs <- as(fcs, "flowFrame")
			exprs(fcs) <- subset(exprs(fcs), subset=T, select=colnames(fcs)[1:(ncol(exprs(fcs)) -1)])
		}else {
			nevents = vector(mode="numeric", length=length(fcs))
			nevents[1] = min(sampleSize, nrow(fcs[[1]]))
			tmp = fcs[[1]]
			ind = sample(1:nrow(fcs[[1]]), nevents[1])
			exprs(tmp) = exprs(tmp)[ind,]
			if (length(fcs) > 1) {
				for(i in 2:length(fcs)) {
					nevents[i] = min(sampleSize, nrow(fcs[[i]]))
					ind = sample(1:nrow(fcs[[i]]), nevents[i])
					exprs(tmp) = rbind(exprs(tmp), exprs(fcs[[i]])[ind,])
				}
			}
			training_events = median(nevents)
			fcs = tmp
		}
	} else {
		if (length(fcs) == 0) {
			stop("Can not create a model from a flowFrame with zero events\n")
		}
		trainingSet = identifier(fcs)
		if(is.null(sampleSize)) {
			training_events = nrow(fcs)
		} else {
			training_events = min(sampleSize, nrow(fcs))
			ind = sample(1:nrow(fcs), training_events)
			exprs(fcs) = exprs(fcs)[ind,]
		}
	} 

	parameters = parse_parameters(colnames(fcs), parameters, excludeTime)
	
	if (nRecursions == 'auto') {
		nRecursions = (log2(training_events) %/% 1) - 3 # at least 8 events / bin
	}

	validate_params(fcs, parameters, nRecursions)
  

	model = new("flowFPModel", name=name, parameters=parameters, 
	            nRecursions=nRecursions, .cRecursions=nRecursions,
	            trainingSet=trainingSet, dequantize=dequantize, 
	            trainingSetParams=colnames(fcs))

	model@.tmp_tags = vector(mode = "integer", length = nrow(fcs))
	model@.tmp_tags[] = as.integer(1)
	for(i in 1:model@nRecursions) {
		model <- bin_level(fcs, model, i)
	}

	model@binBoundary = createBinBoundaries(model, fcs)
	# empty temp tags & hope the GC returns memory.
	model@.tmp_tags = vector(mode = "integer", length = 0)
	
  return (model)
}



## =========================================================================
## Helper Functions for flowFPModel.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


is.flowFPModel <-function(obj) {
	return( is(obj)[1] == "flowFPModel")
}


##
## Function: bin_level
##
## This private function sub-divides the multi-dimensional fcs data into
## some number of bins. It gets called 'nRecursions' times. It keeps track
## of which event belongs in which bin, by way of the model@.tmp_tags.
## There is a one to one relationship between the .tmp_tags int vector
## and the events (rows) in the fcs data matrix. As we move through the
## process, each bin gets divided into 2 bins, thus the number of bins
## equal 2^nRecursions.
bin_level <- function(fcs, model, level) {
  # number of bins at this level
  num_splits = as.integer(2**(level -1))
  model@split_val[[level]] = vector(mode="numeric", length=num_splits)
  model@split_axis[[level]] = vector(mode="integer", length=num_splits)

  param_idx = parameters(model, index=TRUE)

  .Call("bin_level", fcs@exprs, model@.tmp_tags, model@split_axis[[level]], 
        model@split_val[[level]], level, param_idx)
  return (model)
}

##
## Function: createBinBoundaries
##
## This function is a book-keeping nightmare... It figures out in multi-dimensional
## space the bin boundaries of each bin in a model.
createBinBoundaries <-function(model, fcs) {

	nRecursions = model@.cRecursions
	numfeatures = 2^nRecursions
	lenCoord = ncol(fcs)

	## Find the boundary of all of the training data used to create this model.
	## This is the total extent of the model.
	rng = range(fcs)
	ll = unlist(rng["min",])
	ur = unlist(rng["max",])
	
	## create the starting bin 'level 0' contains all of the data.
	tmpBB = new("binBoundary", ll=ll, ur=ur)
	binBList = list();
	binBList[[1]] = tmpBB
	for(i in 1:nRecursions) {
		inplay = seq(1, numfeatures, by=2^(nRecursions - i))
		odd = seq(1, length(inplay), by=2)
		even = seq(2, length(inplay), by=2)
		src = inplay[odd]
		dest = inplay[even]
		for(j in 1:length(src)) {
			binBList[[dest[j]]] = binBList[[src[j]]]
			split_axis = model@split_axis[[i]][j]
			split_val = model@split_val[[i]][j]
			binBList[[dest[j]]]@ll[split_axis] = split_val
			binBList[[src[j]]]@ur[split_axis] = split_val
		}
		
	}

	return (binBList)
}

##
## Function: dequantize
##
## This function adds a small ramp to all of the events in an fcs flowFrame
## to insure that each value is unique. We invented this approach to make
## each event in an fcs flowFrame unique, based on position. 8 bit machines
## for example produce only 256 unique values, if you collect 100,000 events
## you would expect to find lots of repeated values, if you continue to divide
## the space up, you will end up, creating bins with no events in them, this
## fixes that in a reproducible way.
dequantize <- function(x, alpha=1e-8) {

	## Note: This copies the flowFrame or flowSet that got passed in
	##       so we leave the original data as it was.
	fcs <- x
  	if (is(fcs, "flowFrame")) {
  		if (nrow(fcs) > 0) {
	  	  ofs <- seq(alpha, alpha * nrow(fcs), length.out=nrow(fcs))
		  exprs(fcs) <- exprs(fcs) + ofs
		}
	} else {
		for(i in 1:length(fcs)) {
			if (nrow(fcs[[i]]) > 0) {
				ofs <- seq(alpha, alpha * nrow(fcs[[i]]), length.out=nrow(fcs[[i]]))
			    exprs(fcs[[i]]) <- exprs(fcs[[i]]) + ofs
			}
		}
		
	}
  return(fcs)
}


## This private function checks the input parameters for validity. Some
## of the rules are:
##  - The FCS parameters list cannot contain NA's. 
##  - The number of leaf bins, (based on nRecursions) must be smaller
##    than the number of events in the fcs flowFrame.
##
validate_params <-function(x, parameters, nRecursions) {
  if (length(which(is.na(parameters))) > 0) {
    stop("ERROR: 'parameters' is out of bounds for this FCS data\n")
  }

  if (2^nRecursions > nrow(x)) {
    stop("ERROR: Model with too many nRecursions, max num nRecursions = ", 
         as.integer(log2(nrow(x)) - 1), "\n")
  }
  return
}

