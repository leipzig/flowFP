##
## Package: flowFP
## File: flowFP.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

## =========================================================================
## flowFP - constructor
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
flowFP <- function(fcs, model=NULL, sampleClasses=NULL, ...) {

	checkType(fcs, c("flowSet", "flowFrame"), "flowFP")

    if (is.null(model)) {
        model = flowFPModel(fcs, ...)
    }


    fp = as(model, "flowFP")
    
    fp@tags = list();
    numFeatures = 2^fp@.cRecursions

    if(is(fcs,"flowSet")) {
        fp@counts = matrix(as.integer(0), nrow=length(fcs), ncol=numFeatures)
        fp@sampleNames = sampleNames(fcs)
        
        benchmarks<-list()
        for(i in 1:length(fcs)) {
            if(model@dequantize) {
                fcs[[i]] = dequantize(fcs[[i]])
            }
          
          #list by samples, each with all events for that sample with indices indicating which bin they belong to
        	fp@tags[[i]] = tag_events(fcs[[i]]@exprs, model)
        	#matrix of samples by bins, with number of events in that bin
        	fp@counts[i,] = count_events(fp@tags[[i]], numFeatures)
        	
        	
        	#median_events(fcs[[i]]@exprs,fp@tags[[i]])
        }
        if (!is.null(sampleClasses)) {
			sampleClasses(fp) <- sampleClasses
		}	
	} else {
		if(model@dequantize) {
			fcs = dequantize(fcs)
		}

		fp@sampleNames = identifier(fcs)

		fp@counts = matrix(as.integer(0), nrow=1, ncol=numFeatures)
		fp@tags[[1]] = tag_events(fcs@exprs, model)
		fp@counts[1,] = count_events(fp@tags[[1]], numFeatures)
	}

    return (fp)
}

is.flowFP <-function(obj) {
	return( is(obj)[1] == "flowFP")
}
## =========================================================================
## flowFP - Helper Functions...
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tag_events <- function(fcs_events, model) 
{
    tags = vector(mode="integer", length=nrow(fcs_events))

    if (nrow(fcs_events) == 0)
        return (tags)

    tags[] = as.integer(1)
    for(i in 1:model@nRecursions) {
        .Call("tag_events", fcs_events, i, model@split_axis[[i]], 
              model@split_val[[i]], tags)
    }
    return(tags)
}

## This private function acts as an iterface to the C library
## it allocates memory to hold the results of the counting.
##
count_events <- function(tags, numFeatures) {
    counts = vector(mode="integer", length=numFeatures)
    .Call("count_events", counts, tags)
    return (counts)
}

## This private function should only be called by the method 'counts',
## with the signature "flowFP". It takes care of combinding counts to produce
## a lower resolution version, 'view' of a fingerprint.
getFPCounts <- function(object, transformation=c("raw", "normalized", "log2norm")) {

	transformation = match.arg(transformation)

	if (object@nRecursions == object@.cRecursions) {
		counts = object@counts
	} else {
	    ## Make new counts. This needs to be cleaned up. HAH
		step = 2^(object@.cRecursions - object@nRecursions)
		counts = matrix(0, nrow=nrow(object@counts), ncol=ncol(object@counts)/step)

		end = 0 # initialize
		for(i in 1:ncol(counts)) {
			begin = end + 1
			end = begin + step - 1
			counts[,i] = rowSums(object@counts[,begin:end])
		}
	}
	if (transformation == "raw")
		return (counts)
		
	expected_counts = sapply(object@tags, length) / nFeatures(object)
	counts = counts / expected_counts
	
	if (transformation == "normalized")
		return (counts)
	else
		return (log2(counts))
}


# fcs_events is a matrix for that sample
# tags is the bin index for each event
# a flowSet will report tags as a list, this can be safely unlisted regardless
getFPbinCentroids<-function(fp,fcs){
  if(is(fcs,"flowSet")) {
    fcs <- as(fcs, "flowFrame")
  }
  
  
  if (nrow(fcs) != length(unlist(fp@tags))) {
    stop(paste("ERROR: the number of events in the flowFrame or flowSet (",nrow(fcs),") does not match the number of events in the flowFP model (",length(unlist(fp@tags)),") for this FCS data\n"))
  }

  
  binCentroids<-do.call(rbind,lapply(split(seq_len(nrow(fcs@exprs)), as.factor(unlist(fp@tags)), drop = FALSE), function(ind) apply(fcs@exprs[ind, , drop = FALSE], 2, median)))

  return(binCentroids)
}
