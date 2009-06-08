##
## Package: flowFP
## File: flowFPPlex.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

## =========================================================================
## flowFPPlex - constructor
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
flowFPPlex <- function(fingerprints=NULL) {

	plex = new("flowFPPlex")
	
	if (is.null(fingerprints))
		return (plex)

	plex = flowFPPlexAppend(plex, fingerprints)
	
	return(plex)
		
}

is.flowFPPlex <-function(obj) {
	return( is(obj)[1] == "flowFPPlex")
}

## function: flowFPPlexAppend
##
## This basic form was copied from the R base append function
## and was modified to use flowFPPlex objects.
##
flowFPPlexAppend <- function(x, values, after = length(x@fingerprints)) {

	if (is.flowFPPlex(values)) {
		fingerprints = values@fingerprints
	} else if (is.flowFP(values)){
		fingerprints = c(values)
	} else if (is.list(values)) {
		fingerprints = values
	} else {
		stop("You are trying to append the wrong type of object to a flowFPPlex\n", 
		     code.=FALSE)
	}

	validateFingerprintsForPlex(x, fingerprints)

	lengx <- length(x@fingerprints)
	
    if (after <= 0) {
        x@fingerprints = c(fingerprints, x@fingerprints)
    } else if (after >= lengx) {
        x@fingerprints = c(x@fingerprints, fingerprints)
    } else {
    	x@fingerprints = c(x@fingerprints[1:after], fingerprints, 
    	                   x@fingerprints[(after + 1):lengx])
    }
    return (x)
}

validateFingerprintsForPlex <- function(plex, fingerprints) {

	if( ! all(sapply(fingerprints, is.flowFP)) ) {
		wrong_type = which(!sapply(fingerprints, is.flowFP))
		errmsg = paste("You are trying to add/create a flowFPPlex out of object of type\n",
		                sapply(fingerprints[wrong_type], is),
		               " plex only accepts flowFP objects.\n")
			stop(errmsg, code.=FALSE)
	}
	
	if (length(plex) == 0)
		return (TRUE) # first fingerprint is always valid
		
	prototypeFP = plex[[1]]
	tmp = sampleNames(prototypeFP)
	tmpHasClasses = hasClasses(prototypeFP)
	tmpClasses = sampleClasses(prototypeFP)
	for(i in 1:length(fingerprints)) {
		if (!all(tmp == sampleNames(fingerprints[[i]]) ) ) {
			stop("The [", i,"] fingerprint in this set has a different",
			      " collection of sample names\n")
		}
		if (tmpHasClasses) {
			if (nlevels(tmpClasses) != nlevels(sampleClasses(fingerprints[[i]])) ||
			    !all(tmpClasses == sampleClasses(fingerprints[[i]])) ) {
				stop("The [", i,"] fingerprint in this set has a different",
			      " collection of sample classes\n")
			}
		}
	}
	return (TRUE)
}
