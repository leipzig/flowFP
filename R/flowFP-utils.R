##
## Package: flowFP
## File: flowFP-utils.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

checkType <- function(object, appropriate_types, funcName) {

	type = paste(appropriate_types, collapse=", ")
	if (is.null(object)) {
		stop(paste(funcName, "() got passed a NULL object to check for type -->",
		           type, "<--\n", sep=""), call.=FALSE)
		return()
	}
		
	if (any(is(object) %in% appropriate_types)) # object is one of the appropriate_types specified.
		return()
		
	stop(paste ("Wrong type of object specified in call to ", funcName, "().\n",
		        "Only -->", type, "<-- are valid\n" , 
		         sep=""), call.=FALSE)
		   
}


parse_parameters <-function(names, params){

	if(is.null(params))
		return(names)
		
	if(is.numeric(params)) {
		if (length(which(is.na(names[params]))) > 0) {
    		stop("'parameters' is out of bounds for this flowSet or flowFrame.\n", call.=FALSE)
  		}
		return(names[params])
	}
		
	
	if(is.character(params)) {
		if (all(params %in% names)) {
			return(params)
		} else {
			stop("One or more 'parameters' don't match this flowSet or flowFrame.\n", call.=FALSE)
		}
	} else {
		stop("Invalid type for 'parameters', must either be numeric or character.\n",
		     call.=FALSE)
	}
}
