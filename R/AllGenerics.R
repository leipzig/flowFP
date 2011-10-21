##
## Package: flowFP
## File: AllGenerics.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##
##  This file contains the generic methods for the classes used in this package.
##  They are grouped together by name, and ordered in each group by lowest level
##  class and smallest signature to more complicated. My driving principle has
##  been that I should make a separate function, if the code is more than a 
##  couple of lines, so I can test it outside of the package environment. 

## =========================================================================
## append Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("append")) {
	setGeneric("append", function(x, values, after=length(x)) standardGeneric("append"))
}

setMethod("append", signature=signature(x="flowFPPlex", values="flowFP"),
			definition=flowFPPlexAppend)

setMethod("append", signature=signature(x="flowFPPlex", values="flowFPPlex"),
			definition=flowFPPlexAppend)

setMethod("append", signature=signature(x="flowFPPlex", values="list"),
			definition=flowFPPlexAppend)
          
## =========================================================================
## counts Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("counts")) {
	setGeneric("counts", function(object, ...) standardGeneric("counts"))
}

setMethod("counts",
	signature=signature(object="flowFP"),
	definition=getFPCounts)

setMethod("counts",
	signature=signature(object="flowFPPlex"),
	definition=function(object, ...) {
		cnts = matrix(0, nrow=nInstances(object), ncol=0)
		for(fp in object@fingerprints)
			cnts = cbind(cnts, counts(fp, ...))
		return(cnts)
	}
)

## =========================================================================
## Indexing Methods. [ [[
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("[",
	signature=signature(x="flowFP", i="numeric"),
	definition=function(x, i, j="missing", ..., drop=FALSE)  {
      	if(missing(i))
              return(x)
              
    	copy <- x
    	copy@counts <- matrix(x@counts[i,], nrow=length(i), ncol=ncol(x@counts))
        copy@tags <- x@tags[i]
        copy@sampleNames <- x@sampleNames[i]
    	if (length(x@sampleClasses) > 1) {
        	copy@sampleClasses <- x@sampleClasses[i]
        }
        return (copy)

      })

setMethod("[",
	signature=signature(x="flowFP", i="character"),
	definition=function(x, i, j="missing", ..., drop=FALSE)  {
      	if(missing(i))
              return(x)
              
    	if (length(x@sampleClasses) > 0 && !anyMissing(match(i, levels(x@sampleClasses)))) {
    		i = as.vector(sapply(i, function(obj) { which(as.character(x@sampleClasses) == obj) }))
    	} else if (!anyMissing(match(i, x@sampleNames))) {
    		i = match(i, x@sampleNames)
    	} else {
    		stop(paste("Can't use ", i, " to subscript this fingerprint.\nSubscript not found\n") ,
    		     call.=FALSE)
    	}

    	copy <- x
    	copy@counts <- matrix(x@counts[i,], nrow=length(i), ncol=ncol(x@counts))
        copy@tags <- x@tags[i]
        copy@sampleNames <- x@sampleNames[i]
    	if (length(x@sampleClasses) > 1) {
            copy@sampleClasses <- x@sampleClasses[i]
        }
        return (copy)

      })

setMethod("[",
          signature=signature(x="flowFPPlex"),
          definition=function(x, i, j="missing", ..., drop=FALSE)
      {
      	if(missing(i))
              return(x)
       new("flowFPPlex", fingerprints=x@fingerprints[i])
	  }
)

setMethod("[[",
          signature=signature(x="flowFPPlex"),
          definition=function(x, i, j="missing", ..., drop=FALSE)
      {
      	if(missing(i)) {
              stop("Accessing the flowFP out of a plex needs a single value\n", 
                   call.=FALSE)
        }
              
        if(!is.numeric(i)) {
        	stop("you must specify a single value for the [[ index\n", call.=FALSE)
        }
        
       return(x@fingerprints[[i]])
	  }
)

setReplaceMethod("[[",
	signature=signature(x="flowFPPlex", i="numeric", value="flowFP"),
	definition=function(x, i, j="missing", ..., value)
	{
		if(missing(i)) {
			stop("Accessing the flowFP out of a plex needs a single value\n",
			     call.=FALSE)
		}
  
		if(!is.numeric(i)) {
			stop("you must specify a single value for the [[ index\n", 
			     call.=FALSE)
		}
		x@fingerprints[[i]] = value
		return(x)
	  }
)
	
## =========================================================================
## name Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("name")) {
	setGeneric("name", function(object) standardGeneric("name"))
}
setMethod("name",
	signature=signature(object="flowFPModel"),
	definition=function(object) {
		object@name
	}
)

setMethod("name",
	signature=signature(object="flowFP"),
	definition=function(object) {
		object@name
	}
)

setMethod("name", 
	signature=signature(object="flowFPPlex"),
	definition=function(object) {
		sapply(object@fingerprints, name)
	}
)

## =========================================================================
## These methods are used to determine the size of our objects.
## nInstances - returns the number of flowFrames used to create a
##              fingerprint, and a plex can only be constructed with
##              a group of fingerprints share the same instances.
##
## nFeatures - returns the number of features, that will be generated from
##             a model, or that are in an individual fingerprint or the
##             total number of features in a plex.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (!isGeneric("nInstances")) {
	setGeneric("nInstances", function(object) standardGeneric("nInstances"))
}

setMethod("nInstances", 
	signature=signature("flowFP"), 
	definition=function(object) { length(object@sampleNames) } )

setMethod("nInstances",
	signature=signature("flowFPPlex"), 
	definition=function(object) { 
		if (length(object@fingerprints) == 0)
			return (0)
		nInstances(object@fingerprints[[1]]) 
	}
)

if (!isGeneric("nFeatures")) {
	setGeneric("nFeatures", function(object) standardGeneric("nFeatures"))
}

setMethod("nFeatures", 
	signature=signature("flowFPModel"), 
	definition=function(object) { 
		2^object@nRecursions
	}
)

setMethod("nFeatures",
	signature=signature("flowFP"), 
	definition=function(object) { 
		2^object@nRecursions
	}
)

setMethod("nFeatures",
	signature=signature("flowFPPlex"), 
	definition=function(object) { 
		sum(sapply(object@fingerprints, nFeatures))
	}
)

if (!isGeneric("length")) {
	setGeneric("length", function(x) standardGeneric("length"))
}

setMethod("length",
	signature=signature(x="flowFPPlex"), 
	definition=function(x) { 
		return (length(x@fingerprints))
	}
)

setMethod("length<-",
	signature=signature(x="flowFPPlex", value="numeric"), 
	definition=function(x, value) { 
		length(x@fingerprints) = value
		return(x)
	}
)



## =========================================================================
## nRecursions Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("nRecursions")) {
	setGeneric("nRecursions", function(object) standardGeneric("nRecursions"))
}

setMethod("nRecursions", 
	signature=signature("flowFPModel"), 
	definition=function(object) { 
		object@nRecursions
	}
)

setMethod("nRecursions",
	signature=signature("flowFP"), 
	definition=function(object) { 
		object@nRecursions
	}
)

if (!isGeneric("nRecursions<-")) {
	setGeneric("nRecursions<-", function(object, value) standardGeneric("nRecursions<-"))
}

setMethod("nRecursions<-",
	signature=signature("flowFP", value="numeric"), 
	definition=function(object, value) { 
		if (value > object@.cRecursions) {
			stop(paste("The underlying model for this fingerprint only supports",
			            object@.cRecursions, "levels of recursion.\n",
			            " To view more levels of recursion recreate the model/fingerprint with a larger nRecursions value\n"),
			            call.=FALSE)
		}
		object@nRecursions=value
		return(object)
	}
)

setMethod("nRecursions<-",
	signature=signature("flowFPModel", value="numeric"), 
	definition=function(object, value) { 
		if (value > object@.cRecursions) {
			stop(paste("This model was computed with", object@.cRecursions, 
			           "levels of recursion and cannot be increased above this value.\n"),
			            call.=FALSE)
		}
		object@nRecursions=value
		return(object)
	}
)

setMethod("nRecursions",
	signature=signature("flowFPPlex"), 
	definition=function(object) { 
		sapply(object@fingerprints, nRecursions)
	}
)

## =========================================================================
## parameters Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.getparameters <- function(object, used=FALSE, index=FALSE) {
	if (used) {
       	used = sort(unique(unlist(object@split_axis)))
    	params = object@trainingSetParams[used]
    } else {
    	params = object@parameters
    }
    if (index)
    	return (which(object@trainingSetParams %in% object@parameters))
    else
    	return (params)
}

setMethod("parameters", signature=signature("flowFP"),
           definition= .getparameters)

setMethod("parameters", signature=signature("flowFPModel"),
           definition= .getparameters)
                        


## =========================================================================
## Plot Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("plot", 
	signature=signature(x="flowFPModel", y="missing"),
	definition=plotModel)

setMethod("plot", 
	signature=signature(x="flowFPModel", y="flowFrame"),
	definition=plotModelWithDots)
	
setMethod("plot", 
	signature=signature(x="flowFPModel", y="flowSet"),
	definition=plotModelWithDots)

setMethod("plot", 
	signature=signature(x="flowFP", y="missing"), 
	definition=function(x, y, type=c("tangle", "stack", "grid","qc", "plate"), ...) {
		type = match.arg(type)
		switch(type,
			tangle = plotTangleFP(x, ...),
			stack = plotStackFP(x, ...),
			grid = plotGridFP(x, ...),
			qc = plotqcFP(x, ...),
			plate = plotPlateFP(x,...),
			stop(paste("This type of plot =>", type, "<= is not implemented yet\n"), call.=FALSE)
		)
	}
)

setMethod("plot",  
	signature=signature(x="flowFP",  y="flowSet"), 
	definition=function(x, y, showbins=NULL, showfp=TRUE, ...) {
		plotFancy(x, y, showbins=showbins, showfp=showfp, ...)
	}
)

setMethod("plot",
    signature=signature(x="flowFP",  y="flowFrame"),
    definition=function(x, y, ...) {
    	if (nInstances(x) != 1) {
    		stop(paste("If you want to plot a 'flowFrame', you must specify exactly one fingerprint\n",
		     "For example you can subscript you fingerprint 'fp[1]'\n"), call.=FALSE)
		}
		fs = as(y, "flowSet")
		sampleNames(fs) = x@sampleNames
		plotFancy(x, fs, ...)
	}
)

setMethod("plot", 
	signature=signature(x="flowFPPlex", y="missing"), 
	definition=function(x, y, type=c("tangle", "stack", "grid","qc", "plate"), ...) {
		type = match.arg(type)
		switch(type,
			tangle = plotTangleFP(x, ...),
			grid = plotGridFP(x, ...),
			qc = plotqcFP(x, ...),
			stop(paste("This type of plot =>", type, "<= is not implemented yet\n"), call.=FALSE )
		)
	}
)

## =========================================================================
## sampleNames Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("sampleNames",   signature=signature(object="flowFP"),
	definition=function(object) {
		object@sampleNames
	}
)

setMethod("sampleNames",   signature=signature(object="flowFPPlex"),
	definition=function(object) {
		sampleNames(object@fingerprints[[1]])
	}
)

## =========================================================================
## sampleClasses Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("sampleClasses")) {
	setGeneric("sampleClasses", function(object) standardGeneric("sampleClasses"))
}
setMethod("sampleClasses",
	signature=signature(object="flowFP"),
	definition=function(object) {
		object@sampleClasses
	}
)

setMethod("sampleClasses",
	signature=signature(object="flowFPPlex"),
	definition=function(object) {
		sampleClasses(object@fingerprints[[1]])
	}
)

if (!isGeneric("hasClasses")) {
	setGeneric("hasClasses", function(object) standardGeneric("hasClasses"))
}
setMethod("hasClasses",
	signature=signature(object="flowFP"),
	definition=function(object) {
		return (length(sampleClasses(object)) != 0)
	}
)

setMethod("hasClasses",
	signature=signature(object="flowFPPlex"),
	definition=function(object) {
		return (length(sampleClasses(object)) != 0)
	}
)

if (!isGeneric("sampleClasses<-")) {
	setGeneric("sampleClasses<-", function(object, value) standardGeneric("sampleClasses<-"))
}

setMethod("sampleClasses<-",
	signature=signature(object="flowFP"),
	definition=function(object, value) {
		if (!is.factor(value)) {
			value = factor(value)
		}
		if (length(value) == 0 || (length(value) == length(object@sampleNames))) {
			object@sampleClasses = value
		} else {
			stop("The number of sampleClasses must equal the number of fingerprints." ,
    		     call.=FALSE)
		}
		return(object)
	}
)

setMethod("sampleClasses<-",
	signature=signature(object="flowFPPlex"),
	definition=function(object, value) {
		object@fingerprints = lapply(object@fingerprints, "sampleClasses<-", value)
		return (object)
	}
)

## =========================================================================
## Show Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
	signature=signature(object="flowFPModel"),
	definition=function(object)  {
		cat("A flowFPModel:\n\n")
		cat("  Name =", object@name, "\n")
		cat(sprintf("  nRecursions (max) = %d (%d)\n", object@nRecursions, object@.cRecursions))
		cat("  Dequantize =", object@dequantize, "\n")
		cat("  Parameters Considered:\n   ", paste(object@parameters,collapse=", "), "\n")
		tmp = unlist(object@split_axis)
		used = which(object@parameters %in% object@trainingSetParams[tmp])
		cat("  Parameters Used:\n   ", paste(object@parameters[used],collapse=", "), "\n")
		cat("  Training Set:\n   ", paste(object@trainingSet, collapse="\n    "), "\n")
	})

setMethod("show", signature=signature(object="flowFP"),
	definition=function(object) {
		cat("A flowFP containing", nInstances(object), "instances with",
	        nFeatures(object), "features.\n\n")
	    if (length(object@sampleClasses) > 0) {
	    	cat(sprintf("  %20s %20s\n", "Name", "Class"))
		    for(i in 1:nInstances(object)) {
    	    	cat(sprintf("  %20s %20s\n", object@sampleNames[i], object@sampleClasses[i]))
		    }
		} else {
			show(object@sampleNames)
		}
		
		cat("\n")
		cat("Extends ")
		show(as(object, "flowFPModel"))
	})

setMethod("show", signature=signature(object="flowFPPlex"),
	definition=function(object) {
		cat("A flowFPPlex containing ", length(object@fingerprints), "flowFPs\n")
		tmp = name(object)
		show(tmp)
	})


## =========================================================================
## Summary Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary", 
	signature=signature(object="flowFP"),
	definition=function(object) {
		cat("A flowFP containing", nInstances(object), "instances with",
		    nFeatures(object), "features each.\n\n")
		s = matrix(ncol=3, nrow=nInstances(object))
		colnames(s) = c("min", "max", "mean")
		rownames(s) = object@sampleNames
		for(i in 1 : nInstances(object)) {
			counts = object@counts[i,]
			mincnt = min(counts)
			maxcnt = max(counts)
			meancnt = mean(counts)
			s[i,] = c(mincnt, maxcnt, meancnt)
		}
		show(s)
		cat("\n")
	})



## =========================================================================
## tags Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("tags")) {
	setGeneric("tags", function(object) standardGeneric("tags"))
}


setMethod("tags", 
	signature=signature(object="flowFP"),
	definition=function(object) {
		if (object@nRecursions == object@.cRecursions)
			return(object@tags)
			
		f = 2^(object@.cRecursions - object@nRecursions) 
		return( lapply(object@tags, f=f, function(obj, f=1) { (obj - 1) %/% f + 1} ))
	}
)


## =========================================================================
## binBoundry Methods.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isGeneric("binBoundary")) {
	setGeneric("binBoundary", function(object) standardGeneric("binBoundary"))
}

setMethod("binBoundary", 
	signature=signature(object="flowFPModel"),
	definition=function(object) {
		return(object@binBoundary)
	}
)

setMethod("binBoundary", 
	signature=signature(object="flowFP"),
	definition=function(object) {
		if (object@nRecursions == object@.cRecursions)
			return(object@binBoundary)
		
		binBoundaries = list()
		end = 0
		step = 2^(object@.cRecursions - object@nRecursions) 
		for(i in 1:nFeatures(object)) {
			begin = end + 1
			end = begin + step - 1
			binBoundaries[[i]] = new ("binBoundary", ll=object@binBoundary[[begin]]@ll, 
			                        ur=object@binBoundary[[end]]@ur)
		}
		return (binBoundaries)
	}
)


