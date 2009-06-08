##
## Package: flowFP
## File: plot-model.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

plotModel <-function(x, parameters=NULL, alpha=1, border="gray", 
                     showbins=1:nFeatures(x), ylim=NULL, xlim=NULL, main="Model", ...){

	parameters = getDefaultParameters(x, parameters)
	check_plot_args(x, fs=NULL, highlight=NULL, showbins=showbins)
	
  	xlab = x@trainingSetParams[parameters[1]]
  	ylab = x@trainingSetParams[parameters[2]]
  	
  	if (is.null(ylim) || is.null(xlim))
  		limits = getBoundaryLimit(x, showbins, parameters)
  	
  	
  	if(is.null(ylim))
  		ylim = limits$ylim
  		
  	if (is.null(xlim))
  		xlim = limits$xlim
  		
	plot('', type='l', ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, main=main)
	plotBinBoundaries(x, y=NULL, parameters, alpha, border, showbins, ...)
	
}

plotModelWithDots<-function(x, y, parameters=NULL, alpha=.5, border="gray", 
                            showbins=1:nFeatures(x), main="Model", ...){

	parameters = getDefaultParameters(x, parameters)
	
	if (is(y) == "flowSet") {
		y = as(y, "flowFrame")
	}
  
  	xlab = x@trainingSetParams[parameters[1]]
  	ylab = x@trainingSetParams[parameters[2]]
  	plotParamNames = x@trainingSetParams[parameters]
  	length(plotParamNames) = 2
  	ffplot = getMethod("plot", signature("flowFrame", "character"))
	ffplot(y, plotParamNames, smooth=FALSE, main=main, ...)
	plotBinBoundaries(x, y=y, parameters, alpha, border, showbins, ...)
}


plotBinBoundaries <-function(x, y=NULL, parameters=NULL, alpha=1, border="gray", 
                     showbins=1:nFeatures(x), pch='.', cex=1, ...) {

	if (!is.flowFPModel(x))
		all_tags = unlist(tags(x))
	parameters = getDefaultParameters(x, parameters)
	len = length(showbins)
	
	binBoundary = binBoundary(x)
	
	xleft = vector(mode="numeric", length=len)
	ybottom = vector(mode="numeric", length=len)
	xright = vector(mode="numeric", length=len)
	ytop =vector(mode="numeric", length=len)
	
	bin_cols = binColorLut(nFeatures(x), showbins, alpha=alpha)
	point_cols = binColorLut(nFeatures(x), showbins, alpha=1)

	if (!is.null(y))
		y = as(y, "flowFrame")
		
	j = 1
	for(i in showbins) {
		xleft[j] = binBoundary[[i]]@ll[parameters[1]]
		ybottom[j] = binBoundary[[i]]@ll[parameters[2]]
		xright[j] = binBoundary[[i]]@ur[parameters[1]]
		ytop[j] = binBoundary[[i]]@ur[parameters[2]]
		if (!is.flowFPModel(x)) {
			event_idx = which(all_tags %in% i)
            points(exprs(y)[event_idx, parameters[1]], 
                   exprs(y)[event_idx, parameters[2]], pch=pch, 
                   col=point_cols[j], cex=cex)
		}
		j = j + 1
	}
  
	rect(xleft, ybottom, xright, ytop, density=NA, 
		col=bin_cols, border=border, ...)
}

getBoundaryLimit <- function(x, showbins, parameters) {

	len = length(showbins)
	xleft = vector(mode="numeric", length=len)
	ybottom = vector(mode="numeric", length=len)
	xright = vector(mode="numeric", length=len)
	ytop =vector(mode="numeric", length=len)
	j = 1
	for(i in showbins) {
		xleft[j] = x@binBoundary[[i]]@ll[parameters[1]]
		ybottom[j] = x@binBoundary[[i]]@ll[parameters[2]]
		xright[j] = x@binBoundary[[i]]@ur[parameters[1]]
		ytop[j] = x@binBoundary[[i]]@ur[parameters[2]]
		j = j + 1
	}
	xlim = range(c(xleft, xright),  na.rm=TRUE, finite=TRUE)
	ylim = range(c(ybottom, ytop),  na.rm=TRUE, finite=TRUE)
	return (list(xlim=xlim, ylim=ylim))
}
