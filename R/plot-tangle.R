##
## Package: flowFP
## File: plot-tangle.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

plotTangleFP <-function(x, transformation=c("raw", "normalized", "log2norm"), 
                        linecols=NULL, highlight=NULL, ylim=NULL, 
                        useClasses=FALSE, main="Fingerprints", xlab='Feature Index', ...){

	transformation = match.arg(transformation)
	
	check_plot_args(x, fs=NULL, highlight=highlight, showbins=NULL)
	
  	initFPplot(x, ylim=ylim, main=main, 
  	          xlab=xlab,
  	          ylab=paste(transformation, "Fingerprint", sep=' '), 
              transformation=transformation, ...)
                       
  	draw_tangle_lines(x, transformation, linecols, highlight, useClasses=useClasses, xlab=NULL, ...)
}


highlightColors <- function(colors, highlights) {
	for(i in 1:length(colors)) {
		tmp = col2rgb(colors[i])
		tmp = rgb2hsv(r=tmp[1], g = tmp[2], b = tmp[3])
		if (i %in% highlights) {
			colors[i] = hsv(h=tmp[1], s=tmp[2], v=tmp[3], alpha = 1)
		} else {
			colors[i] = hsv(h=tmp[1], s=tmp[2], v=tmp[3], alpha = .25)
		}
	}
	return(colors)
}

draw_tangle_lines<- function(object, transformation, linecols, highlight, 
                             useClasses=FALSE, ...){

	data = counts(object, transformation=transformation)
    	
	data[which(is.infinite(data))] <- NA
	ind <- 1:nFeatures(object)
	nfp <- nInstances(object)

	if (useClasses) {
		if (!hasClasses(object)) {
			stop(paste("This set of fingerprints does not contain class information\n",
			     "Use sampleClasses method to set the classes as a 'factor' for the fingerprint\n"))
		}

		if(!is.null(linecols) && (length(linecols) < length(levels(sampleClasses(object))))) {
			stop(paste("There are more classes in this data set than colors specified\n",
			           "you need to specify at least the number of classes colors\n"))
		}

		if (is.null(linecols)) {
			linecols = as.numeric(object@sampleClasses)
	  	} else {
	  		tmp = as.numeric(object@sampleClasses)
  			linecols = linecols[tmp]
  		}
  	} else if (is.null(linecols)) {
		linecols = rep("black", length.out=nfp)
	} else 	if (length(linecols) < nfp) {
		linecols = rep(linecols, length.out=nfp)
	}
	
	if (!is.null(highlight)) {
		linecols = highlightColors(linecols, highlight)
		zorder = c(which(!(1:nfp %in% highlight)), highlight)
	} else {
		zorder = 1:nfp
	}
	xycoords = par("usr")
	for(i in zorder) {
       	lines(ind, data[i,], col=linecols[i], ...)
		points(ind, data[i,], col=linecols[i], pch=20, cex=0.5, ...)
		if (is.flowFPPlex(object) && length(object@fingerprints) > 1) {
			lens = sapply(object@fingerprints, nFeatures)
			acc = 0;
			for(i in 1:length(lens)) {
				acc = acc + lens[i]
				lens[i] = acc
			}
			length(lens) = length(lens) - 1
			segments(lens, xycoords[3], lens, xycoords[4], col='lightblue')
		}
	}
}


