##
## Package: flowFP
## File: plot-grid.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

#
#	R function to plot a group of fingerprints in a grid.
#
plotGridFP <- function (fp, vert_scale=3, main="Fingerprints", linecols="black",
                        transformation=c("raw", "normalized", "log2norm"),
                        respect=FALSE) {

	transformation = match.arg(transformation)

	nplots <- nInstances(fp)
	linecols <- rep(linecols, nplots)
	ylim = getYlim(fp, transformation, vert_scale)
	xlim = getXlim(fp)

	opar <- par(no.readonly=TRUE)
	on.exit (par(opar))
  
	par(mar=c(0,0,0,0))
  	gridLayout(fp, respect=respect)
                       
  	plot('', xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', ,xlab='', ylab='', bty='n')
  	if (!is.null(main)) {
		legend("center", main, bty="n", cex=3)
	}

	fpMatrix = counts(fp, transformation=transformation)
	fpMatrix[which(is.infinite(fpMatrix))] = NA
	
	
	for (i in 1:nplots) {
		initFPplot(fp, ylim=ylim, xlim=xlim, transformation=transformation, 
                       sampleName=sampleNames(fp)[i],
                       xaxt='n', yaxt='n')
        lines(fpMatrix[i,], col=linecols)
	}

}


