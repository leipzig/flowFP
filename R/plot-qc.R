##
## Package: flowFP
## File: plot-qc.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

#
#	R function to plot a group of fingerprints in a grid with qc colors.
#
plotqcFP <- function (fp, main="Fingerprint Deviation Plot",
                        transformation=c("log2norm", "raw", "normalized"),
                        vert_scale=3, method=c("sd", "max"), red_limit=1.0,
                        respect=FALSE)  {

	transformation = match.arg(transformation)
	method = match.arg(method)

	nplots <- nInstances(fp)
	ylim = getYlim(fp, transformation=transformation , vert_scale=vert_scale)
	xlim = getXlim(fp)

	opar <- par(no.readonly=TRUE)
	on.exit (par(opar))
  
	par(mar=c(0,0,0,0))
  	gridLayout(fp, respect=respect)
  	# plot the title info
	plot ("", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', bty='n', xlab="", ylab="")
	text (x=.5, y=.6, label=main, pos=3, cex=2)
	text (x=.25, y=.3, label=sprintf("method = %s", method), pos=4, cex=1.5)
	if (transformation == "log2norm")
		text (x=.25, y=.1, label=sprintf("vertical scale factor = %.1f", vert_scale), pos=4, cex=1.5)
  	
  	draw_wedge(0.65, 0.2, 0.85, 0.3, red_limit) 
  	
	fpMatrix = counts(fp, transformation=transformation)
	fpMatrix[which(is.infinite(fpMatrix))] = NA
	
	if (method == "sd") {
		qcval = apply(fpMatrix, 1, na.rm=TRUE, sd)
	} else {
		qcval = apply(fpMatrix, 1, na.rm=TRUE, max)
	}
	
	
	

	for (i in 1:nplots) {
		bgcol = getCellColor(qcval[i], red_limit)
		initFPplot(fp, ylim=ylim, xlim=xlim, transformation=transformation, 
                       sampleName=sampleNames(fp)[i],
                       xaxt='n', yaxt='n', qcval=qcval[i], bgcol=bgcol)
        lines(fpMatrix[i,])
	}
}
