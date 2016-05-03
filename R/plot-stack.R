##
## Package: flowFP
## File: plot-stack.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

plotStackFP <-function(x, transformation=c("raw", "normalized", "log2norm"), 
                        linecols=NULL, useClasses=FALSE, vert_scale=3,
                        ylim=NULL, main="Fingerprints",...){

	transformation = match.arg(transformation)
	if (is.null(ylim))
		ylim = getYlim(x, transformation, vert_scale)

	num_plots = nInstances(x)
	opar <- par(no.readonly=TRUE)
	on.exit (par(opar))

	nplot_rows = num_plots + 2 # add a top and bottom slice for title and axis
	plot_grid = matrix(1:nplot_rows, nrow=nplot_rows, ncol=1)
	heights = c(0.75, rep(1.0, num_plots), 0.75)
	nf<-layout(plot_grid, heights=heights)

	# Margins
	mar = par("mar") # get current margins
	# bottom, left, top, right.
	mar[1] = 0.2
	mar[3] = 0.2
	par(mar=mar)

	plot('', xlim=c(0,1), ylim=c(0,1), ylab="",	xlab="", bty="n", xaxt="n",
	     yaxt="n", pch="")
	text(0.5, 0.5, adj = c(0.5, 0.5), labels=main, cex=2)
	
	for(i in 1:num_plots) {
		fp = x[i]
		plot('', xaxt="n", type='l', ylim = ylim, 
		xlim = c(0, nFeatures(fp)), xlab="",  ylab="")
		draw_tangle_lines(fp, transformation, linecols, NULL, useClasses=useClasses, ...)
		legend("top", sampleNames(fp), bty="n", cex=1.0)
	}
	plot('', xaxt="n", yaxt='n', pch="", ylim = c(0, 1),
	     xlim = c(0, nFeatures(fp)), xlab="", bty='n', ylab="")
	axis(side=1, pos=1)
}
