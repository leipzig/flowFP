##
## Package: flowFP
## File: plot-fancy.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

plotFancy <- function(x, y, highlight=NULL, parameters=NULL, 
                      transformation=c("raw", "normalized", "log2norm"), 
                      showbins=NULL, showfp=TRUE, alpha=.1, border='gray', linecols=NULL, 
                      pch='.', cex=1, main="Fingerprints", ...) {

	fp = x
	nfp = nInstances(fp)
	fs = y
	transformation = match.arg(transformation)
	
	check_plot_args(fp, fs=fs, highlight=highlight, showbins=showbins)
	
	opar <- par(no.readonly=TRUE)
	on.exit (par(opar))
	
	if (showfp) {
		layout( matrix(c(1,2), nrow=2, ncol=1), heights=c(1, 0.5))
	}
	
	
	## Figure out if the user wants to see the bin boundaries, and create vector
	## for user short hand 'all'. Note if we are to show bins we plot all of the other
	## points in black.
	if(!is.null(showbins) ) {
		if ( (is.logical(showbins) && showbins == TRUE) || showbins == 'all') {
			showbins = 1:nFeatures(fp)
		}
	}

	# create luts
	satLut=binColorLut(nFeatures(fp))
	if(is.null(highlight))
		semiSatLut = satLut
	else
		semiSatLut = binColorLut(nFeatures(fp), saturation=0.25)

	dotcols = list()
	if (is.numeric(showbins)) {
		for(i in 1:nfp) {
			dotcols[[i]] = "black"
		}

	} else {
		tags = tags(fp)
		for(i in 1:nfp) {
			if (i %in% highlight) {
				dotcols[[i]] = satLut[tags[[i]]]
			} else {
				dotcols[[i]] = semiSatLut[tags[[i]]]
			}
		}
	}

	zorder = c(which(!(1:nfp %in% highlight)), highlight)

	parameters = getDefaultParameters(x, parameters)
		
	# This is a work around, I truly don't understand.
	# I think that the dispatcher should handle calling plot()
	# with a 'flowFrame', and the 'parameter names'  perfectly but it 
	# doesn't work. So instead I force the dispatcher to get me the specific
	# plot function, and then call it.
	ffplot = getMethod("plot", signature("flowFrame", "character"))

	first_time = TRUE
	for(i in zorder) {
		if (first_time) {
			ffplot(fs[[i]], colnames(fs[[i]])[parameters] , smooth=FALSE, 
			       col=dotcols[[i]], main=main)
		} else {
			x = exprs(fs[[i]])[,parameters[1]]
			y = exprs(fs[[i]])[,parameters[2]]
			points(x, y , col=dotcols[[i]], pch='.')
		}
		first_time = FALSE
	}
	
	if(!is.null(showbins) ) {
		if ( (is.logical(showbins) && showbins == TRUE) || showbins == 'all') {
			showbins = 1:nFeatures(fp)
		}
	}
	if (is.numeric(showbins)) {
		plotBinBoundaries(fp, y=fs, parameters, alpha, border, showbins, pch=pch, cex=cex)
	} 
	
	if (showfp) {
		mar = par("mar")
		mar[3] = 0.0
		par(mar=mar)

		plotTangleFP(fp, transformation=transformation, linecols=linecols, 
		             highlight=highlight, main="", ...)
	 	par(xpd=TRUE)
 		usr = par("usr")
	 	color_wedge_h = (usr[4] - usr[3]) * 0.05

		ypos = seq(from=usr[3] - color_wedge_h, to=usr[3] + color_wedge_h, length.out=2)
		xpos = seq(from=0.5, to=length(satLut) + .5, length.out=length(satLut) + 1)
  		wedge = matrix(rep(c(1:length(satLut)), 1), ncol=length(satLut), nrow=1, byrow=TRUE)

		image(xpos, ypos, t(wedge), col=satLut, add=TRUE, xlab=NA, ylab=NA, yaxt="n" )
	}
	
}
