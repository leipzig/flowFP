##
## Package: flowFP
## File: plot-utils.R
## Author: Herb Holyst
##
##  Copyright (C) 2009 by University of Pennsylvania,
##  Philadelphia, PA USA. All rights reserved. 
##

plateLayout <- function (fp) {
	# how many fingerprints are we plotting?
	nplots = nInstances(fp)
	# if nplots is 96, let's layout as 8 x 12
	if (nplots == 96) {
		initialize_96wellPlate()
	} else {
	    # if nplots isn't 96, we'll make the array as square as possible
		nrow = floor (sqrt(nplots))
		ncol = round(nplots/nrow + 0.499)
		plotPos = vector(mode="integer", length=nrow * ncol)
		plotPos[1:nplots] = (1:nplots) + 1
	    layoutMat = matrix (nrow=nrow + 1, ncol=ncol)
	    layoutMat[1,] = 1
	    layoutMat[2:(nrow+1),] = matrix(plotPos, nrow=nrow, ncol=ncol, byrow=TRUE)
	    layout (layoutMat)
	}
	
}

initialize_96wellPlate <-function() {
  	
	layoutMatrix = matrix(
	c( 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
	    0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
	   13, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
	   14, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
	   15, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
	   16, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
	   17, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
	   18, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93,
	   19, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,
	   20,106,107,108,109,110,111,112,113,114,115,116,117), 
	   ncol=13, nrow=10, byrow=TRUE)
	heights = c(1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
	widths = c(0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
	nf<-layout(layoutMatrix, heights=heights, widths=widths, respect=TRUE)
	rowLabels = c("A", "B", "C", "D", "E", "F", "G", "H")
	colLabels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
	labels = c(colLabels, rowLabels)
	
	for(i in 1: length(colLabels)) {
		plot ("", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', bty='n', xlab="", ylab="")
		text (x=0.5, y=0, label=colLabels[i], pos=3, cex=2)
	}
	for(i in 1: length(rowLabels)) {
		plot ("", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', bty='n', xlab="", ylab="")
		text (x=1, y=.5, label=rowLabels[i], pos=2, cex=2)
	}

}


getCellColor <- function(x, red_thresh) {
  # we divide the color space up into three zones
  saturation_low = 0.5
  value = 1.0
  
  if(x >= red_thresh)
  	return (hsv(0,1,1))
  
  t = x / red_thresh
  if (t < .333) { # zone1
  	saturation = 1.0 - (1.0 - saturation_low) *  t / .333
  } else if (t < .666) {
  	saturation = saturation_low
  } else {
  	saturation = saturation_low + (1.0 - saturation_low) * (t - .666) / (.333)
  }
  
  hue = 0.3333 * (1.0 - t) # hue ranges from green(0.333) to red (0)
  return(hsv(hue, saturation, value))
}


draw_wedge <- function(x1, y1, x2, y2, red_limit) {
  nx = 100
  ny = 1
  ypos = seq(from=y1, to=y2, length.out=ny + 1)
  xpos = seq(from=x1, to=x2, length.out=nx + 1)
  wedge = matrix(rep(c(1:nx), ny), ncol=nx, nrow=ny, byrow=TRUE)
  col = vector(len=nx)
 
  for (i in 1:nx) {
    col[i] = getCellColor (red_limit*i/nx, red_limit)
  }
  image(xpos, ypos, t(wedge), col=col, add=TRUE)
  text(x1, y1, "0.0", pos=1)
  text (x2, y1, sprintf("%.1f", red_limit), pos=1)
  rect (x1, y1, x2, y2)
}


binColorLut <- function(ncolors, colsneeded, saturation=1.0, alpha=1.0, negative=FALSE) {

	if (missing(colsneeded))
		colsneeded = 1:ncolors
	num_val = 4
    vlo = 0.35
    hue = vector(length=ncolors)
    val = vector(length=ncolors)
    for(i in 1:(ncolors/num_val)) {
    	begin = (i-1) * num_val + 1
    	end = begin + num_val - 1
    	hue[begin:end] = .333 + (i /(ncolors / num_val)) * .5
    	val[begin:end] = seq(vlo, 1, length.out=num_val)
    }
    if (negative)
    	hue = 1.0 - hue
    return (hsv(h=hue[colsneeded], v=val[colsneeded], s=saturation, alpha=alpha))
}

# currently un-used TODO: delete this if we don't need it.
replaceColorAlpha <- function(colors, alpha=0.25) {

	if (length(colors) > 0) {
	
		for(i in 1:length(colors)) {
			tmp = col2rgb(colors[i])
			tmp = rgb2hsv(r=tmp[1], g = tmp[2], b = tmp[3])
			colors[i] = hsv(h=tmp[1], s=tmp[2], v=tmp[3], alpha = alpha)
		}
	}	
	return (colors)
	
}
##
## Initialize fingerprint/plex plot area
##
initFPplot <- function(object, ylim=NULL, xlim=NULL, xlab="", ylab="", main=NULL, 
                       transformation=c("raw", "normalized", "log2norm"), 
                       bgcol=NULL, qcval=NULL, sampleName=NULL, vert_scale=3, ...) {

	transformation = match.arg(transformation)

	if(is.null(ylim))
		ylim=getYlim(object, transformation=transformation, vert_scale)
		
	if (is.null(xlim))
		xlim=getXlim(object)
		
	plot('', type='l', ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, main=main, ...)

	xycoords = par("usr")
	
	if (!is.null(bgcol)) {
	    rect (xycoords[1], xycoords[3], xycoords[2], xycoords[4], col=bgcol)
	}
	
	if (is.flowFPPlex(object) && length(object) > 1) {
		lens = c(0, sapply(object@fingerprints, nFeatures))
		acc = 0;
		for(i in 1:length(lens)) {
			acc = acc + lens[i]
			lens[i] = acc
		}

		segments(lens, xycoords[3], lens, xycoords[4], col='lightblue')
	}

	horz = NULL
	if (transformation == "normalized") {
		horz = 1
	} else if(transformation == "log2norm") {
		horz = 0
	}

	if(!is.null(horz)) {
		segments(xycoords[1], horz, xycoords[2], horz, col="lightblue")
	}
	
	
	if(!is.null(sampleName)) {
		legend("top", sampleName, bty="n")
	}
	
	if(!is.null(qcval)) {
		ypos = 0.75 * (xycoords[4]-xycoords[3]) + xycoords[3]
		xpos = 0.2 * (xycoords[2]-xycoords[1]) + xycoords[1]
		text (xpos, ypos, sprintf ("%.2f", qcval), cex=1.0)
	}
		
}

getYlim <- function(object, transformation, vert_scale) {
	
	if (transformation == "log2norm")
		return(c(-vert_scale, vert_scale))
		
	data = counts(object, transformation=transformation)
	return( range(data, na.rm = TRUE, finite = TRUE))
}

getXlim <- function(object) {
	return (c(1,nFeatures(object)))
}

getDefaultParameters <- function(x, parameters=NULL) {
	if (is.null(parameters)) {
		tmp = unlist(x@split_axis)
		used = which(x@parameters %in% x@trainingSetParams[tmp])
        parameters = which(x@trainingSetParams %in% x@parameters[used])
	}
	length(parameters) = 2
	return (parameters)
}

gridLayout <- function (fp, respect=TRUE) {
	# how many fingerprints are we plotting?
	nplots = nInstances(fp)
    # we'll make the array as square as possible
	nrow = floor (sqrt(nplots))
	ncol = round(nplots/nrow + 0.499)
	layoutMat = matrix (c(rep(1, ncol), 2:((nrow * ncol) + 1)), 
	                    nrow=nrow + 1, ncol=ncol, byrow=TRUE)
	blank_column = rep(0, nrow + 1)

	layoutMat = cbind(blank_column, layoutMat, blank_column)
	blank_row = rep(0, ncol(layoutMat))
	layoutMat = rbind(layoutMat, blank_row)
	widths = c(0.25, rep(1, ncol), 0.25)
	heights = c(.5, rep(1, nrow), 0.25)
	layout (layoutMat, widths=widths, heights=heights, respect=respect)
	
}

check_plot_args <- function(fp, fs=NULL, highlight=NULL, showbins=NULL) {

	if (!is.null(fs) && !all(sampleNames(fp) == sampleNames(fs))) {
		stop("Fingerprints, and the flowSet must contain the same sample names.\n",
		     call.=FALSE)
	}
	
	if (!is.null(highlight) && !all(highlight %in% 1:nInstances(fp))) {
		hs = paste(highlight, collapse=", ")
		stop(paste("'highlight' is used to highlight instances from a flowFP\n",
		           "There are only ", nInstances(fp), " in this fingerprint\n",
		           "You can't highlight ", hs, "\n", sep=""), call.=FALSE)
	}
	
	if (!is.null(showbins) && !all(showbins %in% 1:nFeatures(fp))) {
		sbs = paste(showbins, collapse=", ")
		stop(paste("'showbins' is used to show selected bins from a flowFP\n",
		           "There are only ", nFeatures(fp), " in this fingerprint\n",
		           "You can't show -> ", sbs, "\n", sep=""), call.=FALSE)
	}
	
	
		
}
