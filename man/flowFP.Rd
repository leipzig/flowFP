\name{flowFP}
\alias{flowFP}
\title{Fingerprint constructor}
\description{
Constructor for \link{flowFP-class}.
}
\usage{
flowFP(fcs, model = NULL, sampleClasses = NULL, ...)
}

\arguments{
  \item{fcs}{ A \link[flowCore:flowFrame-class]{flowFrame} or 
           \link[flowCore:flowSet-class]{flowSet} for which fingerprint(s) are desired.}
  \item{model}{ A model generated with the \link{flowFPModel} constructor, or \emph{NULL}. 
                if \emph{NULL}, a default model will be silently generated from all instances
                in \emph{x}.}
  \item{sampleClasses}{ An optional character vector describing modeling classes.  If
                supplied, there must be exactly one element for each \emph{flowFrame}
                in the \emph{flowSet} in \emph{x} (see Details).}
  \item{\dots}{ If model is \emph{NULL}, additional arguments are passed on to the 
                model constructor. see \link{flowFPModel} for details.}
}

\details{
  A flowFP object is a reduced representation of a 
\link[flowCore:flowFrame-class]{flowFrame} or \link[flowCore:flowSet-class]{flowSet} 
with respect to a \link{flowFPModel}. The model is derived by finding multivariate regions containing
(nearly) equal numbers of events in a training set using one or more user-specified 
parameters (\emph{e.g.} CD45 and Side Scatter).  The resulting \emph{flowFP} represents
the probability (i.e. number of events) of the instance(s) in \emph{x} in each of the bins
in the underlying model.  

This representation is carried in a slot (called \emph{counts})
in the object, a matrix in which each row corresponds ordinally to the \emph{flowFrames}
in \emph{x}.  Note that if \emph{x} is a \emph{flowFrame} and not a \emph{flowSet}, the counts
matrix will have only 1 row.  In any case, each row of the matrix is an individual fingerprint,
where the fingerprint elements (columns) have a one-to-one correspondence with the bins
in the underlying \emph{model}.

Another slot in the \emph{flowFP} object describes the bin index for each
event, allowing the user to select events according to their membership in bins.  Think of this
as "micro-gating".

}
\value{
  An instance of an object of type \emph{flowFP}.
}
\references{
M. Roederer, et. al. (2001) Probability Binning Comparison: A Metric for Quantitating 
Multivariate Distribution Differences, \emph{Cytometry} \bold{45}, 47-55.

W. Rogers et. al. (2008) Cytometric Fingerprinting: Quantitative Characterization of 
Multivariate Distributions, \emph{Cytometry Part A} \bold{73}, 430-441.
}
\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}

\seealso{ 
	\link{flowFP-class}, \link{flowFP}, \link{flowFPModel-class}, \link{flowFPModel}
}
\examples{

	# load a flowSet to use as an example.
	# fs <- read.flowSet(your fcs files...)
	data(fs1)
	model <- flowFPModel(fs1, parameters=c(4,5), nRecursions=6)
	
	fp <- flowFP(fs1, model)
	
	plot(fp)
	
}

\keyword{classes}

