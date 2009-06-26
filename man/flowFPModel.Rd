\name{flowFPModel}
\alias{flowFPModel}
\title{Fingerprint model constructor}
\description{
A constructor for objects of type \link{flowFPModel-class}.
	
}
\usage{
flowFPModel(fcs, name="Default Model", parameters=NULL, nRecursions="auto", 
            dequantize=TRUE, sampleSize=NULL, excludeTime=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcs}{Training data for model, either a \link[flowCore:flowFrame-class]{flowFrame} or 
                    \link[flowCore:flowSet-class]{flowSet}.}
  \item{parameters}{A vector of parameters to be considered during model 
                    construction.}
  \item{nRecursions}{Number of times the FCS training data will be subdivided. 
                    Each level generated \eqn{2^{nRecursions}}{2^nRecursions} bins. A warning will be generated
                    if the number of expected events in each bin is < 1. (e.g. if
                    your training set had 1000 events, and you specified \code{level=10}.)
                    }
  \item{dequantize}{If TRUE, all of the events in the training set will be made 
                    unique by adding a tiny value (proportional to the ordinal
                    position of each event) to the data.}
  \item{sampleSize}{Used to specify the per-\emph{flowFrame} sample size of the data 
                    to use in model generation. If \emph{NULL}, all of the data in
                    \emph{x} is used.  Setting this to a smaller
                    number will speed up processing, at the cost of accuracy.}
  \item{name}{A descriptive name assigned to the model.}
  \item{excludeTime}{If TRUE (default) and no parameters are specified, the model 
                     constructor will atempt to exclude the time parameter from
                     consideration when building the model. The time parameter is
                     identified as have the label 'time', matched ignoring case.}

}
\details{
This function is used to create a \emph{flowFPModel}, which can then be used to
create a set of fingerprints using \link{flowFP}. A model is a
representation of the multivariate probability density function for the training
set \emph{x}.  This representation is in the form of a set of \emph{bins}, each
of which contains (nearly) the same number of events in the training set.  Thus,
the model can be regarded as a multivariate histogram of \emph{x} with fixed
bin count and variable bin size.  Bins will be smallest in regions of high 
density, and largest in regions that are sparsely populated.

Bins are constructed by recursively subdividing the multivariate space specified
by \emph{parameters}.  At the first level, the entire space is divided in
half in the direction of the parameter with the highest variance.  At the next
level, each of the halves from the first level is itself divided in half,
again in the direction of the parameter whose variance (among the events in the
parent bin) is the highest.  Thus, for each level the number of resulting bins
is doubled.  The number of bins in the model is thus \eqn{2^{nRecursions}}{2^nRecursions}. 

A model is of little interest by itself.  Its only utility is to provide a sort
of basis function for forming fingerprints.  It may be interesting in some cases
to use the entire set \emph{x} to build the model, and then to generate fingerprints
for each instance in \emph{x}.  For this reason, calling \link{flowFP} without
a model causes \emph{flowFPModel} to be invoked silently under the hood in order
to create a model of the \emph{flowSet} provided to \link{flowFP}.

}

\note{
If a \link[flowCore:flowSet-class]{flowSet} is provided as input, it is internally
collapsed to a \link[flowCore:flowSet-class]{flowFrame} for binning.
}

\value{
  An object of type \link[flowFPModel-class]{flowFPModel} is returned.
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
	\link{flowFP}, \link{flowFPModel}
}
\examples{

	# load a flowSet to use as an example.
	library(flowFP)
	data(fs1)
	fs1
	model <- flowFPModel(fs1, parameters=c(4,5), nRecursions=6)
	
	fp <- flowFP(fs1, model)
	
	plot(fp)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

