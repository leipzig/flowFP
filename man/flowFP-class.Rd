\name{flowFP-class}
\docType{class}
\alias{flowFP-class}

\title{Fingerprint class description.}

\description{
This class represents fingerprints generated from a model derived
from a \link[flowCore:flowFrame-class]{flowFrame} or 
\link[flowCore:flowSet-class]{flowSet}.
}

\section{Objects from the Class}{
Objects are created by calling the constructor \link{flowFP}.
}

\section{Slots}{
	 \describe{
    \item{\code{counts}:}{This matrix is organized such that rows correspond 
                          to individual \emph{flowFrame}s, and columns to 
                          the bins in the underlying \code{flowFPModel}.  A
                          matrix element thus represents the number of events in
                          a single bin from a single \emph{flowFrame}.
                          The sum of each row will equal the total number of events in
                          the corresponding \emph{flowFrame}.
     }
    \item{\code{tags}:}{This list specifies the bin number to which each event is assigned.
                        There is a list element for each instance in the object.  Each
                        list element is a integer vector whose length is the number of
                        events in the corresponding \emph{flowFrame}.
     }
    \item{\code{sampleNames}:}{Names of instances, taken from the flowSet on
                               construction.}
    \item{\code{sampleClasses}:}{Factor object describing the category of each 
                                 instance. (e.g. "Normal" or "Diseased")}
	\item{\code{name}:}{Inherited from \link{flowFPModel}.}
    \item{\code{parameters}:}{Inherited from \link{flowFPModel}.}
    \item{\code{nRecursions}:}{Inherited from \link{flowFPModel}.}
    \item{\code{trainingSet}:}{Inherited from \link{flowFPModel}.}
    \item{\code{trainingSetParams}:}{Inherited from \link{flowFPModel}.}
    \item{\code{dequantize}:}{Inherited from \link{flowFPModel}.}
    \item{\code{split_val}:}{Inherited from \link{flowFPModel}.}
    \item{\code{split_axis}:}{Inherited from \link{flowFPModel}.}
    \item{\code{binBoundary}:}{Inherited from \link{flowFPModel}.}
    \item{\code{.cRecursions}:}{Inherited from \link{flowFPModel}.}
    \item{\code{.tmp_tags}:}{Inherited from \link{flowFPModel}.}
  }
}
\section{Extends}{
Class \linkS4class{flowFPModel}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{[}}{Used to select a subset of a fingerprint object. }
    \item{\code{counts}}{ Returns the counts matrix of a fingerprint object.  See \link{counts}. }
    \item{\code{length}}{ Returns the number of instances in the fingerprint object. }
    \item{\code{nbins}}{ Returns the number of bins (equal to 2^nRecursions) in the fingerprint object. }
    \item{\code{sampleClasses}}{ Returns the sampleClasses slot of the fingerprint object. }
    \item{\code{sampleClasses<-}}{ Set method for the sampleClasses
					slot of the fingerprint object.  Note that the length of the provided object must be
					equal to the number of instances in the fingerprint object.  If the provided object
					is not a factor it will be coerced to one. }
    \item{\code{tags}}{ Returns the \code{tags} slot.  See \link{tags}. }
    \item{\code{show}}{ Shows the contents of the fingerprint object.}
    \item{\code{summary}}{Gives a summary of the fingerprint object.}
	 }
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
	model <- flowFPModel(fs1, parameters=c(4,5), nRecursions=6)
	
	fp <- flowFP(fs1, model)
	
	plot(fp)
	
}
\keyword{classes}
