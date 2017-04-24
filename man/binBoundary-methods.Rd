\name{binBoundary-methods}
\docType{methods}
\alias{binBoundary}
\alias{binBoundary-methods}
\alias{binBoundary,flowFP-method}
\alias{binBoundary,flowFPModel-method}
\title{Return bin boundaries for a flowFPModel}

\description{
	These methods return the bin boundaries as a list.  The number of list elements
	is equal to the number of features in the fingerprint.
	Each boundary has two points, lower left (called \code{ll}) and upper right
	(called \code{ur}). The dimensionality of \code{ll} and \code{ur} is the
	dimensionality of the \code{flowFrame} or \code{flowSet} used to construct
	the object.
}


\section{Methods}{
\code{binBoundary(object)}
\describe{
\item{object = "flowFP"}{ This method returns the bin boundaries from the underlying
   \code{flowFPModel}. }

\item{object = "flowFPModel"}{ This Method returns the \code{flowFPModel} bin boundaries. }
}}
\keyword{methods}
