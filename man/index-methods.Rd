\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{[,flowFP,numeric-method}
\alias{[,flowFP,character-method}
\alias{[,flowFPPlex,ANY-method}
\alias{[[-methods}
\alias{[[,flowFPPlex,ANY-method}
\alias{[[<-,flowFPPlex,numeric,ANY,flowFP-method}

\title{ Indexing methods for flowFP and flowFPPlex}
\description{
 These methods allow subset indexing for \code{flowFP} and \code{flowFPPlex} objects.
}
\section{Methods}{
\describe{

\item{x = "flowFP", i = "numeric"}{ [ indexing returns a subset of the fingerprints in a flowFP. }

\item{x = "flowFP", i = "character"}{ [ indexing returns a subset of the fingerprints in a flowFP.
If i matches one or more \emph{sampleClasses} in the flowFP, the subset operator returns all fingerprints 
of that class (or classes). If no matches with \emph{sampleClass} are found, then the operator attempts
to match \emph{sampleNames}, and returns all matches found.}

\item{x = "flowFPPlex", i = "numeric"}{ [ indexing returns a new \code{flowFPPlex} with the specified subset of 
                                          fingerprints. }

\item{x = "flowFPPlex", i = "numeric"}{ [[ indexing returns a single \code{flowFP} object from a flowplex.}

}}
\keyword{methods}
