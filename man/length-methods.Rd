\name{length-methods}
\docType{methods}
\alias{length-methods}
\alias{length,flowFPPlex-method}
\alias{length<-,flowFPPlex,numeric-method}

\title{The number of fingerprints in a flowFPPlex }
\description{
 This method sets or gets the number of fingerprints that have been concatinated together
 to form this flowFPPlex. 
}
\section{Methods}{
\code{length(x)}

\describe{
\item{x = "flowFPPlex"}{Sets or gets the number of fingerprints in a flowFPPlex }

}}

\section{Notes}{
The get method returns the number of \emph{flowFP} objects in the plex. \emph{It is
not the number of instances or the number of features in the counts matrix!}

The set method can only truncate the set of \emph{flowFP} objects from the right.
}
\keyword{methods}

