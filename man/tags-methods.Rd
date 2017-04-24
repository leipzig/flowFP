\name{tags-methods}
\docType{methods}
\alias{tags}
\alias{tags-methods}
\alias{tags,flowFP-method}
\title{Returns tags from a flowFP object.}
\description{
 The tags method is a map into a \code{flowFP} object which defines the bin number 
 to which each event in the object has been assigned in fingerprinting process.
}
\section{Methods}{
\code{tags(object)}
\describe{

\item{object = "flowFP"}{ returns a list of vectors.}
}}

\section{Notes}{
This method returns a list of vectors.  There is a one-to-one correspondence between
vectors in the list and \emph{flowFrame}s in the \emph{flowSet} in the \code{flowFP}.

For each vector in the list, there is a one-to-one correspondence between the vector
elements and the events in the corresponding \emph{flowFrame}.  Each vector element
is an integer index that indicates the bin number to which that event was assigned
in the fingerprinting process.
}

\keyword{methods}
