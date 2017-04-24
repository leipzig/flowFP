\name{hasClasses-methods}
\docType{methods}
\alias{hasClasses}
\alias{hasClasses-methods}
\alias{hasClasses,flowFP-method}
\alias{hasClasses,flowFPPlex-method}
\title{ Returns TRUE if the flowFP has a class factor for its instances }
\description{
 This read-only accessor returns TRUE if the data object \link{flowFP} or 
 \link{flowFPPlex} instances have been assigned a factor object representing
 the class membership of the instances.
}
\section{Methods}{
\code{hasClasses(object)}
\describe{
\item{object = "flowFP"}{Returns TRUE if there are class assignments for each
      instance in this \code{flowFP}.}

\item{object = "flowFPPlex"}{ If a \code{flowFPPlex} has classes associated with its 
      instances, they must all be the same.}
}}
\keyword{methods}
