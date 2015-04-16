\name{nRecursions-methods}
\docType{methods}
\alias{nRecursions}
\alias{nRecursions<-}
\alias{nRecursions<-,flowFPModel,numeric-method}
\alias{nRecursions<-,flowFP,numeric-method}
\alias{nRecursions-methods}
\alias{nRecursions,flowFP-method}
\alias{nRecursions,flowFPModel-method}
\alias{nRecursions,flowFPPlex-method}
\title{Methods to set or get nRecursions}
\description{
 Gets or Sets the number of recursions for an object
}
\section{Methods}{
\code{nRecursions(object)}
\describe{

\item{object = "flowFP"}{Sets and Gets the number of recursions for viewing a 
      \code{flowFP} object.}

\item{object = "flowFPModel"}{Gets the number of recursions for viewing a 
      \code{flowFPModel} object. }

\item{object = "flowFPPlex"}{Gets the number of recursions for viewing a 
      \code{flowFPPlex} object. }
}}

\section{Notes}{
When a model is computed (either with the \link{flowFPModel} constructor, or implicitly wit
the \link{flowFP} constructor)
the number of bins into which events are divided is determined.  Because the binning process
is recursive, given a high-resolution binning, any lower-resolution representation can
easily be derived.  The set methods do this derivation, and the get methods return the
current resolution of the object.

Thus, the maximum value to which \code{nRecursions} may be set is determined by the value
of \code{nRecursions} at construction time.  If an attempt is made to set \code{nRecursions}
higher than this value a warning will result.

}
\keyword{methods}

