\name{nFeatures-methods}
\docType{methods}
\alias{nFeatures}
\alias{nFeatures-methods}
\alias{nFeatures,flowFP-method}
\alias{nFeatures,flowFPModel-method}
\alias{nFeatures,flowFPPlex-method}
\title{Returns the number of features in an object.}
\description{
 Theses methods return the number of leaf bins that were generated or will be generated
 by a specific model of plex of models.
}
\section{Methods}{
\code{nFeatures(object)}
\describe{

\item{object = "flowFP"}{The number of bins in a fingerprint.}

\item{object = "flowFPModel"}{The number of bins that will be created by this model.}

\item{object = "flowFPPlex"}{The total number of bins in a plex fingerprint.}
}}
\keyword{methods}
