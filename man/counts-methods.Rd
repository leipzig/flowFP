\name{counts-methods}
\docType{methods}
\alias{counts}
\alias{counts,flowFP-method}
\alias{counts,flowFPPlex-method}
\title{Method for getting fingerprint data.}
\description{
 These read-only methods retrieve data out of a fingerprint in three
 different transformations: \dQuote{raw}, \dQuote{normalized}, \dQuote{log2norm}.
}


\section{Methods}{
\code{counts(object, transformation=c("raw", "normalized", "log2norm")}
\describe{

\item{object = "flowFP"}{ This method returns a feature vector for each instance
      in the \code{flowFP}, the user can change the way the features are represented,
      by using the transformation argument. see below. }

\item{object = "flowFPPlex"}{For a \code{flowFPPlex} the code assembles all of the
      counts matrices from all of the \code{flowFP} contained in the \code{flowFPPlex}.}
}}

\details{
 Transformations: There are three representations of the data that the user has access to,
 "raw" is the number of events in each bin, "normalized" is the ratio of the number of
 events divided by the expected number (the expected number of events is calculated
 by dividing the total number of events in the instance used to make the fingerprint
 by the number of features in the fingerprint). "log2norm" is the log2 of the 
 normalized values.
}
\keyword{methods}
