\name{fs1}
\alias{fs1}
\docType{data}
\title{Example FCS flowSet}
\description{
An example FCS flowSet.
}
\usage{data(fs1)}
\format{
 This is an object of class \emph{flowSet} composed of 7 
  \emph{flowFrames}.
}
\details{
  The \emph{flowSet} corresponds to a Leukemia/Lymphoma panel for
  a sample.  Each \emph{flowFrame} correponds to one tube.  The
  first one is an isotype control, and the rest are stained with
  different antibody cocktails.  Parameter 5 (FL3 Log) is notable
  in that it is CD45-ECD for all of the tubes.  This makes it possible
  to gate the entire collection of tubes using CD45 vs SSC.  This
  set has a very similar distribution of CD45 vs SSC across all 7 tubes.

  column names:
  FS Lin SS Log FL1 Log FL2 Log FL3 Log FL4 Log FL5 Log
}
\source{
 Data were generously provided by Clarient, Inc.
}
\references{
  W.T. Rogers and H. A. Holyst, \dQuote{flowFP: A Bioconductor Package for Fingerprinting
  Flow Cytometric Data}, \emph{Submitted to: Advances in Bioinformatics, April 1, 2009}.
}
\examples{
data(fs1)
}
\keyword{datasets}
