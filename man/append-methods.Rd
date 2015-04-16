\name{append-methods}
\docType{methods}
\alias{append-methods}
\alias{append,flowFPPlex,flowFP-method}
\alias{append,flowFPPlex,list-method}
\alias{append,flowFPPlex,flowFPPlex-method}

\title{ Methods to append flowFP fingerprints to a flowFPPlex }


\description{
  
  This method joins one or more \code{\link{flowFP}} objects into a
  \code{\link{flowFPPlex}}. Conceptually, each instance's fingerprint feature
  vector is extended. An error will occur if the sample names and the class 
  names differ. If we imagine the fingerprints as a matrix where columns are
  features and the rows are instances, this method binds the columns
  onto the existing fingerprint matrix. 

}

\section{Methods}{

  \code{append(x, values, after=length(x))} 

  \code{\link{flowFPPlex}} methods:

  \describe{
  
  \item{x = "flowFPPlex", values = "flowFP", after = length(x@fingerprints)}{
        Appends a fingerprint onto a flowFPPlex. }

  \item{x = "flowFPPlex", values = "list", after = length(x@fingerprints)}{
        Appends a list of fingerprints to a flowFPPlex. }
          
  \item{x = "flowFPPlex", values = "flowFPPlex", after = length(x@fingerprints)}{
        Combines two flowFPPlex-es into one. }	
  }
}

\details{

  These methods append one or more fingerprints into a flowFPPlex.
}


\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}

\seealso{
	\link{flowFPPlex} \link{sampleNames} \link{sampleClasses}
}
\examples{
	data(fs1)
	fp1 = flowFP(fs1, name="fingerprint 1")
	fp2 = flowFP(fs1, name="fingerprint 2")
	plex = flowFPPlex(fp1)
	plex = append(plex, fp2)
	plex
}

\keyword{methods}
