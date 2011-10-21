\name{flowFPPlex-class}
\docType{class}
\alias{flowFPPlex-class}

\title{Fingerprint collection class description.}

\description{This class holds a collection of \link{flowFP}s.}

\section{Objects from the Class}{
Objects should be created by calling the constructor function \link{flowFPPlex}.
}
\section{Slots}{
	 \describe{
    \item{\code{fingerprints}:}{Holds a list of flowFP objects.}
  }
}
\section{Methods}{
  \describe{
    \item{counts}{Concatenates all of the counts for all of the fingerprints in
                    the set.}
    \item{length}{returns the number of fingerprints in the set.}
    \item{nInstances}{returns the number of flowFrames used to create this
                           fingerprint.}
    \item{nFeatures}{returns the total number of features from the collection of
                     fingerprints.}
	 }
}

\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}

\seealso{
	\link{flowFPPlex} - Constructor.
}
\examples{
showClass("flowFPPlex")
}
\keyword{classes}
