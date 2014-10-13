\name{is.flowFP}
\alias{is.flowFP}
\title{ Test to see if object is a flowFP }
\description{
  This is a convenience function to check if the object is a \code{flowFP}.
}
\usage{
is.flowFP(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{returns TRUE if the object is a \code{flowFP}}
}

\value{
  returns TRUE if the object is a \code{flowFP}
}

\seealso{ \code{\link{flowFP}}}
\examples{
	data(fs1)
	fp = flowFP(fs1)
	is.flowFP(fp)
}

