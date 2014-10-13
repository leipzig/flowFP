\name{is.flowFPModel}
\alias{is.flowFPModel}
\title{ Test to see if object is a flowFPModel }
\description{
  This is a convenience function to check if the object is a \code{flowFPModel}.
}
\usage{
is.flowFPModel(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{returns TRUE if the object is a \code{flowFPModel}. }
}

\value{
  Returns TRUE is the object is a \code{flowFPModel}.
}

\seealso{ \code{\link{flowFPModel}} }
\examples{
	data(fs1)
	fp = flowFPModel(fs1)
	is.flowFPModel(fp)
}

