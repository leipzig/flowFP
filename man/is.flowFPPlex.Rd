\name{is.flowFPPlex}
\alias{is.flowFPPlex}
\title{Test to see if object is a flowFPPlex}
\description{
  This is a convenience function to check if the object is a \code{flowFPPlex}.
}
\usage{
is.flowFPPlex(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ returns TRUE if the object is a \code{flowFPPlex}. }
}

\value{
 returns TRUE if the object is a \code{flowFPPlex}.
}

\seealso{ \code{\link{flowFPPlex}}}
\examples{
data(fs1)
	fp = flowFP(fs1)
	plex = flowFPPlex(fp)
	is.flowFPPlex(fp)
}
