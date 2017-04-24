\name{plate}
\alias{plate}
\docType{data}
\title{96 well plate data.}
\description{
 An example of 96 well FCS flow data.
}
\usage{data(plate)}
\format{
 This is an object of class \emph{flowSet} comprising 96 
  \emph{flowFrames}.
  Each \emph{flowFrame} correponds to one sample.
  column names:
  FSC-H SSC-H FL1-H FL2-H FL3-H FL1-A FL4-H Time
}
\details{
  This data set is derived from the references below.  Each well in a 96-well
plate is described by a \emph{flowFrame}.  Data were drastically down-sampled
to 1000 events per \emph{flowFrame} so that they could be included in 
this package.  Note that flourescence parameters 5 (CD4 PerCP Cy5.5) and
7 (CD3 APC) are common across the entire set.
}
\source{
 Original, non-sampled data are available at
\url{http://www.ficcs.org/software.html#Data_Files}.
}
\references{
  M. Inokuma, C. dela Rosa, C. Schmitt, P. Haaland, J. Siebert, D. Petry, M. Tang, M. A.
 Suni, S. A. Ghanekar, D. Gladding, J. F. Dunne, V. C. Maino, M. L. Disis, and H. T.
 Maecker. \dQuote{Functional t cell responses to tumor antigens in breast cancer patients have a
 distinct phenotype and cytokine signature}. \emph{J Immunol}, 179(4):2627-33, Aug 15 2007.

M. Inokuma, C. dela Rosa, C. Schmitt, P. Haaland, J. Siebert, D. Petry, M. Tang, M. A. Suni,
 S. A. Ghanekar, D. Gladding, J. F. Dunne, V. C. Maino, M. L. Disis, and H. T. Maecker,
 2008.  Data available at \url{http://www.ficcs.org/software.html#Data_Files}.

}
\examples{
data(plate)
}
\keyword{datasets}
