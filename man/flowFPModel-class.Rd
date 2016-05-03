\name{flowFPModel-class}
\docType{class}
\alias{flowFPModel-class}

\title{Fingerprint model class description.}
\description{This class is the fundamental clase for the package.  It contains data and 
methods used to construct a model of the probability density function of a prototype 
dataset provided in the form of a \link[flowCore:flowFrame-class]{flowFrame} or 
\link[flowCore:flowSet-class]{flowSet}.
}
\section{Objects from the Class}{
Objects may only be created by calling the constructor function \link{flowFPModel}.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{A user-supplied descriptive name for the model.}
    \item{\code{parameters}:}{List of FCS parameters to use for model creation.  Can 
                              be specified either by the names of the parameters 
                              or the indices of the parameters.}
    \item{\code{nRecursions}:}{Number of levels of recursive subdivision.
                              The number of bins in the model will equal 
                              \eqn{2^{nRecursions}}{2^nRecursions}.}
    \item{\code{trainingSet}:}{Names of \emph{flowFrames} from the \emph{FlowSet} used to
                               construct the model. }
	\item{\code{trainingSetParams}:}{Names of all of the parameters from the \emph{flowFrames}
	                                 from the \emph{FlowSet} used to construct the model. }
    \item{\code{dequantize}:}{If TRUE, all of the event parameter values in the training 
                              set will be made unique by adding a tiny value (proportional 
                              to the ordinal position of each event) to the data.}
    \item{\code{split_val}:}{A hairy array, aka list of vectors.  Each list element is a
                             vector representing the median values at which the data were split.}
    \item{\code{split_axis}:}{A hairy array, aka list of vectors.  Each list element is a
                             vector representing the axis on which the data were split.}
    \item{\code{binBoundary}:}{An object of class \code{binBoundary}, used to hold boundary information
                               used primarily for visualization.}
	\item{\code{.cRecursions}:}{Private value to hold the number of levels of recursion used to 
	                            construct this model. Using \code{nRecursions} the resolution of a fingerprint
	                            can be reduced, but it can never exceed this value.}
    \item{\code{.tmp_tags}:}{Scratch array, total number of events in the 
                             training set long, that keeps track of the event's
                             bin number.  (this exists only to provide the underlying C function with
                             a persistent scratch space it needs for bookkeeping.  Not useful
                             to the user.)}
  }
}
\section{Methods}{
  \describe{
    \item{show}{shows the contents of the model.}
	 }
}

\note{ 
When creating a model you must keep in mind that it doesn't make sense to create more bins 
(which is \eqn{2^{nRecursions}}{2^nRecursions}) than the total number of events used to create the model. 
The constructor checks for this.

When creating a model, you should specify only parameters that are common to all of the instances 
(\emph{flowFrame}s) in training data. For example, it does not make sense to compare PE from one 
\emph{flowFrame} with FITC from another.
 }

\references{
M. Roederer, et. al. (2001) Probability Binning Comparison: A Metric for Quantitating 
Multivariate Distribution Differences, \emph{Cytometry} \bold{45}, 47-55.

W. Rogers et. al. (2008) Cytometric Fingerprinting: Quantitative Characterization of 
Multivariate Distributions, \emph{Cytometry Part A} \bold{73}, 430-441.
}
\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}


\seealso{
	\link{flowFPModel} - Constructor.
}

\examples{

	# load a flowSet to use as an example.
	library(flowFP)
	data(fs1)
	mod = flowFPModel(fs1)
}

\keyword{classes}
