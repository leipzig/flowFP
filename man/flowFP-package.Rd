\name{flowFP-package}
\docType{package}
\alias{flowFP-package}
\title{
Package overview

}
\description{
This package is used to explore flow cytometry data through the use of 
\emph{fingerprints}. The broad aim of the package is to transform flow 
cytometric data into a form amenable to algorithmic analysis tools. Thus, it is
useful to think of \pkg{flowFP} as an intermediate step between the acquisition
of high-throughput flow cytometric data and empirical modeling, machine learning
and knowledge discovery.

A \emph{fingerprint} is a feature vector meant to efficiently represent the multivariate
probability distribution function for a flow cytometry data set.  It is produced
by first creating a data-relevant \emph{model} of a space, and then applying the
model to a dataset, thereby producing fingerprints.  
Model creation is done through the \link{flowFPModel} constructor which can be 
customized via function arguments. After the model is built, it can be applied
to arbitrary \link[flowCore:flowFrame-class]{flowFrame}s or 
\link[flowCore:flowSet-class]{flowSet}s using the \link{flowFP} 
constructor. The resulting S4 object implements plotting and summary methods 
that allow the user to compare and contrast instances, using the \code{flowFPModel}
as a sort of basis representation, akin for example to trigonometric functions in a 
Fourier Transform.

This package is closely integrated with \code{flowCore}.
You will want to become familiar with it in order to effectively use \pkg{flowFP}.
}

\details{
\tabular{ll}{
Package: \tab flowFP\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Depends: \tab R(>= 2.5.0), flowCore, flowViz\cr
Collate: \tab flowFPModel.R flowFP.R\cr
Bioinformatics: \tab Flowcytometry, CellBasedAssays, Clustering, Statistics, Visualization\cr
License: \tab Artistic-2.0\cr
Built: \tab R 2.8.0; unix\cr
}
}

\section{Classes}{
\code{\link{flowFPModel-class}} is the fundamental class for the 
\pkg{flowFP} package.  It represents the multivariate probability distribution 
function for a flow cytometry data set. Information is maintained in a number of
slots, which should only be accessed through methods, described below, \bold{not
by direct use of the @ operator}. For a complete detailed list of slot names and
descriptions look at the \code{\link{flowFPModel-class}} help page.


\code{\link{flowFP-class}} extends the \code{flowFPModel} and contains
additional slots to record the assignment to and number of events in the bins
of a \code{flowFPModel}. Methods are supplied to retrieve and visualize the 
the contents of a \code{flowFP}.


\code{\link{flowFPPlex-class}} is a container for a set of congruent \code{flowFP}
objects (by congruent, we mean that each \code{flowFP} is a description of the same
set of instances). When constructing or appending \code{flowFP}s into a plex,
simple error checking is done to ensure each instance in each \code{flowFP} is
equivalent. Both the \code{sampleNames} and the \code{sampleClasses} slots are consulted
for internal consistency. An error is generated if any of the \code{flowFP}s 
to be joined in a plex contain different \code{sampleNames} and/or s\code{ampleClasses}.
}

\section{Constructors}{

\code{\link{flowFPModel}} is the constructor for this class.
	\preformatted{
flowFPModel(fcs, name="Default Model", parameters=NULL, nRecursions="auto",
            dequantize=TRUE, sampleSize=NULL)
}
	\tabular{ll}{
      \code{fcs} \tab Either an \code{flowFrame} or a \code{flowSet} used to create the model.\cr
	  \code{name} \tab Name given the model.\cr
	  \code{parameters} \tab Parameters to consider when constructing the model.(e.g. c(1,5)). \cr
   	  \code{nRecursions} \tab The number of level of recursive subdivision. \cr
   	  \code{dequantize} \tab Setting this value causes a small incremental value to be added\cr
   	                    \tab to each event starting with 1e-8. This effectively reduces the\cr
   	                    \tab number of duplicate values to break ties when binning.\cr
   	  \code{sampleSize} \tab The max number of events to use out of each fcs file in a\cr
   	                    \tab \code{flowFrame} or \code{flowSet}.\cr
	  }


\code{\link{flowFP}} is the constructor for this class.

	\preformatted{
flowFP(fcs, model=NULL, sampleClasses=NULL, ...)
}
	\tabular{ll}{
	   \code{fcs} \tab Either an \code{flowFrame} or a \code{flowSet} used to create the fingerprint.\cr
	   \code{model} \tab A model created using \code{flowFPModel} or NULL in which case\cr
	                \tab  a model will be created from the fcs data supplied.\cr
	   \code{sampleClasses} \tab List of sample class names to be assigned in order to the instances. \cr
   	   \code{...} \tab If the \code{model} is NULL the parameters are passed through to\cr
   	              \tab \code{flowFPModel}.
	  }


\code{\link{flowFPPlex}} is the constructor for this class.

	\preformatted{
flowFPPlex(fingerprints=NULL)
}
	\tabular{ll}{
	   \code{fingerprints} \tab Either an single \code{flowFP} or a list of \code{flowFP}s\cr
	                       \tab (e.g. \code{fingerprints = c(fp1, fp2)}). All of the fingerprints \cr
	                       \tab must share the same sample names, and class names (or no class names).
	  }
}


\note{For further information please see the vignette.}

\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}

\seealso{ 
flowCore, flowViz
}

\examples{

	# load a flowSet to use as an example.
	library(flowFP)
}
\references{
M. Roederer, et. al. (2001) Probability Binning Comparison: A Metric for Quantitating 
Multivariate Distribution Differences, \emph{Cytometry} \bold{45}, 47-55.

W. Rogers et. al. (2008) Cytometric Fingerprinting: Quantitative Characterization of 
Multivariate Distributions, \emph{Cytometry Part A} \bold{73}, 430-441.
}

\keyword{package}
