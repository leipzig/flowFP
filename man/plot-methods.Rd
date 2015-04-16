\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,flowFP,flowFrame-method}
\alias{plot,flowFP,flowSet-method}
\alias{plot,flowFP,missing-method}
\alias{plot,flowFPModel,flowFrame-method}
\alias{plot,flowFPModel,flowSet-method}
\alias{plot,flowFPModel,missing-method}
\alias{plot,flowFPPlex,missing-method}

\title{ Methods for visualizing flowFP objects}
\description{
 These methods allow the user to plot flowFP objects with a number of options.
}

\section{Methods}{
\code{plot (x, y, ...)}
\describe{

\item{\code{x = "flowFPModel", y = "missing"}}{ Visualize a \code{flowFPModel}.\cr
                     Optional Args: \code{(parameters=NULL, alpha=1, border="gray", 
                     showbins=1:nFeatures(x), ylim=NULL, xlim=NULL, main="Model", ...)}}

\item{\code{x = "flowFPModel", y = "flowFrame"}}{ Visualize a \code{flowFPModel} along with a \code{flowFrame}.\cr
                     Optional Args: \code{(parameters=NULL, alpha=1, border="gray", 
                     showbins=1:nFeatures(x), ylim=NULL, xlim=NULL, main="Model", ...)}}

\item{\code{x = "flowFPModel", y = "flowSet"}}{ Visualize a \code{flowFPModel} along with a \code{flowSet}.\cr
                     Optional Args: \code{(parameters=NULL, alpha=1, border="gray", 
                     showbins=1:nFeatures(x), ylim=NULL, xlim=NULL, main="Model", ...)}}

\item{\code{x = "flowFP", y = "missing"}}{ Visualize Fingerprints.\cr
                     Optional Args: \code{(type=c("tangle", "stack", "grid", "qc", "plate"), ...) }\cr
                     See Notes for additional details.}

\item{\code{x = "flowFP", y = "flowFrame"}}{ Visualize a single fingerprint with a \emph{flowFrame}.  See Notes. }

\item{\code{x = "flowFP", y = "flowSet"}}{ Visualize Fingerprints with a \emph{flowSet}.\cr
                     Optional Args: \code{(x, y, showbins=NULL, showfp=TRUE, ...)}}

\item{\code{x = "flowFPPlex", y = "missing"}}{ Visualize Fingerprints in a \emph{flowFPPlex}.\cr
                     Optional Args: \code{(type=c("tangle", "stack", "grid","qc", "plate"), ...)} }

}}

\section{Arguments}{
\describe{
\item{\code{parameters}}{Flow parameters involved in fingerprinting.}
\item{\code{showbins}}{Which bins to render.  May be an integer, list of integers, or boolean.  Default: 1:nFeatures(x).}
\item{\code{showfp}}{Boolean.  Do we want to see the fingerprints (\emph{TRUE}) or just the parameter
                     plot (\emph{FALSE}).}
\item{\code{xlim}}{Range of values for the x-axis.  Default: range of the data.}
\item{\code{ylim}}{Range of values for the y-axis.  Default: range of the data.}
\item{\code{type}}{Specify the type of plot.  Current values supported are: \dQuote{tangle}, 
\dQuote{stack}, \dQuote{grid}, \dQuote{qc}, \dQuote{plate}.}
\item{\code{transformation}}{Specify transformation.  Values may be: \dQuote{raw}, 
                             \dQuote{normalized}, \dQuote{log2norm}}
\item{\code{useClasses}}{Boolean.  Colors lines by sampleClass.  Default=\emph{FALSE}.}
\item{\code{linecols}}{Specify color(s) for fingerprint lines.  Default = \sQuote{black}.  When
used in conjunction with \code{useClasses}, this can be a list of colors whose length matches
the number of levels in the \code{sampleClasses} factor object.}
\item{\code{alpha}}{Opacity of rectangles representing bins.}
\item{\code{border}}{Color to draw the border of a bin.  Default: `gray'.}
\item{\code{highlight}}{Which fingerprints to highlight.  Default=\code{1:length(x)}.}
\item{\code{main}}{Title for the plot.  Default: various.}
\item{\code{xlab}}{Label for the x-axis.}
\item{\code{vert_scale}}{Sets the vertical scale for some plots.}
\item{\code{respect}}{Boolean.  If \emph{TRUE}, respect a square aspect ratio.}
\item{\code{method}}{Determine QC metric for fingerprint deviation.  Either \code{sd} (standard deviation) 
                     or \code{max} (maximum deviation).}
\item{\code{red_limit}}{Value corresponding to a red color for the QC metric.}
\item{\code{...}}{Other parameters passed through to generic \code{plot}.}


}}

\section{Notes}{
\describe{
\item{\emph{In conjunction with showbins:}}{generic plot args such as \code{pch} and \code{cex} can be 
used to visualize events in specific bins.}
\item{\emph{For type=\dQuote{tangle} methods:}}{
          Optional Args: \code{(transformation=c("raw", "normalized", "log2norm"), 
                        linecols=NULL, highlight=NULL, ylim=NULL, 
                        useClasses=FALSE, main="Fingerprints", xlab='Feature Index', ...)}}
\item{\emph{For type=\dQuote{stack} methods:}}{
          Optional Args: \code{(transformation=c("raw", "normalized", "log2norm"), 
                        linecols=NULL, useClasses=FALSE, vert_scale=3,
                        ylim=NULL, main="Fingerprints",...)}}
\item{\emph{For type=\dQuote{grid} methods:}}{
          Optional Args: \code{(vert_scale=3, main="Fingerprints", linecols="black",
                        transformation=c("raw", "normalized", "log2norm"),
                        respect=FALSE)}}
\item{\emph{For type=\dQuote{qc} methods:}}{
          Optional Args: \code{(main="Fingerprint Deviation Plot",
                        transformation=c("log2norm", "raw", "normalized"),
                        vert_scale=3, method=c("sd", "max"), red_limit=1.0,
                        respect=FALSE)}}
\item{\emph{For type=\dQuote{plate} methods:}}{
          Optional Args: \code{(main="Fingerprint Deviation Plot",
                        transformation=c("log2norm", "raw", "normalized"),
                        vert_scale=3, method=c("sd", "max"), red_limit=1.0)}}

}}
\keyword{methods}
