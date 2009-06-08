\name{flowFPPlex}
\alias{flowFPPlex}

\title{Fingerprint collection constructor.}
\description{
This constructor creates a flowFPPlex, which is a collection of \link{flowFP}s. 
}
\usage{
flowFPPlex(fingerprints = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fingerprints}{ List of \link{flowFP}s, or \emph{NULL}.  If \emph{NULL} an
                       empty \emph{flowFPPlex} will be prepared, to which \emph{flowFP}s
                       may be added with \link{append-methods}. }
}

\details{
A \emph{flowFPPlex} is a container object for a collection of \link{flowFP}s.  It is useful
in several ways.  First, multiple-tube panels are often used to assay more reagents than
can be done on a particular instrument.  Second, it is sometimes useful to represent the
same data with multiple \link{flowFPModel}s.

Suppose that we have collected data from 1,000 patients, using an 8-tube panel.  Imagine
that Tubes 1 and 8 are isotype and viability tubes, respectively (we will ignore these
tubes for now).  The remaining Tubes 2-7 are of interest from a fingerprinting perspective.
We wish to treat them as a unit.  We might then create models that represent parameters
in each of the tubes, across some (or all) of the patient samples (say, mod2, mod3, ..., mod7).
We could then create corresponding \emph{flowFP}s (say, fp2, fp3, ..., fp7).  These can now
be combined and treated as a single object of type \emph{flowFPPlex}, as:

\code{> plex <- flowFPPlex(c(fp2, fp3, fp4, fp5, fp6, fp7))}\cr\cr
or if you prefer,\cr\cr
\code{> plex <- flowFPPlex (fingerprints=NULL)}\cr
\code{> plex <- append(plex, c(fp2, fp3, fp4, fp5, fp6, fp7))} 

The counts or density matrices can then be extracted simply using methods provided in 
\link{flowFPPlex-class}.

The second idea is to use multiple \emph{models} to represent the same data.  In this case we might
create a model from, say, the \emph{"Normal"} instances (call it mod\_norm), and another model
from the \emph{"Cancer"} instances (mod\_cancer).  We might wish to do this to enhance the detection
of regions of the distribution that are characteristically dominated by one type or the other.
If our \emph{flowSet} of all instances is called "fs1", then our two representations would be:

\code{> fp_norm <- flowFP (fs1, mod_norm)}\cr
\code{> fp_cancer <- flowFP (fs2, mod_cancer)}

and the plex is

\code{plex <- flowFPPlex (c(fp_norm, fp_cancer))}


}
\value{
  Returns a flowFPPlex.
  
}
\author{
Herb Holyst <\email{holyst@mail.med.upenn.edu}>, Wade Rogers <\email{rogersw@mail.med.upenn.edu}>
}


\seealso{ 
\link{flowFPPlex-class}, \link{flowFP}
}
\examples{
data(fs1)
data(fs2)
mod1 <- flowFPModel (fs1, parameters=c(2,5))
mod2 <- flowFPModel (fs2, parameters=c(2,5))
fp1_1 <- flowFP (fs1, mod1)
fp1_2 <- flowFP (fs1, mod2)

plex <- flowFPPlex(c(fp1_1, fp1_2))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
