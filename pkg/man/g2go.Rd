\name{g2go}
\alias{g2go}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to extract GO annotation and corresponding ancestor terms
}
\description{
  The function is used to extract annotated GO terms including their
  ancestor terms based on mapping information.
}
\usage{
g2go(gid, tid = NULL, filter = NULL, map = org.Hs.egGO)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gid}{The gene ID to extract GO annotation.}
  \item{tid}{If \code{tid} is a group of GO IDs instead of NULL, the
    ancestors of input GO IDs will return.} 
  \item{filter}{The evidence code to filter GO terms.}
  \item{map}{The mapping between gene IDs to go IDs.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame including GO IDs, Evidence codes, Ontology and their
  counts.
}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Qiang Hu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
g2go("10")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
