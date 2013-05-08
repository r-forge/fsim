\name{calSim}
\alias{calSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to calculate similarity score
  %%  ~~function to do ... ~~
}
\description{
  The function can be used to calculate similarity scores between two
  genes or gene and GO terms.
}
\usage{
calSim(g1, g2 = NULL, tids = NULL, ontology = "ALL", an.go = NULL,
term2ancestor = TRUE, count = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g1}{The first gene ID to compare.}
  \item{g2}{The second gene ID to compare.}
  \item{tids}{The GO term IDs to compare if the second gene ID is NULL.}
  \item{ontology}{The GO ontology. It can be "BP", "MF", "CC" and
    "ALL".}
  \item{an.go}{The R enviroment that provide mappings between gene IDs
    and GO IDs including their ancestors.}
  \item{term2ancestor}{The option works when terms is not NULL. If TRUE,
    the input term and their ancestor terms will be used. If FALSE, the
    input term will be used directly.}
  \item{count}{Some ancestor GO terms can be counted serveral times for
    the same gene. If count is TRUE, the count number will be a weight
    in the calculation. If count is FALSE, the counts will not be used.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The returned values include Kappa value, Approximate Standard Error
  (ASE) and Z value.
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
  \code{\link{SearchGene}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
g1 <- "9"
g2 <- "10"
calSim(g1, g2, an.go=an.Hs.egGO, ontology="ALL", count=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kappa }
