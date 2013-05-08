\name{SearchTerm}
\alias{SearchTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to search keywords for similar GO terms. 
  %%  ~~function to do ... ~~
}
\description{
  Providing a group of functional keywords, the function return a group
  of related GO IDs based on GO terms and their synonym terms.
}
\usage{
SearchTerm(fun, fuzzy = FALSE, method = "union")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{Charactors of a group of function keywords.
%%     ~~Describe \code{fun} here~~
}
  \item{fuzzy}{if TRUE, search for approximate matches of functional
    keywords. if FALSE, only search for matched keywords.}
  \item{method}{The method to combine GO terms from different
    keywords. If it is "union", the function return all matched
    terms. if it is "intersect", only overlapped functions are
    returned.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data.frame object including matched GOID, Ontology and Term.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
SearchTerm("chromatin remodel")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~keyword }
\keyword{ ~search }% __ONLY ONE__ keyword per line
