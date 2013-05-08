\name{Onto2LC}
\alias{Onto2LC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to calculate the Level Coefficients
  %%  ~~function to do ... ~~
}
\description{
  Function to calculate the Level Coefficients for different GO
  ontologies.
}
\usage{
Onto2LC(onto = "BP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{onto}{The GO ontologies, i.e. "BP", "MF" and "CC".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The LC values with GO terms as names. 
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
#BPLC <- Onto2LC("BP")
#MFLC <- Onto2LC("MF")
#CCLC <- Onto2LC("CC")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ LC }

