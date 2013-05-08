\name{ovreGO}
\alias{ovreGO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to analyze over-represented GO terms for a giving gene list
}
\description{
  Given a list of genes, the function can be used to find
  over-represented GO terms. The function is based on topGO package.
}
\usage{
ovreGO(genes, allgenes, ontology = "BP", algorithm = "parentchild",
statistic = "fisher", nterm = 100, plot = FALSE, nword = 30, min.freq =
1, scale = c(2, 0.5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{The gene set to analyze over-represented GO terms.}
  \item{allgenes}{The background gene set compared to the interested
    gene set. It can be all GO annotated genes.}
  \item{ontology}{The GO ontology.}
  \item{algorithm}{The details of algorithm can be find in the
    \code{runTest} function from the topGO package.}
  \item{statistic}{The test method to use. The details can be find in
    the \code{runTest} function.}
  \item{nterm}{The number of top over-represented GO terms.}
  \item{plot}{Wether to plot the returned GO terms. If TRUE, the
    wordcloud will be plotted.}
  \item{nword}{The maximum character length of GO terms in the wordcloud
    plot.}
  \item{min.freq}{words with frequency below min.freq will not be plotted}
  \item{scale}{A vector of length 2 indicating the range of the size of
    the words.}
  \item{...}{More options for \code{wordcloud} function.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame from \code{GenTable} function of topGO package. Two more
  colummn "Score" and "wScore" are added to list the significant scores
  and LC weighted scores.
  If the plot option is TRUE, a wordcloud will be plotted based on the
  over-represented GO terms and its weighted significant scores.
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
#paths <- as.list(KEGGPATHID2EXTID)
#paths <- paths[grep("^hsa", names(paths))]
#allgenes <- unique(unlist(paths))
#ovterms <- ovreGO(genes=paths[[1]], allgenes=allgenes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ over-represent }
\keyword{ wordcloud }
