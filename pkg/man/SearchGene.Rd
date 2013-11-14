\name{SearchGene}
\alias{SearchGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to search for functionally similar genes 
  %%  ~~function to do ... ~~
}
\description{
  For interested gene or a group of GO terms, the function can be used
  to search for functionally similar genes in a target gene list or all
  GO annotated genes. 
}
\usage{
SearchGene(gene = NULL, symbol = NULL, terms = NULL, targets = "ALL", n
= 50, an.go, ontology = "ALL", term2ancestor = TRUE, map =
"org.Hs.eg.db", plot = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{The input gene id for searching functionally similar gene
    list.}
  \item{symbol}{The input gene symbol for searching gene list.}
  \item{terms}{The input GO terms for searching gene list.}
  \item{targets}{The searched gene list. It can be a group of gene IDs
    or "ALL". "ALL" means to search for the most related gene list in
    all GO annotated genes.}
  \item{n}{If targets is "ALL", n is the amount of target genes based on
    the number of shared GO terms between input gene and target genes in
    decreasing order.}
  \item{an.go}{The R enviroment that provide mappings between gene IDs
    and GO IDs including their ancestors.}
  \item{ontology}{The GO ontology to use in the search. It can be "BP",
    "MF", "CC" and "ALL" for the three ontologies.}
  \item{term2ancestor}{The option works when terms is not NULL. If TRUE,
    the input term and their ancestor terms will be used. If FALSE, the
    input term will be used directly.}
  \item{map}{The GO annotation package.}
  \item{plot}{If TRUE, a heatmap will be plotted to show the shared GO
    terms.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  
  A data.frame object including functionally similar genes. The columns
  are
  \item{Symbol}{The gene symbol of similar genes.}
  \item{sharedTerms}{The number of shared GO terms between input gene
    and similar genes.}
  \item{value}{The kappa value between input gene and returned gene.}
  \item{ASE}{The Approximate Standard Error for the kappa value.}
  \item{z.value}{The Z score to standardize the kappa value.}
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
  \code{\link{SearchTerm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# To search similar genes in a gene set.
# SearchGene(gene="10", targets=c("9", "2950"), an.go=an.Hs.egGO)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ search }
\keyword{ gene }% __ONLY ONE__ keyword per line
