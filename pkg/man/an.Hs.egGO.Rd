\name{an.Hs.egGO}
\alias{an.Hs.egGO}
\docType{data}
\title{an.Hs.egGO}
\description{An environment provides mappings between human entrez gene
  identifiers and the GO identifiers including their ancestors that they
  are associate with.
}
\usage{data(an.Hs.egGO)}
\format{An environment}
\source{The object can be produced by the function g2go. For example,
  g2go("10")
}
\examples{  
# #all genes
# hsgo <- as.list(org.Hs.egGO)
# gids <- names(hsgo) #42106

##environment
#gene to GO
# an.Hs.egGO <- new.env()
# for(i in 1:length(gids)){
#   assign(gids[i], g2go(gids[i]), envir=an.Hs.egGO)
# }
}

