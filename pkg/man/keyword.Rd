\name{keyword}
\alias{keyword}
\docType{data}
\title{keyword}
\description{An R object provides mappings between GO IDs and GO terms
  including their synonyms.}
\usage{data(keyword)}
\format{An environment}
\source{
## script to produce keyword.
# goterms <- Term(GOTERM)
# syn <- unlist(lapply(Synonym(GOTERM), function(x)paste(x, collapse="; ")))
# keywords <- paste(goterms, syn, sep="; ")
# names(keywords) <- goterms
# keyword <- new.env()
# for(i in 1:length(goterms)){
#   assign(names(goterms)[i], keywords[i], envir=keyword)
# }
}


