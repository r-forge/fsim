ovreGO <- function(genes, allgenes, ontology="BP", algorithm="parentchild", statistic="fisher", map="org.Hs.eg.db", nterm=100, plot=FALSE, nword=30, min.freq=1, scale=c(2,0.5), ...){
  geneList <- factor(as.integer(allgenes %in% genes))
  names(geneList) <- allgenes
  godata <- new("topGOdata", ontology=ontology, allGenes=geneList, annot = annFUN.org, mapping=map)
  result <- runTest(godata, algorithm = algorithm, statistic = statistic)
  orGO <- GenTable(godata, result, topNodes=nterm)
  lc <- unlist(mget(orGO$GO.ID, LC, ifnotfound=NA))
  pv <- score(result)
  pv <- pv[match(orGO$GO.ID, names(pv))]
  #lc <- unlist(mget(orGO$GO.ID, LC))
  Score <- -log10(pv)
  wScore <- Score * lc
  wScore[is.na(wScore)] <- 0
  orGO <- cbind(orGO, Score, wScore)
  if(plot==TRUE){
    words <- Term(as.character(orGO$GO.ID))
    wn <- nchar(words)
    words[wn>nword] <- substr(words[wn>nword], start=1, stop=nword-2)
    words[wn>nword] <- paste(words[wn>nword],"..",sep="")
    colors <- brewer.pal(8,"Dark2")
    wordcloud(words, freq=round(orGO[,"wScore"]), min.freq=min.freq, scale=scale, colors=colors, ...)
  }
  orGO
}
