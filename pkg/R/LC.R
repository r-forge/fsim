#buildGOgraph from ontoTools
buildGOgraph <- function (useenv = GOMFPARENTS){
  nds <- unique(ls(useenv))
  eds <- mget(nds, useenv, ifnotfound = NA)
  eds <- lapply(eds, function(x) {
    lk <- match(x, nds)
    if (length(lk) == 0) 
      return(list(edges = character(0)))
    else if (length(lk) == 1 && is.na(lk)) 
      return(list(edges = character(0)))
    else if (any(is.na(lk))) 
      return(list(edges = lk[!is.na(lk)]))
    return(list(edges = lk))
  })
  tmp <- new("graphNEL", nodes = nds, edgeL = eds, edgemode = "directed")
  attr(tmp, "toolInfo") <- library(help = GO.db)$info[[1]][c(1, 4)]
  updateGraph(tmp)
}

Onto2LC <- function(onto="BP"){
  eval(parse(text=paste("graph <- buildGOgraph(GO",onto,"PARENTS)", sep="")))
  levels <- buildLevels(graph)
  treetop <- leaves(graph, degree.dir="in")
  M <- levels$noOfLevels
#define LC(treetop)=1
  treetopLC <- rep(1, length(treetop))
  names(treetopLC) <- treetop
#the children of terms are in different level:
  LC <- treetopLC
  for(l in (M-1):1){ #15-1
    tl <- mget(as.character(l), levels$level2nodes)[[1]]
    tt <- tl %in% names(LC)
    ttl <- tl[!tt]
    ttLC <- NULL
    CHILDREN <- NULL
    for(i in 1:length(ttl)){
      eval(parse(text=paste("CHILDREN <- GO",onto,"CHILDREN", sep="")))
      tc <- mget(ttl[i], CHILDREN)[[1]]
      tcLC <- LC[match(tc, names(LC))]
    #ttLC[i] <- (l/(l+1))*mean(tcLC)
      tcLCc <- c()
      for(j in 1:length(tc)){
        tcLCc[j] <- tcLC[j]*(l/(get(tc[j], levels$nodes2level)))
      }
      ttLC[i] <- mean(tcLCc)
    }
    names(ttLC) <- ttl
    LC <- c(LC, ttLC)
  }
  return(LC)
}
