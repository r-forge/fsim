SearchTerm <- function(fun, fuzzy=FALSE, method="union"){
  keywords <- as.list(keyword)
  if(fuzzy){
    idx <- lapply(fun, function(x)agrep(x, unlist(keywords), ignore.case=TRUE))
  }else{
    idx <- lapply(fun, function(x)grep(x, unlist(keywords), ignore.case=TRUE))
  }
  idx <- Reduce(method, idx)
  terms <- names(keywords)[idx]
  Terms <- data.frame(GOID=terms, Ontology=Ontology(terms), Term=Term(terms), row.names=NULL)
  Terms <- Terms[order(Terms$Ontology, Terms$GOID),]
  return(Terms)
}

#SearchGene.R
SearchGene <- function(gene=NULL, symbol=NULL, terms=NULL, targets="ALL", n=50, an.go, ontology="ALL", term2ancestor=TRUE, map="org.Hs.eg.db", plot=FALSE){
  require(package=map, character.only=TRUE)
  eval(parse(text=paste("mapGO <- ",sub(".db$", "", map),"GO", sep="")))
  eval(parse(text=paste("mapSYMBOL2EG <- ",sub(".db$", "", map),"SYMBOL2EG", sep="")))
  eval(parse(text=paste("mapSYMBOL <- ",sub(".db$", "", map),"SYMBOL", sep="")))
  eval(parse(text=paste("mapGO2EG <- ",sub(".db$", "", map),"GO2EG", sep="")))
  
  if(!is.null(gene)){
    t1 <- names(mget(gene, mapGO, ifnotfound=NA)[[1]])    
  }else if(!is.null(symbol)){
    gene <- mget(symbol, mapSYMBOL2EG, ifnotfound=NA)[[1]]
    t1 <- names(mget(gene, mapGO, ifnotfound=NA)[[1]])    
  }else if(!is.null(terms)){
    t1 <- as.character(terms)
  }
  
  if(any(targets=="ALL")){
    g1 <- lapply(mget(t1, mapGO2EG, ifnotfound=NA), unique)
    sg <- sort(table(unlist(g1)), decreasing=TRUE)
    if(length(sg)>n){
      if(is.null(gene)){
        sg1 <- head(sg, n=n)
      }else{
        sg1 <- head(sg[names(sg)!=gene],n=n)
      }
    }else{
      sg1 <- sg
    }
    sim <- t(sapply(names(sg1), function(x)calSim(g1=x, tids=t1, ontology=ontology, an.go=an.go, term2ancestor=term2ancestor)))

  }else{
    g1 <- lapply(mget(t1, mapGO2EG, ifnotfound=NA), unique)
    sg <- sort(table(unlist(g1)), decreasing=TRUE)
    sg1 <- sg[names(sg) %in% targets]
    sg1 <- sg[match(targets, names(sg))]
    sg1[is.na(sg1)] <- 0
    names(sg1) <- targets
    sim <- t(sapply(names(sg1), function(x)calSim(g1=x, tids=t1, ontology=ontology, an.go=an.go, term2ancestor=term2ancestor)))
  }
  Symbol <- unlist(mget(names(sg1), mapSYMBOL, ifnotfound=NA))
  Sim <- data.frame(Symbol, sharedTerms=sg1, sim)
  Sim <- Sim[order(Sim[,"z.value"], decreasing=TRUE),]

  if(plot==TRUE){
    gs <- rownames(Sim)
    termlist <- mget(gs, mapGO, ifnotfound=NA)
    termlist <- lapply(termlist, function(x)unique(names(x)))
    tf <- do.call("cbind", lapply(termlist, function(x)t1 %in% x))
    #tn <- Term(t1)
    #nword=25
    #tl <- nchar(tn)
    #tn[tl>nword] <- substr(tn[tl>nword], start=1, stop=nword-2)
    #tn[tl>nword] <- paste(tn[tl>nword],"..",sep="")
    rownames(tf) <- t1
    #library(ggplot2)
    #library(reshape2)
    tfm <- melt(tf)
    tfm$Var1 <- factor(tfm$Var1)
    tfm$Var2 <- factor(tfm$Var2, levels=rev(gs))
    #tfm$value <- factor(tfm$value, levels=c("TRUE", "FALSE"))
    tfm$lc <- unlist(mget(as.character(tfm$Var1), LC))
    tfm$lc[tfm$value==FALSE] <- 0
    colnames(tfm) <- c("Term", "Gene", "Annotated", "LC")
    p <- ggplot(data=tfm, aes_string(x="Term", y="Gene")) + geom_tile(aes(fill=LC), colour="white") + scale_fill_gradient(low="grey", high="steelblue")
    print(p)
  }
  return(Sim)
}
