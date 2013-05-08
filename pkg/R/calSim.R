t2ancestor <- function(t1a){
  ancestor <- NULL
  eval(parse(text=paste("ancestor <- GO",t1a[3],"ANCESTOR", sep="")))
  t2a <- rbind(t1a, cbind(get(t1a[1], ancestor), t1a[2], t1a[3]))
  t2a <- t2a[t2a[,1]!="all",]
  t2a
}

#gid to GO including ancestor
g2go <- function(gid, tid=NULL, filter=NULL){
  if(is.null(tid)){
    terms1 <- get(gid, org.Hs.egGO)
    if(is.null(names(terms1))){
      return(NA)
    }
    terms1 <- do.call("rbind", lapply(terms1, function(x)unlist(x)))
  }else{
    terms1 <- cbind(GOID=tid, Evidence=NA, Ontology=Ontology(tid))
  }
  if(!is.null(filter)){
    terms1 <- terms1[!(terms1[,"Evidence"] %in% filter),,drop=FALSE]
  }
  if(any(duplicated(terms1[,1])==TRUE)){
    d1 <- tapply(terms1[,2], terms1[,1], function(x)paste(sort(x), collapse=","))
    d2 <- terms1[match(names(d1), terms1[,1]),3]
    terms1 <- cbind(names(d1), d1, d2)
  }
  if(nrow(terms1)>1){
    t2 <- do.call("rbind", lapply(1:nrow(terms1), function(x)t2ancestor(terms1[x,])))
  }else{
    t2 <- t2ancestor(terms1)
  }
  e1 <- tapply(t2[,2], t2[,1], function(x)paste(sort(x), collapse=";"))
  w1 <- table(t2[,1])
  o1 <- t2[match(names(e1), t2[,1]),3]
  go <- data.frame(names(e1), e1, o1, as.numeric(w1))
  colnames(go) <- c(colnames(t2), "Count")
  go <- go[order(go[,3], go[,1]),]
  #go <- apply(go, 1, as.list)
  return(go)
}

calSim <- function(g1, g2=NULL, tids=NULL, ontology="ALL", an.go=NULL, term2ancestor=TRUE, count=FALSE){
  #gid to GO
  if(is.null(an.go)){
    t1 <- g2go(g1)
    if(is.null(g2)){
      if(term2ancestor==TRUE){
        t2 <- g2go(tid=tids)
      }else{
        t2 <- cbind(GOID=tids, Evidence=NA, Ontology=Ontology(tids))
      }
    }else{
      t2 <- g2go(g2)
    }
  }else{
    t1 <- mget(g1, an.go, ifnotfound=NA)[[1]]
    if(is.null(g2)){
      if(term2ancestor==TRUE){
        t2 <- g2go(tid=tids)
      }else{
        t2 <- cbind(GOID=tids, Evidence=NA, Ontology=Ontology(tids))
      }
    }else{
      t2 <- mget(g2, an.go, ifnotfound=NA)[[1]]
    }
  }
  if(is.null(nrow(t1)) | is.null(nrow(t2))){
    return(rep(NA, 3))
  }
  #contingency table
  if(ontology=="ALL"){
    terms <- names(as.list(LC))
    t1g <- as.character(t1[,1])
    t2g <- as.character(t2[,1])  
  }else{
    Terms <- NULL
    eval(parse(text=paste("Terms <- GO",ontology,"Term", sep="")))
    terms <- names(as.list(Terms))
    t1g <- as.character(t1[t1[,3]==ontology,1])
    t2g <- as.character(t2[t2[,3]==ontology,1])
  }
  goint <- intersect(t1g, t2g)
  gount <- unique(c(t1g, t2g))
  gdiff1 <- setdiff(t1g, t2g)
  gdiff2 <- setdiff(t2g, t1g)
  con2 <- sum(unlist(mget(setdiff(terms, gount), LC, ifnotfound=NA)), na.rm=TRUE)
  if(count){
    count1 <- pmin(t1[match(goint, t1g), "Count"], t2[match(goint, t2g), "Count"])
    con1 <- sum(unlist(mget(goint, LC, ifnotfound=NA))*count1, na.rm=TRUE)
    count2 <- t1[match(gdiff1, t1g), "Count"]
    diff1 <- sum(unlist(mget(gdiff1, LC, ifnotfound=NA))*count2, na.rm=TRUE)
    count3 <- t2[match(gdiff2, t2g), "Count"]
    diff2 <- sum(unlist(mget(gdiff2, LC, ifnotfound=NA))*count3, na.rm=TRUE)
   }else{
    con1 <- sum(unlist(mget(goint, LC, ifnotfound=NA)), na.rm=TRUE)
    diff1 <- sum(unlist(mget(gdiff1, LC, ifnotfound=NA)), na.rm=TRUE)
    diff2 <- sum(unlist(mget(gdiff2, LC, ifnotfound=NA)), na.rm=TRUE)   
  }
  Kp <- Kappa(matrix(c(con1, diff1, diff2, con2), ncol=2))
  c(Kp$Unweighted, z=Kp$Unweighted[1]/Kp$Unweighted[2])
}
