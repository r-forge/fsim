### R code from vignette source 'FSim.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(FSim)


###################################################
### code chunk number 2: calSim1
###################################################
calSim(g1="9", g2="10", ontology="BP", an.go=an.Hs.egGO)
calSim(g1="9", g2="10", ontology="ALL", an.go=an.Hs.egGO)


###################################################
### code chunk number 3: calSim2
###################################################
terms <- names(get("10", org.Hs.egGO))
calSim(g1="9", tids=terms, an.go=an.Hs.egGO)


###################################################
### code chunk number 4: SearchGene
###################################################
SearchGene(symbol="NAT2", an.go=an.Hs.egGO, targets="ALL", n=10)


###################################################
### code chunk number 5: SearchGeneSet
###################################################
SearchGene(gene="9", targets=c("10", "100", "124"), 
           an.go=an.Hs.egGO)


###################################################
### code chunk number 6: SearchTerm
###################################################
t1 <- names(get("9", org.Hs.egGO))
t1
SearchGene(terms=t1, an.go=an.Hs.egGO, n=5)


###################################################
### code chunk number 7: SearchKeywords
###################################################
t2 <- SearchTerm(fun=c("chromatin remodeling", "histone binding"))
t2


###################################################
### code chunk number 8: SearchKeywords
###################################################
SearchGene(terms=t2$GOID[c(1,3,6)], an.go=an.Hs.egGO, n=5)


###################################################
### code chunk number 9: SearchGeneSet
###################################################
library(KEGG.db)
geneset <- get("hsa00232", KEGGPATHID2EXTID)
geneset


###################################################
### code chunk number 10: SearchGeneSet
###################################################
paths <- as.list(KEGGPATHID2EXTID)
paths <- paths[grep("^hsa", names(paths))]
allgenes <- unique(unlist(paths))
BPterms <- ovreGO(genes=geneset, allgenes=allgenes, 
                  ontology="BP", nterm=10)
BPterms


###################################################
### code chunk number 11: SearchGeneSet
###################################################
SearchGene(terms=BPterms$GO.ID, targets="ALL", 
           an.go=an.Hs.egGO, n=5)


###################################################
### code chunk number 12: evaluation
###################################################
MFterms <- ovreGO(genes=geneset, allgenes=allgenes, 
                  ontology="MF", nterm=10)
CCterms <- ovreGO(genes=geneset, allgenes=allgenes, 
                  ontology="CC", nterm=10)
allterms <- c(BPterms$GO.ID, MFterms$GO.ID, CCterms$GO.ID)


###################################################
### code chunk number 13: evaluation
###################################################
score1 <- SearchGene(terms=allterms, targets=geneset, 
                     an.go=an.Hs.egGO, ontology="ALL", 
                     term2ancestor=FALSE)
score1


###################################################
### code chunk number 14: evaluation
###################################################
set.seed(1)
ctlgene <- sample(setdiff(allgenes, geneset), 50)
score2 <- SearchGene(terms=allterms, targets=ctlgene, 
                     an.go=an.Hs.egGO, ontology="ALL", 
                     term2ancestor=FALSE)
head(score2)
pv <- suppressWarnings(ks.test(score1$z.value, score2$z.value))
pv


###################################################
### code chunk number 15: evaluation1
###################################################
boxplot(score1$z.value, score2$z.value, 
        names=c("score1", "score2"), 
        main=paste("KS test: ", format(pv$p.value, digits=4)))


###################################################
### code chunk number 16: heatmap
###################################################
res1 <- SearchGene(symbol="NAT2", an.go=an.Hs.egGO, 
                   targets="ALL", n=10, plot=TRUE)


###################################################
### code chunk number 17: wordcloud
###################################################
res2 <- ovreGO(genes=geneset, allgenes=allgenes, 
               ontology="BP", plot=TRUE, scale=c(1,0.5))


###################################################
### code chunk number 18: session
###################################################
sessionInfo()


