library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",
  host="sep2011.archive.ensembl.org")

genesRnaSeq = list(
  over = scan("over.txt",what=character()), ##1193 ensembl_id
  under = scan("under.txt",what=character()) ##75 ensembl_id
  ) 

#### Ici je fais l'analyse GO pour trouver les termes enrichis:
library(topGO)
genome = getBM(attr=c("ensembl_gene_id","go_id","name_1006"),
  filter=c("with_go","biotype"),
  values=list(TRUE,"protein_coding"),mart=ensembl)
allGenes = unique(genome[,"ensembl_gene_id"])
gene2GO = split(genome[,"go_id"],genome[,"ensembl_gene_id"])
tab = list()
goterms = list()
gonames = list()
for (set in names(genesRnaSeq)) {
    geneList = factor(as.integer(allGenes %in% genesRnaSeq[[set]]))
    names(geneList) = allGenes
    data = new("topGOdata",description=set, ontology="BP",
      allGenes=geneList,annot=annFUN.gene2GO,gene2GO=gene2GO)
    result = list(classicFisher=runTest(data, statistic="fisher", algo="classic"),
      elimFisher=runTest(data, statistic="fisher", algo="elim"))
    tab[[set]] = GenTable(data, classicFisher=result$classicFisher, 
      elimFisher=result$elimFisher, 
      orderBy="elimFisher", ranksOf="elimFisher", topNodes=5)
    for (go in tab[[set]][,"GO.ID"]) {#tab[[set]]$elimFisher<.1 # important point
        goterms[[go]] = genesInTerm(data,go)[[go]]
        gonames[[go]] = tab[[set]]$Term[which(tab[[set]]$GO.ID == go)]
    }
#    showSigOfNodes(data,score(result$classicFisher),first=10,useInfo="all")
}
########## Maintenant les ratios pour le barplot 
full_table = data.frame()
for (go in names(goterms)) {
    all_atGo = goterms[[go]]
    full_table["genome",go] = length(all_atGo)/length(allGenes)
    full_table["underexpressed",go] = sum(genesRnaSeq[["under"]] %in% all_atGo)/length(genesRnaSeq[["under"]])
    full_table["overexpressed",go] = sum(genesRnaSeq[["over"]] %in% all_atGo)/length(genesRnaSeq[["over"]])
}
par(mfrow=c(2,2))
for (g in names(full_table)) {
    main = gonames[[g]]
    barplot(full_table[,g],names=row.names(full_table),main=main,las=2)
}

options(stringsAsFactors=F)       
exprdata = read.delim("genes_differential_Hes1_Cre-Hes1_ctrl.txt",skip=1)
fullexpr = data.frame(log2ratio=exprdata$log2FoldChange,row.names=sapply(strsplit(exprdata$id,split="|",fixed=T),"[",1))
for (go in names(goterms)) {
    term = gonames[[go]]
    genes = goterms[[go]]
    x1 = fullexpr$log2ratio
    x2 = fullexpr[genes,"log2ratio"]
    plot(density(x1,na.rm=T),ylim=c(0,.5),main=term,xlab="log2 expression")
#    hist(x2,add=T,col='red',border='red',freq=F,br=100)
    lines(density(x2,na.rm=T),col='red')
    pv = wilcox.test(x1,x2,alternative="less")$p.value
    text(0.1,0.4,paste("p=",round(pv,6),sep=''))
}

