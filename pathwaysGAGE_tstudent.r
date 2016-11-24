library(biomaRt)
library(sqldf)
 require(gage)
 data(kegg.gs)


# # # # # # # # # # # # # #pathways LL-TT # # # #    
data_expr = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genes.htseq.normalized.txt")

top=data_expr[,c("X12NT", "X9NT","X10L100","X7L100")]
top$ensembl_gene_id=data_expr$Gene_ID
 
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 

eattr<-c("ensembl_gene_id","entrezgene" )
bpres<-getBM(attributes=eattr,filters="ensembl_gene_id",values=top$ensembl_gene_id,mart= mart)
results1  <-bpres[bpres$entrezgene!="",]
   
res1  <- sqldf("select *
              from results1 
              group by ensembl_gene_id")

topbis <- sqldf("SELECT *
FROM top
LEFT  JOIN res1 
ON top.ensembl_gene_id=res1.ensembl_gene_id
WHERE res1.ensembl_gene_id!=\"NA\" ;")
  
topbis2<-as.matrix(topbis[,c("X12NT", "X9NT","X10L100","X7L100")] )
row.names(topbis2)<- topbis$entrezgene 
   gage_Treated_NonTreated <- gage( as.matrix(topbis2  ),
   gsets = kegg.gs, 
      ref = grep ("NT",colnames(topbis2  )), 
     samp = grep ("L100",colnames(topbis2 )),
    compare='paired') #   ,test4up=TRUE)
# rien de significatif










