 library(sqldf)
library(SPIA)
library(biomaRt)
library(limma)
# library(hgu133plus2.db) 
setwd("/Users/benmoham/Desktop/InessBenmessaoud/")
data1 = read.delim("genelist_Treated_3d-Ctrlneg_3d_21082015.txt" )

data1 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genelist_Treated_3d-Ctrlneg_3d.txt" )
data2=data1[grep ("TRUE",data1$P.Value<0.05),]
data1=data2
# names(data1) =
#   "ID"          "logFC"       "P.Value"     "adj.P.Val"   "Fold_Change"
# head(data1)
#                ID     logFC  P.Value  adj.P.Val Fold_Change
# 1 ENSG00000108556  4.282201 2.89e-06 0.01928665   19.456784
# 2 ENSG00000150551 -2.156129 6.52e-06 0.01928665   -4.457173
# 3 ENSG00000157570 -2.185073 1.04e-05 0.01928665   -4.547498
# 4 ENSG00000185101  3.067245 7.33e-06 0.01928665    8.381710
# 5 ENSG00000214021  2.111385 4.81e-06 0.01928665    4.321060
# 6 ENSG00000235703  2.330643 8.15e-06 0.01928665    5.030293

# # # # # # # # # # # # # #pathways (LL+TT)/2-(HEC) # # # #    

top_Treated_NonTreated <- data1
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 
# 32          hsapiens_gene_ensembl              Homo sapiens genes (GRCh38.p5)

eattr<-c("ensembl_gene_id","entrezgene" )
bpres<-getBM(attributes=eattr,filters="ensembl_gene_id",values= top_Treated_NonTreated$ID,mart= mart)
results1_Treated_NonTreated <-bpres[bpres$entrezgene!="",]
   
res1_top_pLNcells_Hepatocytes <- sqldf("select *
              from top_Treated_NonTreated 
              LEFT OUTER JOIN   results1_Treated_NonTreated
              ON top_Treated_NonTreated.ID=results1_Treated_NonTreated.ensembl_gene_id;")
res1_top_pLNcells_Hepatocytes <-res1_top_pLNcells_Hepatocytes [grep("TRUE",res1_top_pLNcells_Hepatocytes$entrezgene!="NA"),]
res1_top_pLNcells_Hepatocytes = res1_top_pLNcells_Hepatocytes[!duplicated(res1_top_pLNcells_Hepatocytes$entrezgene),]
#remove duplicates
# dim(res1_top_pLNcells_Hepatocytes)
# [1] 219   7
  setwd("/Users/benmoham/Desktop/InessBenmessaoud/")
# _UP
sig_genes_diff_Treated_NonTreated <-  res1_top_pLNcells_Hepatocytes$logFC 
names(sig_genes_diff_Treated_NonTreated) <-  res1_top_pLNcells_Hepatocytes$entrezgene 
# length(sig_genes_UP_Treated_NonTreated)
# [1] 177
length(names(sig_genes_diff_Treated_NonTreated))
# 219
length(unique(names(sig_genes_diff_Treated_NonTreated)))
# [1] 218

# length((all_genes_UP_Treated_NonTreated))
# [1] 24553
# > length(unique(all_genes_UP_Treated_NonTreated))
# [1] 24553
 sig_genes_diff_Treated_NonTreated <-  as.vector(res1_top_pLNcells_Hepatocytes$logFC )
names(sig_genes_diff_Treated_NonTreated) <- as.vector(as.character (res1_top_pLNcells_Hepatocytes$entrezgene ) )
# sig_genes_diff_Treated_NonTreated

all.entrezgene <- unique( getBM(attributes = "entrezgene",
                    values = "*", 
                    mart = mart) )
# all_genes_diff_Treated_NonTreated   <-   unique(all.entrezgene[,1]) 

# run SPIA.
 setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffSPIApval005")
 # setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2M2/SPIA")
spia_result_diff_Treated_NonTreated <- spia(de= sig_genes_diff_Treated_NonTreated , 
  all= as.vector (as.character (all.entrezgene[,1])), 
  organism="hsa", plots=FALSE ,nB=200,combine="fisher" ) 
data(kegg.gs)

unique(unlist(kegg.gs[spia_result_diff_Treated_NonTreated ]))
# norminv

# mydir=system.file("extdata/keggxml/hsa",package="SPIA")
# dir(mydir)
# makeSPIAdata(kgml.path=mydir,organism="hsa",out.path="./")


# de=as.vector(as.character (res1_top_pLNcells_Hepatocytes$entrezgene ) )
# require(ReactomePA)
# x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T, organism = "human")
# head(summary(x))


# write.table( summary(x), file = "ReactomePA_result_diff_Treated_NonTreatedP.Value_0.05.txt", sep = "\t", row.names = FALSE)

pdf("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffSPIApval005/P.Value_0.05_diff.pdf")
plotP(spia_result_diff_Treated_NonTreated, threshold=0.05 )
dev.off()


# head(spia_result_diff_Treated_NonTreated)
# jpeg("_diff.jpeg")
# plotP(spia_result_diff_Treated_NonTreated, threshold=0.05) 
# dev.off()

# setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffSPIApval005")
# pdf("_diff.pdf")
# plotP(spia_result_diff_Treated_NonTreated, threshold=1) 
# dev.off()
 
# vec_hsa_diff_Treated_NonTreated =paste("hsa",spia_result_diff_Treated_NonTreated$ID,sep="")

#  require(pathview)
#  out.suffix="spiadiff_Treated_NonTreated"
#   setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffSPIApval005")
# pv.out.list_diff <- sapply(vec_hsa_diff_Treated_NonTreated, 
#   function(pid) pathview(gene.data = sig_genes_diff_Treated_NonTreated, 
#     pathway.id = pid,species = "hsa",
#     low = list(gene = "green", cpd = "blue"), 
#     mid = list(gene = "gray", cpd = "gray"), 
#     high = list(gene = "red", cpd = "yellow"), 
#     out.suffix=out.suffix,gene.idtype="entrez"))
#  setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffSPIApval005")

# write.table( spia_result_diff_Treated_NonTreated, file = "spia_result_diff_Treated_NonTreateddiff.txt", sep = "\t", row.names = FALSE)
# write.table( vec_hsa_diff_Treated_NonTreated, file = "vec_hsa_diff_Treated_NonTreateddiff.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


# # http://uhts-lgtf.vital-it.ch/account/logout 
# ibenmessaoud
# dU8ge64R


