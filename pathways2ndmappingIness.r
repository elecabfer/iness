library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
 library(limma)
 library("RDAVIDWebService")

 david = DAVIDWebService$new(email="sara.benmohammed@epfl.ch")

setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/")
data1 = read.delim("genes_expression_GEO_HEC.txt" )
data2 = read.delim("genes_expression_LL_TT_12.txt" )
data.rpkm = data.frame(merge(data1[,c(1,grep("rpkm",names(data1)))],
  data2[,c(1,grep("rpkm",names(data2)))]),row.names=1)
I = which(apply(data.rpkm,1,min)>1e-2)
# comment est fait le choix 1e-2 ???
data.rpkm = data.rpkm[I,grep ("HEC|LL|TT",colnames( data.rpkm))]
 
# attention ça change le résultat de LIMMA
names( data.rpkm) =
  c( "HEC_B1.1", "HEC_B1.2", "HEC_B1.3",   "HEC_B1.4", "HEC_B1.5",  
 "LL_B2.1",  "LL_B2.2",  "LL_B2.3",  "LL_B2.4",  "LL_B2.5",  "LL_B2.6",   "LL_B2.7",  "LL_B2.8",  
 "TT_B2.1",  "TT_B2.2",  "TT_B2.3",  "TT_B2.4",   "TT_B2.5",  "TT_B2.6",  "TT_B2.7", 
 "TT_B1.1",  "TT_B1.2",  "TT_B1.3",   "TT_B1.4",  "TT_B1.5",  
 "LL_B1.1",  "LL_B1.2",  "LL_B1.3",  "LL_B1.4",   "LL_B1.5",  "LL_B1.6",  "LL_B1.7",   "LL_B1.8")

#### Median  normalization
meds = apply(data.rpkm,2,median,na.rm=T)
data.norm = sweep(data.rpkm,2,meds,"/")
data.norm = log2(data.norm)



# analyse limma

condition.factor <- factor(c("HEC","HEC","HEC","HEC","HEC",
 "LL","LL","LL","LL","LL","LL" ,
"LL","LL","TT","TT","TT","TT"   ,
"TT","TT","TT"  ,
"TT","TT","TT" ,
 "TT","TT","LL","LL","LL","LL" ,
"LL","LL","LL","LL"), levels = c("HEC","LL" ,"TT" ))
 
batch.factor <- factor(c("B1","B1","B1","B1","B1",
 "B2","B2","B2","B2","B2","B2" ,
"B2","B2","B2","B2","B2","B2"   ,
"B2","B2","B2"  ,
"B1","B1","B1" ,
 "B1","B1","B1","B1","B1","B1" ,
"B1","B1","B1","B1"), levels = c("B1","B2"  ))
 
exp.eset.rm.batch<-removeBatchEffect(data.norm,batch.factor)
 # only consider normal and disease conditions in the model matrix 
 design<-model.matrix(~0+condition.factor) # fit the linear model
colnames(design) <- levels(condition.factor)
 fit<-lmFit(exp.eset.rm.batch,design) 


contrast.matrix <- makeContrasts( 
contrasts=c("LL-TT",
"(LL+TT)/2-(HEC)"   ), 
     levels=design)
contrast.matrix

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
 fit2$Ameans=rowMeans(exp.eset.rm.batch)

# # # # # # # # # # # # # #pathways (LL+TT)/2-(HEC) # # # #    

top_LLTT_HEC <- topTable(fit2, coef="(LL+TT)/2-(HEC)", number=nrow(fit2) )

 
top_LLTT_HECup= rownames(   top_LLTT_HEC [(( top_LLTT_HEC$logFC >  0  ) & (top_LLTT_HEC$P.Value <0.05)), ] )
top_LLTT_HECdown= rownames(  top_LLTT_HEC [(( top_LLTT_HEC$logFC <  0   ) & (top_LLTT_HEC$P.Value <0.05)) ,] )
  
 
# # # # # # # # # # # # # #pathways LL/TT # # # #    

top_LLTT <- topTable(fit2, coef="LL-TT", number=nrow(fit2) )
  
top_LLTTup= rownames(  top_LLTT[(( top_LLTT$logFC >  0  ) & (top_LLTT$P.Value <0.05)) ,] )
top_LLTTdown= rownames(  top_LLTT[(( top_LLTT$logFC <  0   ) & (top_LLTT$P.Value <0.05)) ,] )
 

# # # # # # # # # # # # # #pathways   # # # #    



# DAVID pathway research


getAllAnnotationCategoryNames (david)
#  [1] "BBID"                                
#  [2] "BIND"                                
#  [3] "BIOCARTA"                            
#  [4] "BLOCKS"                              
#  [5] "CGAP_EST_QUARTILE"                   
#  [6] "CGAP_SAGE_QUARTILE"                  
#  [7] "CHROMOSOME"                          
#  [8] "COG_NAME"                            
#  [9] "COG_ONTOLOGY"                        
# [10] "CYTOBAND"                            
# [11] "DIP"                                 
# [12] "EC_NUMBER"                           
# [13] "ENSEMBL_GENE_ID"                     
# [14] "ENTREZ_GENE_ID"                      
# [15] "ENTREZ_GENE_SUMMARY"                 
# [16] "GENETIC_ASSOCIATION_DB_DISEASE"      
# [17] "GENERIF_SUMMARY"                     
# [18] "GNF_U133A_QUARTILE"                  
# [19] "GENETIC_ASSOCIATION_DB_DISEASE_CLASS"
# [20] "GOTERM_BP_2"                         
# [21] "GOTERM_BP_1"                         
# [22] "GOTERM_BP_4"                         
# [23] "GOTERM_BP_3"                         
# [24] "GOTERM_BP_FAT"                       
# [25] "GOTERM_BP_5"                         
# [26] "GOTERM_CC_1"                         
# [27] "GOTERM_BP_ALL"                       
# [28] "GOTERM_CC_3"                         
# [29] "GOTERM_CC_2"                         
# [30] "GOTERM_CC_5"                         
# [31] "GOTERM_CC_4"                         
# [32] "GOTERM_MF_1"                         
# [33] "GOTERM_MF_2"                         
# [34] "GOTERM_CC_FAT"                       
# [35] "GOTERM_CC_ALL"                       
# [36] "GOTERM_MF_5"                         
# [37] "GOTERM_MF_FAT"                       
# [38] "GOTERM_MF_3"                         
# [39] "GOTERM_MF_4"                         
# [40] "HIV_INTERACTION_CATEGORY"            
# [41] "HIV_INTERACTION_PUBMED_ID"           
# [42] "GOTERM_MF_ALL"                       
# [43] "HIV_INTERACTION"                     
# [44] "KEGG_PATHWAY"                        
# [45] "HOMOLOGOUS_GENE"                     
# [46] "INTERPRO"                            
# [47] "OFFICIAL_GENE_SYMBOL"                
# [48] "NCICB_CAPATHWAY_INTERACTION"         
# [49] "MINT"                                
# [50] "PANTHER_MF_ALL"                      
# [51] "PANTHER_FAMILY"                      
# [52] "PANTHER_BP_ALL"                      
# [53] "OMIM_DISEASE"                        
# [54] "PFAM"                                
# [55] "PANTHER_SUBFAMILY"                   
# [56] "PANTHER_PATHWAY"                     
# [57] "PIR_SUPERFAMILY"                     
# [58] "PIR_SUMMARY"                         
# [59] "PIR_SEQ_FEATURE"                     
# [60] "PROSITE"                             
# [61] "PUBMED_ID"                           
# [62] "REACTOME_INTERACTION"                
# [63] "REACTOME_PATHWAY"                    
# [64] "PIR_TISSUE_SPECIFICITY"              
# [65] "PRINTS"                              
# [66] "PRODOM"                              
# [67] "PROFILE"                             
# [68] "SMART"                               
# [69] "SP_COMMENT"                          
# [70] "SP_COMMENT_TYPE"                     
# [71] "SP_PIR_KEYWORDS"                     
# [72] "SCOP_CLASS"                          
# [73] "SCOP_FAMILY"                         
# [74] "SCOP_FOLD"                           
# [75] "SCOP_SUPERFAMILY"                    
# [76] "UP_SEQ_FEATURE"                      
# [77] "UNIGENE_EST_QUARTILE"                
# [78] "ZFIN_ANATOMY"                        
# [79] "UP_TISSUE"                           
# [80] "TIGRFAMS"                            
# [81] "SSF"                                 
# [82] "UCSC_TFBS"                           


setAnnotationCategories(david, c("KEGG_PATHWAY"))

 
result1 = addList(david, top_LLTT_HECup,idType="ENSEMBL_GENE_ID", listName="LTT_HECup", listType="Gene")
enriched_pathways1=getFunctionalAnnotationChart(david)
#  enriched_pathways1
# DAVID Result object
# Result type:  FunctionalAnnotationChart 
#        Category                                               Term Count
# 1  KEGG_PATHWAY hsa04650:Natural killer cell mediated cytotoxicity    10
# 2  KEGG_PATHWAY                   hsa04330:Notch signaling pathway     6
# 3  KEGG_PATHWAY         hsa04660:T cell receptor signaling pathway     8
# 4  KEGG_PATHWAY                    hsa04370:VEGF signaling pathway     6
# 5  KEGG_PATHWAY                        hsa05213:Endometrial cancer     5
# 6  KEGG_PATHWAY                     hsa04310:Wnt signaling pathway     8
# 7  KEGG_PATHWAY                           hsa05215:Prostate cancer     6
# 8  KEGG_PATHWAY hsa04130:SNARE interactions in vesicular transport     4
# 9  KEGG_PATHWAY                    hsa04720:Long-term potentiation     5
# 10 KEGG_PATHWAY       hsa04610:Complement and coagulation cascades     5
# 11 KEGG_PATHWAY         hsa04662:B cell receptor signaling pathway     5
# 12 KEGG_PATHWAY                               hsa04144:Endocytosis     8
# 13 KEGG_PATHWAY                         hsa05210:Colorectal cancer     5

# patient_healthy up
 setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/up/figures")
 pdf("LLTT_HEC_up_KEGG_DAVID.pdf")
par(mar=c(5, 20, 4, 2) + 0.1)
barplot( enriched_pathways1$PValue , names.arg= enriched_pathways1$Term ,col="red",horiz=T,las=2,cex.lab=1,cex.names=0.8,xlab="P value",xlim=c(0,(0.02+max(enriched_pathways1$PValue))))
title("KEGG pathways enriched\nLLTT/HEC up-regulated genes")
dev.off()
 setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/up")
write.table(enriched_pathways1, file = "enriched_pathways1patient_healthy.txt",  sep = "\t" )

top_LLTT_HECup_annotation=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description","chromosome_name" ,"start_position","end_position"),
  filters="ensembl_gene_id",values=top_LLTT_HECup, mart=ensembl)
write.table(top_LLTT_HECup_annotation , file = "up_Genes_patient_healthy.txt",  sep = "\t" )
 

# patient_healthy down
result2 = addList(david, top_LLTT_HECdown,idType="ENSEMBL_GENE_ID", listName="LTT_HECdown", listType="Gene")
enriched_pathways2=getFunctionalAnnotationChart(david)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/down/figures")
 pdf("LLTT_HEC_down_KEGG_DAVID.pdf")
 par(mar=c(5, 20, 4, 2) + 0.1)
barplot( enriched_pathways2$PValue , names.arg= enriched_pathways2$Term ,col="blue",horiz=T,las=2,cex.lab=1,cex.names=0.8,xlab="P value",xlim=c(0,(0.02+max(enriched_pathways2$PValue))))
title("KEGG pathways enriched\nLLTT/HEC down-regulated genes")
dev.off()
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/down")
write.table(enriched_pathways2, file = "enriched_pathways2patient_healthy.txt",  sep = "\t" )
top_LLTT_HECdown_annotation=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description","chromosome_name" ,"start_position","end_position"),
  filters="ensembl_gene_id",values=top_LLTT_HECdown, mart=ensembl)
write.table(top_LLTT_HECdown_annotation , file = "down_Genes_patient_healthy.txt",  sep = "\t" )

# LL12_TT12 up
 result3 = addList(david, top_LLTTup,idType="ENSEMBL_GENE_ID", listName="LLTT_up", listType="Gene")
enriched_pathways3=getFunctionalAnnotationChart(david)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/up/figures")
pdf("LL_TT_up_KEGG_DAVID.pdf")
par(mar=c(5, 20, 4, 2) + 0.1)
barplot( enriched_pathways3$PValue , names.arg= enriched_pathways3$Term ,col="red",horiz=T,las=2,cex.lab=1,cex.names=0.8,xlab="P value",xlim=c(0,(0.02+max(enriched_pathways3$PValue))))
title("KEGG pathways enriched\nLL/TT up-regulated genes")
dev.off()
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/up")
write.table(enriched_pathways3, file = "enriched_pathways3.txt",  sep = "\t" )
top_LLTTup_annotation=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description","chromosome_name" ,"start_position","end_position"),
  filters="ensembl_gene_id",values=top_LLTTup, mart=ensembl)
write.table(top_LLTTup_annotation , file = "up_Genes_LL_TT.txt",  sep = "\t" )

# LL12_TT12 down
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/down/figures")
result4 = addList(david, top_LLTTdown,idType="ENSEMBL_GENE_ID", listName="LLTT_down", listType="Gene")
enriched_pathways4=getFunctionalAnnotationChart(david)
pdf("LL_TT_down_KEGG_DAVID.pdf")
par(mar=c(5, 20, 4, 2) + 0.1)
barplot( enriched_pathways4$PValue , names.arg= enriched_pathways4$Term ,col="blue",horiz=T,las=2,cex.lab=1,cex.names=0.8,xlab="P value",xlim=c(0,(0.02+max(enriched_pathways4$PValue))))
title("KEGG pathways enriched\nLL/TT down-regulated genes")
dev.off()
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/down")
write.table(enriched_pathways4, file = "enriched_pathways4.txt",  sep = "\t" )
top_LLTTdown_annotation=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description","chromosome_name" ,"start_position","end_position"),
  filters="ensembl_gene_id",values=top_LLTTdown, mart=ensembl)
write.table(top_LLTTdown_annotation , file = "down_Genes_LL_TT.txt",  sep = "\t" )


 david2 = DAVIDWebService$new(email="sara.benmohammed@epfl.ch")
setAnnotationCategories(david2, c("GENETIC_ASSOCIATION_DB_DISEASE_CLASS"))

 resultDISEASE1 = addList(david2, top_LLTT_HECup,idType="ENSEMBL_GENE_ID", listName="LTT_HECup", listType="Gene")
enriched_DISEASE1=getFunctionalAnnotationChart(david2)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/up")
write.table(enriched_DISEASE1, file = "GENETIC_ASSOCIATION_DB_DISEASE_CLASS1patient_healthy.txt",  sep = "\t" )

resultDISEASE2 = addList(david2, top_LLTT_HECdown,idType="ENSEMBL_GENE_ID", listName="LTT_HECdown", listType="Gene")
enriched_DISEASE2=getFunctionalAnnotationChart(david2)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/patient_healthy/down")
write.table(enriched_DISEASE2, file = "GENETIC_ASSOCIATION_DB_DISEASE_CLASS2patient_healthy.txt",  sep = "\t" )

 resultDISEASE3 = addList(david2, top_LLTTup,idType="ENSEMBL_GENE_ID", listName="LLTT_up", listType="Gene")
enriched_DISEASE3=getFunctionalAnnotationChart(david2)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/up")
write.table(enriched_DISEASE3, file = "GENETIC_ASSOCIATION_DB_DISEASE_CLASS3LL12_TT12.txt",  sep = "\t" )

resultDISEASE4 = addList(david2, top_LLTTdown,idType="ENSEMBL_GENE_ID", listName="LLTT_down", listType="Gene")
enriched_DISEASE4=getFunctionalAnnotationChart(david2)
setwd("/Users/benmoham/Desktop/annotation_Charlotte/NewDataDec2014/2ndmapping/Annotations/pathwaysSansGEO/DAVID_pvalue0.05/LL12_TT12/down")
write.table(enriched_DISEASE4, file = "GENETIC_ASSOCIATION_DB_DISEASE_CLASS4LL12_TT12.txt",  sep = "\t" )



#           1. *Category*: factor with the main categories under used in
#               the present analysis.

#            2. *Term*: character with the name of the term in format
#               id~name (if available).
# :

#            3.  *Count*: integer with the number of ids of the gene list
#               that belong to this term.

#            4. *X.*: after converting user input gene IDs to
#               corresponding DAVID gene ID, it refers to the percentage
#               of DAVID genes in the list associated with a particular
#               annotation term.  Since DAVID gene ID is unique per gene,
#               it is more accurate to use DAVID ID percentage to present
#               the gene-annotation association by removing any
#               redundancy in user gene list, i.e. two user IDs represent
#               same gene.

#            5. *PValue*: numeric with the EASE Score of the term (see
#               DAVID Help page).

#            6. *Genes*: character in comma separated style with the
#               genes present in the term.

#            7. *List.Total, Pop.Hits, Pop.Total*: integers (in addition
#               to Count) to build the 2x2 contingency table in order to
#               compute the EASE Score (see DAVID Help page).

#            8. *Fold.Enrichment*: numeric with the ratio of the two
#               proportions. For example, if 40/400 (i.e. 10%) of your
#               input genes involved in "kinase activity" and the
#               background information is 300/30000 genes (i.e. 1%)
#               associating with "kinase activity", roughly 10% / 1% = 10
#               fold enrichment.

#            9.  *Bonferroni, Benjamini, FDR*: numerics with p-value
#               adjust different criteria (see p.adjust).






setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data")
data1 = read.delim("genes.htseq.normalized.txt" )

par(mfrow=c(2,1), mar=c(6,2,1,1))
 pdf("boxplots.pdf")
boxplot( log2(data1[,2:13]),pch=20,outline=F,lty=1,las=2,  
col=c( "blue",   "blue",
	"red", "red",
	"grey40",   "grey40",
  "violet" , "violet", 
  "orange", "orange" ,
  "grey50" ,"grey50" )   ,
 cex.axis=0.7, main="log2(normalized)")
boxplot(data1[,2:13],pch=20,outline=F,lty=1,las=2,  
col=c( "blue",   "blue",
	"red", "red",
	"grey40",   "grey40",
  "violet" , "violet", 
  "orange", "orange" ,
  "grey50" ,"grey50" )    ,
 cex.axis=0.7, main=" normalized ")
dev.off()

 "X12NT"   "X9NT"       # negative control 3 days
 "X11DXR"  "X8DXR"      # positive control 3 days  
    "X10L100" "X7L100"  # treated 3 days  
    "X3NT"   "X6NT"     # negative control 6h 
    "X2DXR"   "X5DXR"   # positive control 6h 
     "X1L100"  "X4L100" # Treated 6h
# 
12NT	negative control 3 days b	Ctrlneg_3d
9NT	negative control 3 days a	Ctrlneg_3d
11DXR	positive control 3 days b	Ctrlpos_3d
8DXR	positive control 3 days a	Ctrlpos_3d
10L100	treated 3 days b	Treated_3d
7L100	treated 3 days a	Treated_3d
3NT	negative control 6h a	Ctrlneg_6h
6NT	negative control 6h b	Ctrlneg_6h
2DXR	positive control 6h a	Ctrlpos_6h
5DXR	positive control 6h b	Ctrlpos_6h
1L100	Treated 6h a	Treated_6h
4L100	Treated 6h b	Treated_6h


setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data")
data1 = read.delim("genes.htseq.normalized.txt" )


 pdf("boxplots.pdf")
boxplot(data1[,2:13],pch=20,outline=F,lty=1,las=2,  
col=c( "blue",   "blue",
	"red", "red",
	"grey40",   "grey40",
  "violet" , "violet", 
  "orange", "orange" ,
  "grey50" ,"grey50" )    ,
 cex.axis=0.7, main="TMM Normalized counts ")
dev.off()

pdf("pairsplot.pdf")
pairs(data1[,c("X12NT", "X9NT","X10L100","X7L100")],pch=20,outline=F,lty=1,las=2   ,
 cex.axis=0.7, main="TMM Normalized counts ")
dev.off()


pca=prcomp(t (data1[,2:13]) , cor=TRUE, scores=TRUE)
pdf("pca_norm.pdf")
plot (pca$x[,1],pca$x[,2], pch=20 , col=c( "blue",   "blue",
	"red", "red",
	"grey40",   "grey40",
  "violet" , "violet", 
  "orange", "orange" ,
  "grey50" ,"grey50" )     )
   
text (pca$x[,1],pca$x[,2],labels=substr(names(pca$x[,1]),2,18), col=c( "blue",   "blue",
	"red", "red",
	"grey40",   "grey40",
  "violet" , "violet", 
  "orange", "orange" ,
  "grey50" ,"grey50" )     )
dev.off()


pca=prcomp(t (data1[,c("X12NT", "X9NT","X10L100","X7L100")]) , cor=TRUE, scores=TRUE)
pdf("pca2_norm.pdf")
plot (pca$x[,1],pca$x[,2], pch=20 , 
	col=c( "blue",   "blue",
	"red", "red")     )
   
text (pca$x[,1],pca$x[,2],labels=substr(names(pca$x[,1]),2,18), 
	col=c( "blue",   "blue",
	"red", "red")     )
dev.off()

data_expr = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genelist_Treated_3d-Ctrlneg_3d.txt")

names(data_expr)
# [1] "ID"                   "logFC"                "AveExpr"             
#  [4] "AveExp_Treated_3d"    "AveExp_Ctrlneg_3d"    "t"                   
#  [7] "P.Value"              "adj.P.Val"            "Fold_Change"         
# [10] "Associated_Gene_Name" "Gene_Biotype"         "Description"         
# [13] "Chromosome_Name"      "Gene_Start"           "Gene_End"            
# [16] "Strand"               "Transcript_count"     "Gene_Status"         
# [19] "Gene_Source"         

# 12NT	negative control 3 days b	Ctrlneg_3d
# 9NT	negative control 3 days a	Ctrlneg_3d

# 10L100	treated 3 days b	Treated_3d
# 7L100	treated 3 days a	Treated_3d

A_Treated_Ctrlneg= (data_expr$AveExpr)
M_Treated_Ctrlneg=data_expr$logFC

A_Treated_Ctrlneg_diff= (data_expr$AveExpr[data_expr$adj.P.Val<0.05])
M_Treated_Ctrlneg_diff=data_expr$logFC[data_expr$adj.P.Val<0.05]

# ENSG00000108556
# "ID"	"AveExp_Treated_3d"	"AveExp_Ctrlneg_3d"	
# "ENSG00000108556"	5.25895257726863	0.976751190196135		

pdf("MAplot.pdf")
plot(A_Treated_Ctrlneg,M_Treated_Ctrlneg,xlab="AveExpr", ylab="logFC" ,  
	col = gray.colors(1), pch=20,cex =0.5)
points(A_Treated_Ctrlneg_diff,M_Treated_Ctrlneg_diff,   col = "red",pch=20,cex =0.5)
abline(h =1, bg = "black")
abline(h =-1, bg = "black") 
abline(h =0, bg = "black") 
dev.off()

library(gplots)

input.i  <-list(OverExpressedGenes=data_expr$ID[(data_expr$adj.P.Val<0.05)& (data_expr$logFC>0 )],
	UnderExpressedGenes=data_expr$ID[(data_expr$adj.P.Val<0.05)& (data_expr$logFC<0 )])
pdf("VennNumber_of_genes_analysed.pdf")
venn(input.i )
dev.off()



data1spia = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffspia/vec_hsa_diff_Treated_NonTreateddiff.txt",header=FALSE )
data1gage_greater = substr(read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffgage/fc.kegg.pgreater_Treated_NonTreated.txt" )$pathway,1,8)
data1gage_less = substr(read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffgage/fc.kegg.pless_Treated_NonTreated.txt" )$pathway,1,8)
data1spiapval005 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffSPIApval005/vec_hsa_diff_Treated_NonTreateddiff.txt",header=FALSE )
 
input.i  <-list("spia\npval<0.05"=data1spia,"spia\nadjpval<005"=data1spiapval005,
	gage_greater=data1gage_greater,
	gage_less=data1gage_less)

pdf("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffSPIApval005/VennNumber_of_pathways.pdf")
venn(input.i , small=1 )
dev.off()








