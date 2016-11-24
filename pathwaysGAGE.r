 library(sqldf)
library(biomaRt)
 require(gage)
 data(kegg.gs)
 
# # # # # # # # # # # # # # #pathways (LL+TT)/2-(HEC) # # # #    

# top_Treated_NonTreated <- data1
# mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 


# eattr<-c("ensembl_gene_id","entrezgene" )
# bpres<-getBM(attributes=eattr,filters="ensembl_gene_id",values= top_Treated_NonTreated$ID,mart= mart)
# results1_Treated_NonTreated <-bpres[bpres$entrezgene!="",]
    
# res1_top_pLNcells_Hepatocytes <- sqldf("select *
#               from top_Treated_NonTreated 
#               LEFT OUTER JOIN   results1_Treated_NonTreated
#               ON top_Treated_NonTreated.ID=results1_Treated_NonTreated.ensembl_gene_id;")
# res1_top_pLNcells_Hepatocytes <-res1_top_pLNcells_Hepatocytes [grep("TRUE",res1_top_pLNcells_Hepatocytes$entrezgene!="NA"),]
# res1_top_pLNcells_Hepatocytes = res1_top_pLNcells_Hepatocytes[!duplicated(res1_top_pLNcells_Hepatocytes$entrezgene),]
# #remove duplicates
# # dim(res1_top_pLNcells_Hepatocytes)
# # [1] 219   7
#   setwd("/Users/benmoham/Desktop/InessBenmessaoud/diffgage")
# # 
# sig_genes_diff_Treated_NonTreated <-  res1_top_pLNcells_Hepatocytes$logFC 
# names(sig_genes_diff_Treated_NonTreated) <-  res1_top_pLNcells_Hepatocytes$entrezgene 
# # length(sig_genes_diff_Treated_NonTreated)
# # [1] 177
# length(names(sig_genes_diff_Treated_NonTreated))
# # 219
# length(unique(names(sig_genes_diff_Treated_NonTreated)))
# # [1] 218
 setwd("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/")

data1 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genelist_Treated_3d-Ctrlneg_3d.txt" )
ensembl <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 
table1 = getBM(
  attr=c("ensembl_gene_id","entrezgene","description","hgnc_symbol"),
  filter=c("ensembl_gene_id"),
  values=data1$ID ,
  mart=ensembl)

# deseq.fc=data1$logFC[grep("TRUE",data1$logFC<0)]
# gnames=data1$ID[grep("TRUE",data1$logFC<0)]
deseq.fc=data1$logFC 
gnames=data1$ID 

 gnames=as.character(gnames )
     names(deseq.fc)=gnames
 library(pathview)

gnames.eg=pathview::id2eg(gnames, category =gene.idtype.list[3], org = "Hs")

 sel2=gnames.eg[,2]>""
  deseq.fc=deseq.fc[sel2]
 names(deseq.fc)=as.character(gnames.eg[sel2,2])
exp.fc=deseq.fc

gnames.eg[,2]
# exp.fc<-exp.fc [grep(";", names(deseq.fc) ,invert="TRUE") ]
# exp.fc
out.suffix="T_NT"
 require(gage)
 data(kegg.gs)
 fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
 # go.hs=go.gsets(species="human")
 # gage(exp.fc, gsets = go.gs, ref = NULL, samp = NULL,go.subs="BP")
#  data(kegg.gs)
# data(kegg.gs.dise)
# data(go.gs)
# data(carta.gs)
# BioCarta had not been updating its pathways. The information provided might have been outdated. As a result, we have discontinued offering pathway information online. You may view our pathway figures at http://cgap.nci.nih.gov/Pathways/BioCarta_Pathways. If you are interested in using some of its pathway figures, please contact info@biocarta.com for permission.

# Format
   
 sel <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "p.val"])
 path.ids <- rownames(fc.kegg.p$greater)[sel]
 sel.l <- fc.kegg.p$less[, "p.val"] < 0.05 & !is.na(fc.kegg.p$less[, "p.val"])
 path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids2 <- substr(c(path.ids[1:3], path.ids.l[1:3]), 1, 8)
 path.ids2_Treated_NonTreated <-  c(path.ids , path.ids.l )
  path.ids2_Treated_NonTreated
#  
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGELogFC/up")
write.table(fc.kegg.p$greater[sel,3:4], file = "fc.kegg.p$greater_Treated_NonTreated.txt",  sep = "\t")
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGELogFC/down")
write.table( fc.kegg.p$less[sel.l,3:4], file = "fc.kegg.p$less_Treated_NonTreated.txt", sep = "\t")
 

 require(pathview)
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGELogFC/up")
pv.out.listup <- sapply(substr(c(path.ids  ), 1, 8), function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa",
  low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd = "yellow"),  
    out.suffix=out.suffix,gene.idtype="entrez", kegg.native = F  ) )
 
 data1 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genelist_Treated_3d-Ctrlneg_3d.txt" )


  ensembl <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org")

  TABLE_CC= 
getBM( attr=c("ensembl_gene_id", "name_1006","namespace_1003"),
  filter=c("ensembl_gene_id"),values= unique(data1$ID)   ,
  mart=ensembl) 

ensembl_gene_id_CC=sqldf("select ensembl_gene_id,name_1006
  from TABLE_CC 
  WHERE TABLE_CC.namespace_1003==\"cellular_component\"
  group by ensembl_gene_id")
dim(ensembl_gene_id_CC)
# ensembl_gene_id_CC
 colnames(ensembl_gene_id_CC)=c("ensembl_id", "name_1006" )  

 table1 = getBM( attr=c("ensembl_gene_id","entrezgene","description","hgnc_symbol"),
  filter=c("ensembl_gene_id"),values=data1$ID ,mart=ensembl)
 
ensembl_gene_id_annotation=sqldf("SELECT  *
FROM table1
LEFT OUTER JOIN ensembl_gene_id_CC
ON table1.ensembl_gene_id=ensembl_gene_id_CC.ensembl_id;
")

for (i in c(1:(length(do.call("cbind",pv.out.listup)["plot.data.gene", ]) ))) { 
    print( names(do.call("cbind",pv.out.listup)["plot.data.gene", ]) [i])
    if ( (do.call("cbind",pv.out.listup)["plot.data.gene",i ]) [[1]] !=0){ 
    write.table(
        merge( 
    as.data.frame(do.call("cbind",pv.out.listup)["plot.data.gene",i ]    ) ,
     as.data.frame(ensembl_gene_id_annotation) ,
     by.x = "labels"  , by.y = "hgnc_symbol", 
     all.x=TRUE, all.y=FALSE)
      ,sep="\t"  ,append=TRUE,
      file=paste( names(do.call("cbind",pv.out.listup)["plot.data.gene", ]) [i],"_genes.txt",sep=""),
      row.names=FALSE
    )   }
}
  

setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGELogFC/down")
 
 pv.out.listdown <- sapply(substr(c(path.ids.l   ), 1, 8), function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa",
  low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd = "yellow"),  
    out.suffix=out.suffix,gene.idtype="entrez", kegg.native = F  ) )
 

 
for (i in c(1:(length(pv.out.listdown[1,] ) ))) { 
    print( names(pv.out.listdown[1,]) [i])
    if ( length(pv.out.listdown[1,]  [i] )>0){ 
    write.table(
        merge( 
    as.data.frame(pv.out.listdown[1,i ]    ) ,
     as.data.frame(ensembl_gene_id_annotation) ,
     by.x = "labels"  , by.y = "hgnc_symbol", 
     all.x=TRUE, all.y=FALSE)
      ,sep="\t"  ,append=TRUE,
      file=paste( names(pv.out.listdown[1,]) [i],"_genes.txt",sep=""),
      row.names=FALSE
    )   }
}
  
 
 
 
# data1spia = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffspia/vec_hsa_diff_Treated_NonTreateddiff.txt",header=FALSE )
# data1gage_greater = substr(read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffgage/fc.kegg.pgreater_Treated_NonTreated.txt" )$pathway,1,8)
# data1gage_less = substr(read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/diffgage/fc.kegg.pless_Treated_NonTreated.txt" )$pathway,1,8)




# packageVersion("gage")
# [1] ‘2.16.0’
# Note that kegg.gs has been updated since gage version
#      2.9.1. From then, kegg.gs only include the subset of canonical
#      signaling and metabolic pathways from KEGG pathway database, and
#      kegg.gs.dise is the subset of disease pathways.

 

###################################################
 
#####MASS PROTEIN##################################
############
Treated_NT <- read.table("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genelist_Treated_3d-Ctrlneg_3d.txt", header=TRUE,sep="\t")
Treated_NT_genes <- Treated_NT$ID
head(Treated_NT_genes ); 
length(Treated_NT_genes )
# "ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org"
ensembl <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 
 
sequence = getBM(c(seqType=c("peptide"), type="ensembl_gene_id"), 
  filters = c("ensembl_gene_id"), values = Treated_NT_genes, mart = ensembl, 
  checkFilters = FALSE )

# what have we got here
head(sequence)
head(sequence[1])
head(sequence[[2]])

# reformat result and calculate protein lengths
result <- as.data.frame(cbind(sequence[1], sequence[2]))
length <- nchar(result$peptide) - 1
masseDA=110 *(nchar(result$peptide) - 1)
result <- cbind(result,length,masseDA )
I = which( (result$masseDA >=30000)& (result$masseDA <=34000))
data.result = result[I,]

PI=sapply(strsplit(data.result$peptide, "") , computePI)

result <- cbind(data.result,PI)
head(result)
J = which( (result$PI >=5)& (result$PI <=7))
 resultat = result[J ,]
 dim(resultat )
# [1] 995   5
# write out
library(sqldf)

resultatT_NT<- sqldf("select *
                                      from resultat
                                      LEFT OUTER JOIN  Treated_NT 
                                      ON resultat.ensembl_gene_id=Treated_NT.ID;")

write.table(x=resultatT_NT, sep='\t', file="/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/EMBLgene_Prot_ID_with_length_all.txt", row.names=FALSE)
 


# grep '[^,].*,[^,]' /Users/benmoham/Desktop/InessBenmessaoud/FirstFilter281Genes/mart_exportUniProtKB.txt > /Users/benmoham/Desktop/InessBenmessaoud/FirstFilter281Genes/mart_exportEMBLUniProtKB.txt
# grep '[^,].*,[^,]' /Users/benmoham/Desktop/InessBenmessaoud/FirstFilter281Genes/EMBLgene_Prot.txt > /Users/benmoham/Desktop/InessBenmessaoud/FirstFilter281Genes/EMBLgene_Prot_ID.txt

library(biomaRt)
library(seqinr)
# load ensp table
ensp <- read.table("/Users/benmoham/Desktop/InessBenmessaoud/FirstFilter281Genes/EMBLgene_Prot_ID.txt", header=TRUE,sep=",")
ensp <- ensp$Ensembl.Protein.ID
head(ensp); 
length(ensp)
# "ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org"
ensembl <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org") 
ensembl=useMart("ensembl_mart_49",dataset="hsapiens_gene_ensembl",archive=FALSE)
# see how long it takes to fetch sequences of ensps in list
start_time <- proc.time()

sequence = getBM(c(seqType=c("peptide"), type="ensembl_peptide_id"), 
  filters = c("ensembl_peptide_id"), values =  ensp, mart = ensembl, 
  checkFilters = FALSE )

proc.time() - start_time

# what have we got here
head(sequence)
head(sequence[1])
sequence[[2]]

# reformat result and calculate protein lengths
result <- as.data.frame(cbind(sequence[1], sequence[2]))
length <- nchar(result$peptide) - 1
masseDA=110 *(nchar(result$peptide) - 1)
PI=sapply(strsplit(result$peptide, "") , computePI)

result <- cbind(result,length,masseDA,PI)
head(result)


# write out
write.table(x=result, sep='\t', file="/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/EMBLgene_Prot_ID_with_length.txt", row.names=FALSE)
# R CMD install /Users/benmoham/Downloads/seqinr_3.1-3.tar.gz

# write out
write.table(x=result, sep='\t', file="/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/EMBLgene_Prot_ID_with_length_all.txt", row.names=FALSE)

# gage tstudent
/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/EMBLgene_Prot_ID_with_length_all.txt





 
data1 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genes.htseq.normalized.txt" )
   
top=data1 
 
 mart =   useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org")

eattr<-c("ensembl_gene_id","entrezgene" )
bpres<-getBM(attributes=eattr,filters="ensembl_gene_id",values=top$Gene_ID,mart= mart)
results1  <-bpres[bpres$entrezgene!="",]
  

ensembl_gene_id_annotation=sqldf("SELECT  *
FROM top
INNER JOIN results1
ON top.Gene_ID=results1.ensembl_gene_id;
")
# INNER interesction
 out.suffix="Treated_NonTreated"
 require(gage)
 data(kegg.gs)
  
 ensembl_gene_id_annotation_sel<-as.matrix(ensembl_gene_id_annotation[,grep ("X10L100|X7L100|X12NT|X9NT",colnames(ensembl_gene_id_annotation ))])
 row.names( ensembl_gene_id_annotation_sel)<- ensembl_gene_id_annotation$entrezgene

 gage_Treated_NonTreated <- gage( as.matrix( ensembl_gene_id_annotation_sel),
   gsets = kegg.gs, 
      ref = grep ("X10L100|X7L100",colnames( ensembl_gene_id_annotation_sel)), 
     samp = grep ("X12NT|X9NT",colnames( ensembl_gene_id_annotation_sel)),
    compare='paired')  

 
 sel <- gage_Treated_NonTreated$greater[, "p.val"] < 0.05 & !is.na(gage_Treated_NonTreated$greater[, "p.val"])
 path.ids_Treated_NonTreated <- rownames(gage_Treated_NonTreated$greater)[sel]
 sel.l <- gage_Treated_NonTreated$less[, "p.val"] < 0.05 & !is.na(gage_Treated_NonTreated$less[, "p.val"])
 path.ids_Treated_NonTreated.l <- rownames(gage_Treated_NonTreated$less)[sel.l]
 #  path.ids_Treated_NonTreated2 <- substr(c(path.ids_Treated_NonTreated[1:3], path.ids_Treated_NonTreated.l[1:3]), 1, 8)
 # path.ids_Treated_NonTreated2_Treated_NonTreated <-  c(path.ids_Treated_NonTreated , path.ids_Treated_NonTreated.l )
 #  path.ids_Treated_NonTreated2_Treated_NonTreated
#  
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGEStudent/up")
write.table(gage_Treated_NonTreated$greater[sel,3:4], file = "gage_Treated_NonTreated$greater_Treated_NonTreated.txt",  sep = "\t")
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGEStudent/down")
write.table( gage_Treated_NonTreated$less[sel.l,3:4], file = "gage_Treated_NonTreated$less_Treated_NonTreated.txt", sep = "\t")
 

 require(pathview)
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGEStudent/up")
# pv.out.list_Treated_NonTreatedup <- sapply(substr(c(path.ids_Treated_NonTreated  ), 1, 8), function(pid) 
#   pathview(gene.data = ensembl_gene_id_annotation_sel, pathway.id = pid,species = "hsa",
#   low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd = "yellow"),  
#     out.suffix=out.suffix,gene.idtype="entrez", kegg.native = F  ) )

pv.out.list_Treated_NonTreatedup <- sapply(substr(c( path.ids_Treated_NonTreated ), 1, 8), function(pid) 
  pathview(gene.data = ensembl_gene_id_annotation_sel, pathway.id = pid,
  species = "hsa",
  low = list(gene = "green", cpd = "blue"), 
  mid = list(gene = "gray", cpd = "gray"), 
  high = list(gene = "red", cpd = "yellow"),  
    out.suffix=out.suffix,gene.idtype="entrez", kegg.native = F ))
  
   # annotation des genes avec CC
  data1 = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/Benmessaoud_data/genes.htseq.normalized.txt" )
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="www.ensembl.org")

  TABLE_CC= 
getBM( attr=c("ensembl_gene_id", "name_1006","namespace_1003"),
  filter=c("ensembl_gene_id"),values= unique(data1$Gene_ID)   ,
  mart=ensembl) 

ensembl_gene_id_CC=sqldf("select ensembl_gene_id,name_1006
  from TABLE_CC 
  WHERE TABLE_CC.namespace_1003==\"cellular_component\"
  group by ensembl_gene_id")
dim(ensembl_gene_id_CC)
  colnames(ensembl_gene_id_CC)=c("ensembl_id", "name_1006" )  

 table1 = getBM( attr=c("ensembl_gene_id","entrezgene","description","hgnc_symbol"),
  filter=c("ensembl_gene_id"),values=data1$Gene_ID ,mart=ensembl)
 
ensembl_gene_id_annotation=sqldf("SELECT  *
FROM table1
LEFT OUTER JOIN ensembl_gene_id_CC
ON table1.ensembl_gene_id=ensembl_gene_id_CC.ensembl_id;
")
setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGEStudent/up")


 
for (i in c(1:(length(pv.out.list_Treated_NonTreatedup[1,] ) ))) { 
    print( names(pv.out.list_Treated_NonTreatedup[1,]) [i])
    if ( length(pv.out.list_Treated_NonTreatedup[1,]  [i] )>0){ 
    write.table(
        merge( 
    as.data.frame(pv.out.list_Treated_NonTreatedup[1,i ]    ) ,
     as.data.frame(ensembl_gene_id_annotation) ,
     by.x = "labels"  , by.y = "hgnc_symbol", 
     all.x=TRUE, all.y=FALSE)
      ,sep="\t"  ,append=TRUE,
      file=paste( names(pv.out.list_Treated_NonTreatedup[1,]) [i],"_genes.txt",sep=""),
      row.names=FALSE
    )   }
}
  
      
 
  
   

setwd("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/diffgage2/GAGEStudent/down")
 pv.out.list_Treated_NonTreateddown <- sapply(substr(c( path.ids_Treated_NonTreated.l ), 1, 8), function(pid) 
  pathview(gene.data = ensembl_gene_id_annotation_sel, pathway.id = pid,
  species = "hsa",low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), 
  high = list(gene = "red", cpd = "yellow"),  
    out.suffix=out.suffix,gene.idtype="entrez", kegg.native = F ))
 

 
for (i in c(1:(length(do.call("cbind",pv.out.list_Treated_NonTreateddown)["plot.data.gene", ]) ))) { 
    print( names(do.call("cbind",pv.out.list_Treated_NonTreateddown)["plot.data.gene", ]) [i])
    if ( (do.call("cbind",pv.out.list_Treated_NonTreateddown)["plot.data.gene",i ]) [[1]] !=0){ 
    write.table(
        merge( 
    as.data.frame(do.call("cbind",pv.out.list_Treated_NonTreateddown)["plot.data.gene",i ]    ) ,
     as.data.frame(ensembl_gene_id_annotation) ,
     by.x = "labels"  , by.y = "hgnc_symbol", 
     all.x=TRUE, all.y=FALSE)
      ,sep="\t"  ,append=TRUE,
      file=paste( names(do.call("cbind",pv.out.list_Treated_NonTreateddown)["plot.data.gene", ]) [i],"_genes.txt",sep=""),
      row.names=FALSE
    )   }
}
  
 
 GAGELogFCdown = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/GAGELogFC/down/less_Treated_NonTreatedSIGNIFICATIF.txt",header=FALSE )
GAGELogFCup = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/GAGELogFC/up/greater_Treated_NonTreatedSIGNIFICATIF.txt",header=FALSE )
GAGEStudentdown = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/GAGEStudent/down/less_Treated_NonTreatedSIGNIFICATIF.txt",header=FALSE )
GAGEStudentup = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/GAGEStudent/up/greater_Treated_NonTreatedSIGNIFICATIF.txt",header=FALSE )
SPIA = read.delim("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/SPIA/diffSPIApval005/vec_hsa_diff_Treated_NonTreatedSIGNIFICATIF.txt",header=FALSE )

input.i  <- list("GAGE\nLogFC-"= GAGELogFCdown, "GAGE\nLogFC+"=GAGELogFCup,
    "GAGE\nStudent-"=GAGEStudentdown , "GAGE\nStudent+"= GAGEStudentup,  "SPIA"=SPIA)
library(gplots)

 pdf("/Users/benmoham/Desktop/InessBenmessaoud/pvalue_adjusted0.05_2May/pathways_25May2016/VennNumber_pathways.pdf")
venn(input.i )
dev.off()

intersect ( GAGELogFCup,  GAGEStudentdown)
intersect ( GAGELogFCdown,  GAGEStudentup)
 





