#########################2019-04-19##########################################
##find first junction position (exon1end-exonstart) in hg38 annotation file##
#############################################################################
library(dplyr)
library(tidyr)

##read data##
hg38<-read.delim("hg38_genelnRNA_sorted.gtf",header=FALSE,stringsAsFactors = FALSE)
hg38exon<-hg38[hg38$V3=="exon",]
hg38exon<-hg38exon%>%separate(V9,c(paste0("V",9:26)),sep=";")
write.table(hg38exon,"hg38exon.txt",col.names = FALSE,row.names=FALSE,quote=FALSE,sep="\t")

##different format##
unique(hg38exon$V2) #"HAVANA"-manual annotate    "ENSEMBL"--pipeline annotate   "StringTie"--self-annotate

##subset by format##
##StringTie##UNKNOWN lncRNA##
hg38exonS<-hg38exon[hg38exon$V2=="StringTie",]
for (i in 1:nrow(hg38exonS)){
  hg38exonS$gene_name[i]<-strsplit(hg38exonS$V9[i]," ")[[1]][2]
  
}

hg38exonS<-hg38exonS[,-c(12:26)]
hg38exonS<-as_tibble(hg38exonS)

hg38exonS1<-hg38exonS[hg38exonS$V11==" exon_number 1",]
hg38exonS2<-hg38exonS[hg38exonS$V11==" exon_number 2",]

hg38exonS1.2<-full_join(hg38exonS1,hg38exonS2,by="V10")
hg38exonS1.2<-na.omit(hg38exonS1.2)
hg38exonS1.2$V1.y<-hg38exonS1.2$V2.y<-hg38exonS1.2$V3.y<-hg38exonS1.2$V6.y<-hg38exonS1.2$V8.y<-hg38exonS1.2$V9.y<-NULL

hg38juncS<-hg38exonS1.2[,c(1:3,5,7,9,10,13)]
hg38juncS<-hg38juncS%>%rename(chrom=V1.x,start=V5.x,end=V4.y,junc=V3.x,method=V2.x,strand=V7.x,gene_id=V9.x,transcript_id=V10)
hg38juncS$junc<-"junction"
hg38juncS<-hg38juncS[,c(1:4,8,5:7)]
write.table(hg38juncS,"hg38_junction_stringtie.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

##HAVANA or ENSEMBLE##KNOWN##
hg38exonH<-hg38exon[hg38exon$V2=="HAVANA"|hg38exon$V2=="ENSEMBL",]
hg38exonH1<-hg38exonH[hg38exonH$V12==" exon_number 1"|hg38exonH$V13==" exon_number 1",]
hg38exonH2<-hg38exonH[hg38exonH$V12==" exon_number 2"|hg38exonH$V13==" exon_number 2",]

hg38exonH1.2<-full_join(hg38exonH1,hg38exonH2[,c("V4","V10")],by="V10")
hg38juncH<-hg38exonH1.2[,c(1:3,5,27,7,13:16)]
hg38juncH<-na.omit(hg38juncH)
hg38juncH<-hg38juncH%>%rename(chrom=V1,method=V2,junc=V3,start=V5,end=V4.y,gene_id=V13,gene_name=V14,gene_status=V15,gene_type=V16,strand=V7)
hg38juncH$junc<-"junction"
write.table(hg38juncH,"hg38_junction_HAVANA.ENSEMBL.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")










