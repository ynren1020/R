##################################2019-04-26#############################################
##top expressed known lncRNA DANCR and MALAT1 in BMI1 data###############################
##summarize their first junction:exon1end-exon2start numbers,normalize by total juncs####
##TCGA(PRAD:normal and tumor),SU2C metastic cancer#######################################
##boxplot,data:/data/ryang/tywang/projects/Exitron/SU2C/jannos###########################
##copy hg38_junction_HAVANA.ENSEMBL.txt from RNAseq produced for nat comm study(scratch)#
#########################################################################################

#library(tidyverse)
#library(ggpubr)
library(dplyr)
library(tidyr)

args <- commandArgs(TRUE)
#args[1]<-"hg38_junction_HAVANA.ENSEMBL.txt"

##junction annotation##
hg38junc<-read.delim(args[1],header = TRUE,stringsAsFactors = FALSE)
##subset by lncRNA##
dancr<-filter(hg38junc,gene_name==" gene_name DANCR")
malat1<-filter(hg38junc,gene_name==" gene_name MALAT1")
##junction file of a patient from SU2C##
#su2c<-read.delim(args[2],header=TRUE,stringsAsFactors = FALSE)
##for prad##
su2c<-read.table(gzfile(args[2]),header = TRUE)
total<-sum(su2c$score)
su2c_dancr<-filter(su2c,genes=="DANCR")
su2c_malat1<-filter(su2c,genes=="MALAT1")

##join dancr##
su2c_dancr_join<-left_join(dancr,su2c_dancr,by=c("chrom","start","end","strand"))
su2c_dancr_join<-na.omit(su2c_dancr_join)
su2c_dancr_join<-unique(su2c_dancr_join)
dancr_score<-sum(su2c_dancr_join$score)
dancr_juncNor<-data.frame(sample=args[2],gene="DANCR",score=dancr_score,total=total,scoreN=round(dancr_score/total*10^6,0))
write.table(dancr_juncNor,paste0(args[2],".dancr.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")

##join malat1##
su2c_malat1_join<-left_join(malat1,su2c_malat1,by=c("chrom","start","end","strand"))
su2c_malat1_join<-na.omit(su2c_malat1_join)
su2c_malat1_join<-unique(su2c_malat1_join)
malat1_score<-sum(su2c_malat1_join$score)
malat1_juncNor<-data.frame(sample=args[2],gene="MALAT1",score=malat1_score,total=total,scoreN=round(malat1_score/total*10^6,0))
write.table(malat1_juncNor,paste0(args[2],".malat1.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")
