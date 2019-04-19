############################2019-04-16###########################################
##RNAseq data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120741##
##regtools exact and annotate results of sort and indexed bam file###############
##normalize # junctions at each position by total junctions######################
#################################################################################

library(dplyr)
library(tidyr)

dat<-read.delim("SRR7949386.junction.txt")
dat$scoreN<-(dat$score/sum(dat$score))*(10^6)
junc.lncRNA<-read.delim("hg38_junction_stringtie.txt",stringsAsFactors = FALSE)
junc.gene<-read.delim("hg38_junction_HAVANA.ENSEMBL.txt",stringsAsFactors = FALSE)

####lncRNA####
##full_join##
dat.lncRNA<-full_join(dat,junc.lncRNA,by=c("chrom","start","end"))
dat.lncRNAsub<-dat.lncRNA[!is.na(dat.lncRNA$method)&!is.na(dat.lncRNA$score),]
##filter by median(dat.lncRNA$score)==5##
dat.lncRNAsubfilter<-dat.lncRNAsub[dat.lncRNAsub$scoreN>=median(dat.lncRNAsub$scoreN),]

####genes##
dat.gene<-full_join(dat,junc.gene,by=c("chrom","start","end"))
dat.genesub<-dat.gene[!is.na(dat.gene$method)&!is.na(dat.gene$score),]
##filter by median(dat.gene$score)##
dat.genesubfilter<-dat.genesub[dat.genesub$scoreN>=median(dat.genesub$scoreN),]

##lncRNA/genes have higher junctions keeped for update annotation file##
write.table(dat.lncRNAsubfilter,"dat.lncRNAsubfilter.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(dat.genesubfilter,"dat.genesubfilter.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")





