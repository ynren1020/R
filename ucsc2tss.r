##read su2c_refseq.ucsc file in##
library(dplyr)
library(tidyr)

su2c <- read.delim("hg38_genelnRNA_sorted.ucsc",header=FALSE,sep = "\t")
#View(su2c)
##if there is no col.name##
#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
su2c<-su2c%>%rename("#bin"="V1","name"="V2","chrom"="V3","strand"="V4","txStart"="V5","txEnd"="V6")


for (i in 1:nrow(su2c)) {
  ifelse(su2c$strand[i]=="-",su2c$tssstart[i]<-su2c$txEnd[i]-500,su2c$tssstart[i]<-su2c$txStart[i]-500)
  ifelse(su2c$strand[i]=="-",su2c$tssend[i]<-su2c$txEnd[i]+500,su2c$tssend[i]<-su2c$txStart[i]+500)
}

su2c$trans<-rep("transcript",nrow(su2c))
su2c$length<-rep(1000,nrow(su2c))
su2c$dot<-rep(".",nrow(su2c))

su2csub<-su2c[,c("chrom","name","trans","tssstart","tssend","length","strand","dot","name")]
#View(su2csub)

write.table(su2csub,"hg38_genelnRNA_sorted.tss",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)