#####plot scatter pie plot function####
#####Jan 28,2019####

##load packages##
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(reshape2)
library(tidyr)

#test if pieplotfunctionOne.R works##
##produce y##
source("pieplotfunctionOne.R")
pie5<-pieplotfunction('ucsf_su2c_common.txt')
exon5<-pieplotfunction('with_exon_all.filtered.txt')
##read new broken data and transform to "gene, Broken.samples"##produce x##
broken_new<-broken('sort_by_count_new.txt')
##NO DEL NEW BROKEN DATA##
broken_new2<-broken('sort_by_count_new_no_del.txt')
##JOIN X AND Y##
pie<-join4plot(broken_new2,pie5)
exon<-join4plot(broken_new,exon5)

##define vertical and horizontal line intercept##
testplot(pie,10,20)
#################selected genes#########################
###1)only pie area I###
#slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")
slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FOXP1")
#####slgene in gene and in first area##use this for final plot##
pie.plot<-labelfunction(slgene,pie,10,20)
##plot regular scatter plot##
regscatter(pie.plot,10,20,'testreg')
#####plot pie plot area I and slgenes###############
piescatter(pie.plot,10,20,'test')

write.table(pie,file="pie-with_exon_all.filtered.txt",row.names = FALSE,col.names = TRUE,quote = FALSE,sep="\t")
#############################################################






