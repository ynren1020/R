###non-linear splicing###
###Nov 09 2018###
###for len100###
###new broken sample###

##load packages##
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

##read data in##
##len1000:gene_list_nl.csv##
##len100:sort_by_count_100.txt##
##len500:sort_by_count_500.txt##
broken<-read.delim2('gene_list_broken.csv',sep = ',', skip = 1) #old broken sample#
nonlinear<-read.delim2('gene_list_nl.csv',sep = ',', skip = 1)
prostate<-read.delim2('prostate.txt',sep = '\t',header = FALSE)
all<-read.delim2('combined_gene_list.txt', sep = '\t')
##rename##
broken<-rename(broken,Broken.samples=X..of.samples)  #newname=oldname
nonlinear<-rename(nonlinear,Nonlinear.samples=X..of.samples)
prostate$prostate<-prostate$V1
all$all<-all$GENE
prostate<-rename(prostate,gene=V1)
all<-rename(all,gene=GENE)
##join two df into one##
dat<-full_join(broken,nonlinear,by="gene")
dat[is.na(dat)] <- 0
dat<-left_join(dat,all,by="gene")
dat<-left_join(dat,prostate,by="gene")
dat$all<-as.character(dat$all)
dat$prostate<-as.character(dat$prostate)

##define vertical and horizontal line intercept##
plot(dat$Broken.samples,dat$Nonlinear.samples) #x=20,y=40 (old,new) 
abline(h = 20, v = 10, col = "gray60")
abline(h = 20, v = 15, col = "blue")
abline(h = 30, v = 16, col = "red")  #x=10 y=20




#################selected genes#########################

#slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")
slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FOXP1","CHD9")


#####slgene in gene and in first area##use this for final plot##

##modify##try this##
dat$colorg<-NULL
for (i in 1:nrow(dat)){
  dat$colorg[i]<-ifelse(dat$Broken.samples[i]>10&dat$Nonlinear.samples[i]<20,4,ifelse(dat$Nonlinear.samples[i]>20&dat$Broken.samples[i]<10,2,ifelse(dat$Broken.samples[i]<=10&dat$Nonlinear.samples[i]<=20,3,1)))
  ifelse(dat$colorg[i]==1|dat$gene[i]%in%slgene,dat$labelg[i]<-dat$gene[i],dat$labelg[i]<-'')
}
dat$colorg<-as.factor(dat$colorg)


#######selected genes using ggrepel########



##overlap:color size shape#USE THIS ONE Final one##
pt46<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelg))+
  geom_point(size=ifelse(dat$colorg==1,3,2),color=ifelse(dat$colorg==1,"red",ifelse(dat$colorg==3,"grey90","grey")),shape=ifelse(dat$colorg==1,21,1),fill=ifelse(dat$colorg==1,"red",NA))+
  labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")+
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("old_relselected_shape_1000_final_final.pdf",width=8,height=8)


