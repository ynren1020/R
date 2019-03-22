###non-linear splicing###
###Nov 08 2018###
###for len500###

##load packages##
library(dplyr)
library(ggplot2)
library(ggrepel)

##read data in##
##len1000:gene_list_nl.csv##
##len100:sort_by_count_100.txt##
##len500:sort_by_count_500.txt##
#broken<-read.delim2('gene_list_broken.csv',sep = ',', skip = 1) #old broken sample#
broken<-read.delim2('broken_new_list.txt',sep = '\t', header = FALSE,stringsAsFactors = FALSE) #new broken sample#
broken <-separate(broken,V2,c("Broken.samples","de"),sep = " ") ##split V2 into 2columns
broken<-broken[,c(1,2)]
broken$Broken.samples<-as.numeric(broken$Broken.samples) ##VERY IMPORTANT TO PAY ATTENTION TO DATA TYPE "CHARACTER TO NUMERIC"#
nonlinear<-read.delim2('sort_by_count_500.txt',sep = '\t', header = FALSE)
prostate<-read.delim2('prostate.txt',sep = '\t',header = FALSE)
all<-read.delim2('combined_gene_list.txt', sep = '\t')
##rename##
broken<-rename(broken,gene=V1)  #newname=oldname
nonlinear<-rename(nonlinear,Nonlinear.samples=V2)
nonlinear<-rename(nonlinear,gene=V1)
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
plot(dat$Broken.samples,dat$Nonlinear.samples) #x=20,y=40


for (i in 1:nrow(dat)){
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelall[i]<-as.character(dat$all[i]),' ')
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
  ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>40,dat$plot[i]<-1,dat$plot[i]<-0)
}



for (i in 1:nrow(dat)){
  
  ifelse(dat$plot[i]==1&!is.na(dat$all[i]), dat$labelall[i]<-dat$all[i],dat$labelall[i]<-"")
  ifelse(dat$plot[i]==1&!is.na(dat$prostate[i]), dat$labelpro[i]<-dat$prostate[i],dat$labelpro[i]<-"")
}

for (i in 1:nrow(dat)){
  ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>40,dat$labelcolor[i]<-1,dat$labelcolor[i]<-0)
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
  
}


######use this for label####
######all using repel#######
pt22<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelall))+
  geom_point(size=ifelse(dat$labelcolor==0,2,3),color=ifelse(dat$labelcolor==0,"grey","red"),shape=ifelse(dat$labelall=="",1,19))+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=40,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("new_repall_500.pdf",width=8,height=8)


####prostate###
########use ggrepel for prostate##
pt33<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelpro))+
  geom_point(size=ifelse(dat$labelcolor==0,2,3),color=ifelse(dat$labelcolor==0,"grey","red"),shape=ifelse(dat$labelpro=="",1,19))+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=40,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("new_reppro_500.pdf",width=8,height=8)

#################selected genes#########################

slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")

for (i in 1:nrow(dat)){
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelall[i]<-as.character(dat$all[i]),' ')
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
  ifelse(dat$all[i]%in%slgene | dat$prostate[i]%in%slgene,dat$colorg[i]<-1,dat$colorg[i]<-0)
  ifelse(dat$colorg[i]==1,dat$labelg[i]<-dat$gene[i],dat$labelg[i]<-'')
}
dat$colorg<-as.factor(dat$colorg)

#######selected genes using ggrepel########

pt44<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelg))+
  geom_point(size=ifelse(dat$colorg==0,2,5),color=ifelse(dat$colorg==0,"grey","red"))+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=40,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("new_relselected_500.pdf",width=8,height=8)

##overlap:color size shape#
pt46<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelg))+
  geom_point(size=ifelse(dat$colorg==0,2,4),color=ifelse(dat$colorg==0,"grey","red"),shape=dat$colorg)+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=40,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("new_relselected_shape_500.pdf",width=8,height=8)


