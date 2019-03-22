###non-linear splicing###
###Nov 06 2018###

##load packages##
library(dplyr)
library(ggplot2)
library(ggrepel)

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
plot(dat$Broken.samples,dat$Nonlinear.samples) #x=20,y=20
abline(h = 10, v = 10, col = "gray60")
abline(h = 20, v = 20, col = "gray60")


#for (i in 1:nrow(dat)){
#dat$area[i]<-ifelse(dat$Broken.samples[i]>20,4,ifelse(dat$Nonlinear.samples[i]>25,2,ifelse(dat$Broken.samples[i]<=20&dat$Nonlinear.samples[i]<=25,3,1)))
#}
#for (i in 1:nrow(dat)){dat$plot[i]<-ifelse(dat$area[i]==4|dat$area[i]==2 & dat$all[i]==dat$gene[i], dat$plot[i]<-1,ifelse(dat$area[i]==4|dat$area[i]==2 & dat$prostate[i]==dat$gene[i],dat$plot[i]<-2,dat$plot[i]<-0))}

#for (i in 1:nrow(dat)){
#  dat$size[i]<-ifelse(dat$area[i]==4 | dat$area[i] == 2,ifelse(dat$all[i]==dat$gene[i]&dat$prostate[i]==dat$gene[i],dat$size[i]<-3,ifelse(dat$all[i]==dat$gene[i]&dat$prostate[i]!=dat$gene[i],dat$size[i]<-2,dat$size[i]<-1)),dat$size[i]<-NA)
#}

for (i in 1:nrow(dat)){
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelall[i]<-as.character(dat$all[i]),' ')
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
  ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>20,dat$plot[i]<-1,dat$plot[i]<-0)
}



for (i in 1:nrow(dat)){
 
  ifelse(dat$plot[i]==1&!is.na(dat$all[i]), dat$labelall[i]<-dat$all[i],dat$labelall[i]<-"")
  ifelse(dat$plot[i]==1&!is.na(dat$prostate[i]), dat$labelpro[i]<-dat$prostate[i],dat$labelpro[i]<-"")
}


##plot##
set.seed(00)
p<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=ifelse(plot==1,dat$all,'')))+
  geom_point(color="red")+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

pd<-p+geom_text_repel()

ggsave("all.pdf",width=8,height=8) 

########color#####
pt<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,color=as.factor(color)))+
  geom_point(size = 3)+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")+
  
ggsave("nonlabel.pdf",width=8,height=8)

ptl<-pt+geom_label(label=ifelse(dat$color==1,dat$all,""), color="black", size=3)
ggsave("label.pdf",width=8,height=8)

##label again##
for (i in 1:nrow(dat)){
  ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>20,dat$labelcolor[i]<-1,dat$labelcolor[i]<-0)
  #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
  
}

pt2<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,color=as.factor(labelcolor)))+
  geom_point(size = 1)+ 
  geom_vline(xintercept=20,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text(aes(label=dat$labelall),size=3,nudge_x = 0.75,nudge_y = 0.75)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")
  
  ggsave("labelall2.pdf",width=8,height=8)
  
######all using repel##
  pt22<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelall))+
    geom_point(size=ifelse(dat$labelcolor==0,2,3),color=ifelse(dat$labelcolor==0,"grey","red"),shape=ifelse(dat$labelall=="",1,19))+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text_repel()+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("new_repall_1000.pdf",width=8,height=8)
  
  
  ####prostate###
  pt3<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,color=as.factor(labelcolor)))+
    geom_point(size = 1)+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text(aes(label=dat$labelpro),size=3,nudge_x = 0.75,nudge_y = 0.75)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("labelpro.pdf",width=8,height=8)
########use ggrepel for prostate##
  pt33<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelpro))+
    geom_point(size=ifelse(dat$labelcolor==0,2,3),color=ifelse(dat$labelcolor==0,"grey","red"),shape=ifelse(dat$labelpro=="",1,19))+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text_repel()+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("new_reppro_1000.pdf",width=8,height=8)

  #################selected genes#########################
  
  slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")
  
  for (i in 1:nrow(dat)){
    #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelall[i]<-as.character(dat$all[i]),' ')
    #ifelse(dat$Broken.samples[i]>20 | dat$Nonlinear.samples[i]>25,dat$labelpro[i]<-as.character(dat$prostate[i]),' ')
    ifelse(dat$all[i]%in%slgene | dat$prostate[i]%in%slgene,dat$colorg[i]<-1,dat$colorg[i]<-0)
    ifelse(dat$colorg[i]==1,dat$labelg[i]<-dat$gene[i],dat$labelg[i]<-'')
  }
  dat$colorg<-as.factor(dat$colorg)
  
  pt4<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,color=as.factor(colorg)))+
    geom_point(size = 2)+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text(aes(label=dat$labelg),size=4,nudge_x = 1,nudge_y = 1)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("labelselected.pdf",width=8,height=8)
  
#######selected genes using ggrepel########
  
  pt44<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelg))+
    geom_point(size=ifelse(dat$colorg==0,2,5),color=ifelse(dat$colorg==0,"grey","red"))+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text_repel()+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("new_relselected_1000.pdf",width=8,height=8)
  
  ##overlap:color size shape#
  pt46<-ggplot(dat,aes(x=Broken.samples,y=Nonlinear.samples,label=dat$labelg))+
    geom_point(size=ifelse(dat$colorg==0,2,4),color=ifelse(dat$colorg==0,"grey","red"),shape=dat$colorg)+ 
    geom_vline(xintercept=20,linetype="dashed",size=0.1)+
    geom_hline(yintercept=20,linetype="dashed",size=0.1)+
    geom_text_repel()+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "none")
  
  ggsave("new_relselected_shape_1000.pdf",width=8,height=8)
  
  
