#####plot scatter pie plot function####
#####Jan 28,2019####

##load packages##
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(reshape2)
library(tidyr)

##FUNCTION 1##
readdat<-function(temp){
  ##load data in ##
  temp<-read.delim2(temp,sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  ##select column by name##
  temp<-temp[,c("Break1_Gene","Break2_Gene","SVTYPE","Samples")]
}

##FUNCTION 2##
##GET UNIQUE OBSERVATIONS##
modifygene<-function(temp){
  ##split to get gene name##
  for (i in 1:nrow(temp)){
    temp$gene1[i]<-strsplit(temp$Break1_Gene[i],"_")[[1]][1]
    temp$gene2[i]<-strsplit(temp$Break2_Gene[i],"_")[[1]][1]
  }
  ##subset temp without INTERGENIC in gene1 AND gene2##
  temp<-filter(temp, !grepl('INTERGENIC',gene1,gene2))
  ##rearrange column##
  temp<-temp[,c(5,6,3,4)]
  ##unique rows##
  temp<-unique(temp)
}



##FUNCTION 3 WIDE TO LONG UNIQUE GENE SVTYPE AND SAMPLE##

wide2long<-function(temp){
  ##wide to long:combine gene1 and gene2 into one column##
  temp<-melt(temp,id.vars=c("SVTYPE","Samples"))
  temp<-temp[,-3]
  temp<-rename(temp,gene=value)
  ##unique rows for counting gene and sample##
  temp<-unique(temp)
  ##separate Samples into multiple rows##IMPORTANT STEP##
  temp<-separate_rows(temp,Samples,sep=";",convert = TRUE) 
  temp<-unique(temp) #23223 unique genes and samples and SVTYPE
  
}

##FUNCTION 4##
##ONE SAMPLE TO ONE TYPE##
choose1type<-function(temp){
  ##create one new variable gene:sample##
  temp$sample<-temp$Samples
  temp$genes<-temp$gene
  temp<-unite(temp,"sample_gene",c("genes","sample"))
  ##final dataset to work on statistics##
  temp<-subset(temp,!duplicated(temp[,4]))
  
}


##FUNCTION 5##
##summarize/STATISTICS##
STcounts<-function(temp){
  temp<-temp %>% group_by(gene,SVTYPE) %>% summarize(count = n())
  ##long to wide##
  temp<-spread(temp,SVTYPE,count)
  temp[is.na(temp)] <- 0
  temp1<-temp %>%group_by(gene)%>%summarize(total=sum(INV,TDUP,TRA))
  temp_plot<-left_join(temp,temp1,by="gene")
}

#test if pieplotfunctionOne.R works##
source("pieplotfunctionOne.R")
pie5<-pieplotfunction('ucsf_su2c_common.txt')
######apply functions#########
pie1<-readdat('ucsf_su2c_common.txt')  ##first step output##
pie2<-modifygene(pie1)  ##OUTPUT FROM FUNCTION 2##
pie3<-wide2long(pie2) ##OUTPUT FROM FUNCTION 3##
pie4<-choose1type(pie3)  ##OUTPUT FROM FUNCTION 4##
pie5<-STcounts(pie4)

exon1<-readdat('with_exon_all.filtered.txt')  ##first step output##
exon2<-modifygene(exon1)  ##OUTPUT FROM FUNCTION 2##
exon3<-wide2long(exon2) ##OUTPUT FROM FUNCTION 3##
exon4<-choose1type(exon3)  ##OUTPUT FROM FUNCTION 4##
exon5<-STcounts(exon4)
##########


##FUNCTION 6##
broken<-function(file){
  xtemp<-read.delim2(file,sep = '\t',header=FALSE,stringsAsFactors = FALSE)
  xtemp<-separate(xtemp,V2,c("Broken.samples","V3"))
  xtemp<-rename(xtemp,gene=V1)
  xtemp$Broken.samples<-as.numeric(xtemp$Broken.samples)
  xtemp<-xtemp[,c("gene","Broken.samples")]
  return(xtemp)
  }

##read new broken data and transform to "gene, Broken.samples"##
broken_new<-broken('sort_by_count_new.txt')
##NO DEL NEW BROKEN DATA##
broken_new2<-broken('sort_by_count_new_no_del.txt')


##FUNCTION 7##
##JOIN DATA TO PLOT##
join4plot<-function(xfile,yfile){
  ##join broken and pie data for plot##
  jointemp<-full_join(xfile,yfile,by="gene")
  jointemp[is.na(jointemp)] <- 0
  jointemp<-rename(jointemp,Nonlinear.samples=total)
  
}

pie<-join4plot(broken_new2,pie5)


pie<-join4plot(broken_new,exon5)

##join broken and pie data for plot##
#pie<-full_join(broken_new2,pie5,by="gene")
#pie[is.na(pie)] <- 0
#pie<-rename(pie,Nonlinear.samples=total)
##add all cancer or prostate cancer gene note##
#prostate<-read.delim2('prostate.txt',sep = '\t',header = FALSE)
#all<-read.delim2('combined_gene_list.txt', sep = '\t')
#prostate$prostate<-prostate$V1
#all$all<-all$GENE
#prostate<-rename(prostate,gene=V1)
#all<-rename(all,gene=GENE)
#pie<-left_join(pie,all,by="gene")
#pie<-left_join(pie,prostate,by="gene")
#pie$all<-as.character(pie$all)
#pie$prostate<-as.character(pie$prostate)

##plot##
##define vertical and horizontal line intercept##
plot(pie$Broken.samples,pie$Nonlinear.samples) #x=10,y=20 #uscf 5 15
abline(h = 20, v = 10, col = "gray60")

######TRY THIS USE THIS####somepie#########
#################selected genes#########################
###1)only pie area I###
#slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")
slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FOXP1")
#####slgene in gene and in first area##use this for final plot##
##modify##try this##
#pie$colorg<-NULL
#for (i in 1:nrow(pie)){
#  pie$colorg[i]<-ifelse(pie$Broken.samples[i]>10&pie$Nonlinear.samples[i]<20,4,ifelse(pie$Nonlinear.samples[i]>20&pie$Broken.samples[i]<10,2,ifelse(pie$Broken.samples[i]<=10&pie$Nonlinear.samples[i]<=20,3,1)))
#  ifelse(pie$colorg[i]==1,pie$labelg[i]<-pie$gene[i],pie$labelg[i]<-'')
#}
#pie$colorg<-as.factor(pie$colorg)

##label gene in slgene## and area I USE THIS##
##area I II III IV##
pie$colorg<-NULL
pie$labelg<-NULL
for (i in 1:nrow(pie)){
  pie$colorg[i]<-ifelse(pie$Broken.samples[i]>=10&pie$Nonlinear.samples[i]<20,4,ifelse(pie$Nonlinear.samples[i]>=20&pie$Broken.samples[i]<10,2,ifelse(pie$Broken.samples[i]<10&pie$Nonlinear.samples[i]<20,3,1)))
  ifelse(pie$colorg[i]==1|pie$gene[i]%in%slgene,pie$labelg[i]<-pie$gene[i],pie$labelg[i]<-'')
}
pie$colorg<-as.factor(pie$colorg)
for (i in 1:nrow(pie)){
  
  ifelse(pie$colorg[i]==1|pie$gene[i]%in%slgene,pie$size[i]<-1,pie$size[i]<-0)
}




##plot regular scatter plot##
##overlap:color size shape#USE THIS ONE Final one##
pt46<-ggplot(pie,aes(x=Broken.samples,y=Nonlinear.samples,label=pie$labelg))+
  geom_point(size=ifelse(pie$colorg==1,3,2),color=ifelse(pie$colorg==1,"red",ifelse(pie$colorg==3,"grey90","grey")),shape=ifelse(pie$colorg==1,21,1),fill=ifelse(pie$colorg==1,"red",NA))+
  labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")+
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("ucsf_su2c_common2.pdf",width=8,height=8)


#####plot pie plot area I and slgenes###############
#######FINAL PIE PLOT USE THIS!################

pt44<-ggplot(pie,aes(x=Broken.samples,y=Nonlinear.samples,label=pie$labelg))+
  geom_point(color="grey70")+ 
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_text(aes(x=9, label="x=10", y=90), colour="blue", angle=90, text=element_text(size=11))+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text(aes(x=54, label="y=20", y=21), colour="blue", angle=0, text=element_text(size=11))+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = c(0.85,0.85))

pt44pie<-pt44+geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene,r=size), data=pie[pie$labelg!='',],
                              cols=c("INV","TDUP","TRA")) + coord_equal()+labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")

ggsave("testfunction.pdf",width=8,height=8)

write.table(pie,file="pie-with_exon_all.filtered.txt",row.names = FALSE,col.names = TRUE,quote = FALSE,sep="\t")
#############################################################






