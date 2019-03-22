#####plot scatter pie plot function####
#####Nov 8,2018####

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

pie1<-readdat('RNA_sv_all_100.txt')  ##first step output##
pieF1<-readdat('RNA_sv_all_500.txt')

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

pie2<-modifygene(pie1)  ##OUTPUT FROM FUNCTION 2##
pieF2<-modifygene(pieF1)

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

pie3<-wide2long(pie2) ##OUTPUT FROM FUNCTION 3##
pieF3<-wide2long(pieF2)
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
pie4<-choose1type(pie3)  ##OUTPUT FROM FUNCTION 4##
pieF4<-choose1type(pieF3) 

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

pie5<-STcounts(pie4)
pieF5<-STcounts(pieF4)


##need gene_list.xlsx data sheet##
brokenp<-read.delim2('gene_list_broken.csv',sep = ',', skip = 1)
#nonlinear<-read.delim2('gene_list_nl.csv',sep = ',', skip = 1)
broken<-rename(broken,Broken.samples=X..of.samples) 
pie<-full_join(broken,pie5,by="gene")
pie[is.na(pie)] <- 0
pie<-rename(pie,Nonlinear.samples=total)

prostate<-read.delim2('prostate.txt',sep = '\t',header = FALSE)
all<-read.delim2('combined_gene_list.txt', sep = '\t')
prostate$prostate<-prostate$V1
all$all<-all$GENE
prostate<-rename(prostate,gene=V1)
all<-rename(all,gene=GENE)
pie<-full_join(pie,all,by="gene")
pie<-full_join(pie,prostate,by="gene")
pie$all<-as.character(pie$all)
pie$prostate<-as.character(pie$prostate)

##plot##
######TRY THIS USE THIS####somepie#########
#################selected genes#########################
###1)only pie area I###
#slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FAT1","FOXP1")
slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FOXP1")
#####slgene in gene and in first area##use this for final plot##
##modify##try this##
pie$colorg<-NULL
for (i in 1:nrow(pie)){
  pie$colorg[i]<-ifelse(pie$Broken.samples[i]>10&pie$Nonlinear.samples[i]<20,4,ifelse(pie$Nonlinear.samples[i]>20&pie$Broken.samples[i]<10,2,ifelse(pie$Broken.samples[i]<=10&pie$Nonlinear.samples[i]<=20,3,1)))
  ifelse(pie$colorg[i]==1,pie$labelg[i]<-pie$gene[i],pie$labelg[i]<-'')
}
pie$colorg<-as.factor(pie$colorg)

##label gene in slgene## and area I USE THIS##
pie$colorg<-NULL
for (i in 1:nrow(pie)){
  pie$colorg[i]<-ifelse(pie$Broken.samples[i]>10&pie$Nonlinear.samples[i]<20,4,ifelse(pie$Nonlinear.samples[i]>20&pie$Broken.samples[i]<10,2,ifelse(pie$Broken.samples[i]<=10&pie$Nonlinear.samples[i]<=20,3,1)))
  ifelse(pie$colorg[i]==1|pie$gene[i]%in%slgene,pie$labelg[i]<-pie$gene[i],pie$labelg[i]<-'')
}
pie$colorg<-as.factor(pie$colorg)
for (i in 1:nrow(pie)){
  
  ifelse(pie$colorg[i]==1|pie$gene[i]%in%slgene,pie$size[i]<-1,pie$size[i]<-0)
}


#####plot pie area I and slgenes###############
#######FINAL PIE PLOT USE THIS!################

pt44<-ggplot(pie,aes(x=Broken.samples,y=Nonlinear.samples,label=pie$labelg))+
  geom_point(color="grey70")+ 
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = c(0.85,0.85))

pt44pie<-pt44+geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene,r=size), data=pie[pie$labelg!='',],
                              cols=c("INV","TDUP","TRA")) + coord_equal()+labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")

ggsave("pie-100-sub_final-check.pdf",width=8,height=8)

write.table(pie,file="pie-100.txt",row.names = FALSE,col.names = TRUE,quote = FALSE,sep="\t")
#############################################################

###plot area I PIE#####
pt44<-ggplot(pie,aes(x=Broken.samples,y=Nonlinear.samples,label=pie$labelg))+
  geom_point()+ 
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

pt44pie500<-pt44+geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene), data=pie[pie$Broken.samples>=10&pie$Nonlinear.samples>=20,],
                              cols=c("INV","TDUP","TRA")) + coord_equal()+labs(y="Nonlinear-100")

ggsave("pie-500-sub.pdf",width=8,height=8)






#################all pie#################
pie_500<-ggplot() + geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene), data=pie,
                                    cols=c("INV","TDUP","TRA")) + coord_equal()+labs(y="Nonlinear-500")+
  labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")+
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave("pie-500-final.pdf",width=8,height=8)


###################100##########
pie_100<-ggplot() + geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene), data=pie,
                                  cols=c("INV","TDUP","TRA")) + coord_equal()+labs(y="Nonlinear-100")

ggsave("pie-100.pdf",width=8,height=8)




#####test#######



