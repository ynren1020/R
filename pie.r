#####plot scatter pie plot####
#####Nov 7,2018####

##load packages##
library(dplyr)
library(ggplot2)
library(scatterpie)
library(reshape2)
library(tidyr)

##load data in ##
pie <- read.delim2('RNA_sv_all_1000.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)   #variables are character#
##select column 5,6,7,12##
pie <-pie[,c(5,6,7,12)]

##split to get gene name##
for (i in 1:nrow(pie)){
pie$gene1[i]<-strsplit(pie$Break1_Gene[i],"_")[[1]][1]
pie$gene2[i]<-strsplit(pie$Break2_Gene[i],"_")[[1]][1]
}

##subset pie without INTERGENIC in gene1 AND gene2##
pie1<-filter(pie, !grepl('INTERGENIC',gene1,gene2))

##rearrange column##
pie2<-pie1[,c(5,6,3,4)]
##unique rows##
pie3<-unique(pie2) #nrow=15814
##wide to long:combine gene1 and gene2 into one column##
pie4<-melt(pie3,id.vars=c("SVTYPE","Samples"))
pie4<-pie4[,-3]
pie4<-rename(pie4,gene=value)
##unique rows for counting gene and sample##
pie5<-unique(pie4)
##separate samples in to several single sample##
#samples<-pie5$Samples
#for (i in 1:length(samples)){
#a[i]<-length(strsplit(samples[i],";")[[1]]) #17329
#}
#max(a) #91
##separate Samples into multiple rows##
pie6<-separate_rows(pie5,Samples,sep=";",convert = TRUE) #28014
pie7<-unique(pie6) #23223 unique genes and samples and SVTYPE
#length(unique(pie7$Samples)) #100
#nrow(unique(pie7[,c(2,3)])) #21349
##create one new variable sample:gene##
pie7$sample<-pie7$Samples
pie7$genes<-pie7$gene
pie7<-unite(pie7,"sample_gene",c("genes","sample"))
##final dataset to work on statistics##
pie8<-subset(pie7,!duplicated(pie7[,4]))

##summarize##
pie8 %>% group_by(gene) %>% summarize(count = n())
pie9<-pie8 %>% group_by(gene,SVTYPE) %>% summarize(count = n())
##long to wide##
pie10<-spread(pie9,SVTYPE,count)
pie10[is.na(pie10)] <- 0
pie11<-pie10 %>%group_by(gene)%>%summarize(total=sum(INV,TDUP,TRA))
pie_plot<-left_join(pie10,pie11,by="gene")
##create data file for plot##
##need gene_list.xlsx data sheet##
broken<-read.delim2('gene_list_broken.csv',sep = ',', skip = 1)
#nonlinear<-read.delim2('gene_list_nl.csv',sep = ',', skip = 1)
broken<-rename(broken,Broken.samples=X..of.samples) 

##join pie_plot and broken together##
pie_plot2<-left_join(pie_plot,broken,by="gene")
pie_plot2[is.na(pie_plot2)] <- 0
View(pie_plot2)
##plot##
pie_p<-ggplot() + geom_scatterpie(aes(x=Broken.samples, y=total, group=gene), data=pie_plot2,
                           cols=c("INV","TDUP","TRA")) + coord_equal()
ggsave("pie-1000.pdf",width=8,height = 8)

##use this##
pt44<-ggplot(pie,aes(x=Broken.samples,y=Nonlinear.samples,label=pie$labelg))+
  geom_point(color="grey70")+ 
  geom_vline(xintercept=10,linetype="dashed",size=0.1)+
  geom_hline(yintercept=20,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = c(0.85,0.85))

pt44pie<-pt44+geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene), data=pie[pie$labelg!='',],
                              cols=c("INV","TDUP","TRA")) + coord_equal()+labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")

ggsave("pie-1000-sub_final.pdf",width=8,height=8)


















