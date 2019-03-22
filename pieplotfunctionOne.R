##pie plot functions united into ONE##
##Jan28,2019##

##FUNCTION FOR Y##
pieplotfunction<-function(temp){

##FUNCTION 1##
readdat<-function(temp){
  ##load data in ##
  temp<-read.delim2(temp,sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  ##select column by name##
  temp<-temp[,c("Break1_Gene","Break2_Gene","SVTYPE","Samples")]
  }
temp1<-readdat(temp)
##FUNCTION 2##
##GET UNIQUE OBSERVATIONS##
modifygene<-function(temp1){
  ##split to get gene name##
  for (i in 1:nrow(temp1)){
    temp1$gene1[i]<-strsplit(temp1$Break1_Gene[i],"_")[[1]][1]
    temp1$gene2[i]<-strsplit(temp1$Break2_Gene[i],"_")[[1]][1]
  }
  ##subset temp without INTERGENIC in gene1 AND gene2##
  temp1<-filter(temp1, !grepl('INTERGENIC',gene1,gene2))
  ##rearrange column##
  temp1<-temp1[,c(5,6,3,4)]
  ##unique rows##
  temp1<-unique(temp1)
  
}
temp2<-modifygene(temp1)


##FUNCTION 3 WIDE TO LONG UNIQUE GENE SVTYPE AND SAMPLE##

wide2long<-function(temp2){
  ##wide to long:combine gene1 and gene2 into one column##
  temp2<-melt(temp2,id.vars=c("SVTYPE","Samples"))
  temp2<-temp2[,-3]
  temp2<-rename(temp2,gene=value)
  ##unique rows for counting gene and sample##
  temp2<-unique(temp2)
  ##separate Samples into multiple rows##IMPORTANT STEP##
  temp2<-separate_rows(temp2,Samples,sep=";",convert = TRUE) 
  temp2<-unique(temp2) #23223 unique genes and samples and SVTYPE
}
temp3<-wide2long(temp2)

##FUNCTION 4##
##ONE SAMPLE TO ONE TYPE##
choose1type<-function(temp3){
  ##create one new variable gene:sample##
  temp3$sample<-temp3$Samples
  temp3$genes<-temp3$gene
  temp3<-unite(temp3,"sample_gene",c("genes","sample"))
  ##final dataset to work on statistics##
  temp3<-subset(temp3,!duplicated(temp3[,4]))
  
}
temp4<-choose1type(temp3)

##FUNCTION 5##
##summarize/STATISTICS##
STcounts<-function(temp4){
  temp4<-temp4 %>% group_by(gene,SVTYPE) %>% summarize(count = n())
  ##long to wide##
  temp4<-spread(temp4,SVTYPE,count)
  temp4[is.na(temp4)] <- 0
  temp5<-temp4 %>%group_by(gene)%>%summarize(total=sum(INV,TDUP,TRA))
  temp_plot<-left_join(temp4,temp5,by="gene")
}
final<-STcounts(temp4)
}

##FUNCTION for X##
broken<-function(file){
  xtemp<-read.delim2(file,sep = '\t',header=FALSE,stringsAsFactors = FALSE)
  xtemp<-separate(xtemp,V2,c("Broken.samples","V3"))
  xtemp<-rename(xtemp,gene=V1)
  xtemp$Broken.samples<-as.numeric(xtemp$Broken.samples)
  xtemp<-xtemp[,c("gene","Broken.samples")]
  return(xtemp)
}

##FUNCTION TO JOIN DATA (X AND Y) TO PLOT##
join4plot<-function(xfile,yfile){
  ##join broken and pie data for plot##
  jointemp<-full_join(xfile,yfile,by="gene")
  jointemp[is.na(jointemp)] <- 0
  jointemp<-rename(jointemp,Nonlinear.samples=total)
  
}

##FUNCTION OF PLOT TO FIND LINE INTERCEPT##
testplot<-function(file,X,Y){
##define vertical and horizontal line intercept##
plot(file$Broken.samples,file$Nonlinear.samples) #x=10,y=20 #uscf 5 15
abline(h = Y, v = X, col = "gray60")
}

##FUNCTION TO CREATE LABEL GENES##
labelfunction<-function(slgene,file,X,Y){
##label gene in slgene## and area I USE THIS##
##area I II III IV##
#slgene<-c("ACPP","MYC","PTEN","CDKN1B","RB1","FOXA1","TP53","AR","TMPRSS2","ERG","ELK4","ETV1","FOXP1")
file$colorg<-NULL
file$labelg<-NULL
for (i in 1:nrow(file)){
  file$colorg[i]<-ifelse(file$Broken.samples[i]>=X&file$Nonlinear.samples[i]<Y,4,ifelse(file$Nonlinear.samples[i]>=Y&file$Broken.samples[i]<X,2,ifelse(file$Broken.samples[i]<X&file$Nonlinear.samples[i]<Y,3,1)))
  ifelse(file$colorg[i]==1|file$gene[i]%in%slgene,file$labelg[i]<-file$gene[i],file$labelg[i]<-'')
}
file$colorg<-as.factor(file$colorg)
for (i in 1:nrow(file)){
  
  ifelse(file$colorg[i]==1|file$gene[i]%in%slgene,file$size[i]<-1,file$size[i]<-0)
}
return(file)
}

##FUNCTION FOR REGULAR SCATTER PLOT##
regscatter<-function(file,X,Y,name){
##plot regular scatter plot##
##overlap:color size shape#USE THIS ONE Final one##
pt46<-ggplot(file,aes(x=Broken.samples,y=Nonlinear.samples,label=file$labelg))+
  geom_point(size=ifelse(file$colorg==1,3,2),color=ifelse(file$colorg==1,"red",ifelse(file$colorg==3,"grey90","grey")),shape=ifelse(file$colorg==1,21,1),fill=ifelse(file$colorg==1,"red",NA))+
  labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")+
  geom_vline(xintercept=X,linetype="dashed",size=0.1)+
  geom_hline(yintercept=Y,linetype="dashed",size=0.1)+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = "none")

ggsave(paste0(name,".pdf"),width=8,height=8)

}


##FUNCTION TO PLOT PIE SCATTER PLOT##
piescatter<-function(file,X,Y,name){
#####plot pie plot area I and slgenes###############
#######FINAL PIE PLOT USE THIS!################

pt44<-ggplot(file,aes(x=Broken.samples,y=Nonlinear.samples,label=file$labelg))+
  geom_point(color="grey70")+ 
  geom_vline(xintercept=X,linetype="dashed",size=0.1)+
  #geom_text(aes(x=9, label="x=10", y=90), colour="blue", angle=90, text=element_text(size=11))+
  geom_hline(yintercept=Y,linetype="dashed",size=0.1)+
  #geom_text(aes(x=54, label="y=20", y=21), colour="blue", angle=0, text=element_text(size=11))+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(legend.position = c(0.85,0.85))

pt44file<-pt44+geom_scatterpie(aes(x=Broken.samples, y=Nonlinear.samples, group=gene,r=size), data=file[file$labelg!='',],
                              cols=c("INV","TDUP","TRA")) + coord_equal()+labs(x="# of Samples from SV (DNA)",y="# of Samples from Nonlinear Splicing (RNA)")

ggsave(paste0(name,".pie.pdf"),width=8,height=8)

}