####################2019-03-20#######################################################
##create https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30306-4 Figure 4C###
##33 cohort with 33 defined color in len_boxplot.py##################################
#####################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

exitron<-read.delim("TCGA.tumor-specific.exitron.txt",header=TRUE,stringsAsFactors = FALSE)

cohort<-c('ACC','BLCA','BRCA', 'CESC', 'CHOL','COAD','DLBC','ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',  'LUAD', 'LUSC',
          'MESO', 'OV', 'PAAD','PCPG', 'PRAD','READ','SARC','SKCM','STAD', 'TGCT', 'THCA', 'THYM','UCEC', 'UCS', 'UVM')
color<-c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
         '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82')
cohortcolor<-data_frame(cohort=cohort,color=color)

exitronplot<-left_join(exitron,cohortcolor,by="cohort")

##create rank of number and generate median of it within each cohort##
exitronplot<-exitronplot %>%
  group_by(cohort) %>%
  mutate(my_ranks = order(order(number, decreasing=FALSE)))%>%
  mutate(my_median = median(number))%>%
  mutate(my_mean_ranks=mean(my_ranks))
##order by median##
exitronplot<-as.data.frame(exitronplot)
exitronplot<-exitronplot[order(exitronplot$my_median),]
##create factor variable of cohort, rename color##
exitronplot$cohort_f<-factor(exitronplot$cohort,levels=unique(exitronplot$cohort))
exitronplot<-rename(exitronplot,my_color=color)

##plot##
p<-ggplot(exitronplot,aes(x=my_ranks, y=number,color=cohort_f))+
  geom_point(size=0.5)+
  scale_colour_manual(values=unique(exitronplot$my_color))+
  labs(y="Number of Exitron")


p1<-p+facet_grid(.~cohort_f,switch="x",scales = "free_x")+
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  geom_segment(aes(x=median(my_ranks)-100,xend=median(my_ranks)+100,y=my_median,yend=my_median),colour="red",size=1)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(0, "mm"),panel.border = element_blank(),                      
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside")

p1
#p2<-p1+scale_x_continuous(breaks=exitronplot$my_mean_ranks, labels=exitronplot$cohort)+
#  theme(axis.text.x=element_text(angle=90, vjust=.5))
#p2



