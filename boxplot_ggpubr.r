##https://www.rdocumentation.org/packages/ggpubr/versions/0.2/topics/ggboxplot##
##http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/82-ggplot2-easy-way-to-change-graphical-parameters/##
##http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/##

p <- ggboxplot(df, x = "type", y = "scoreN",
               color = "type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "type",facet.by = "gene")+
               rotate_x_text(90)
p
# Change the plot orientation: horizontal
ggpar(p, orientation = "horiz")
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("Solid Tissue Normal", "Primary Tumor"), c("Primary Tumor", "Metastatic"), c("Solid Tissue Normal", "Metastatic") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 500)                   # Add global p-value
p1
# Violin plots with box plots inside
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change fill color by groups: dose
# add boxplot with white fill color
ggviolin(df, x = "type", y = "scoreN", fill = "type",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)                                      # Add global the p-value 
