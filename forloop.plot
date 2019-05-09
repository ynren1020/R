##QC project##
##for loop plot four type survival plot##
for (i in 13:26){
  status<-prad.surv.status[,i]

fit1<-survfit(Surv(DSS_time,DSS)~status,data=prad.surv.status)
res1<-ggsurvplot(fit1,data=prad.surv.status,
           xlab = "Days",
           ylab = "Disease Specific Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low"))
           )
ggsave(paste0(names(prad.surv.status[i]),".DS.png"), plot = print(res1), width = 8, height = 8, dpi = 500)

##overall survival##
fit2<-survfit(Surv(OS_time,OS)~status,data=prad.surv.status)
res2<-ggsurvplot(fit2,data=prad.surv.status,
                 xlab = "Days",
                 ylab = "Overall Survival Probability (%)",
                 conf.int=TRUE,
                 pval = TRUE,
                 fun="pct",
                 risk.table = TRUE,
                 size=1,
                 linetype = "strata",
                 legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))
ggsave(paste0(names(prad.surv.status[i]),".OS.png"), plot = print(res2), width = 8, height = 8, dpi = 500)

##disease free##
fit3<-survfit(Surv(DFI_time,DFI)~status,data=prad.surv.status)
res3<-ggsurvplot(fit3,data=prad.surv.status,
           xlab = "Days",
           ylab = "Disease Free Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))

ggsave(paste0(names(prad.surv.status[i]),".DF.png"), plot = print(res3), width = 8, height = 8, dpi = 500)
#AC092171.4,MCF2L-AS1,GATA2-AS1,
##progress free##
fit4<-survfit(Surv(PFI_time,PFI)~status,data=prad.surv.status)
res4<-ggsurvplot(fit4,data=prad.surv.status,
          xlab = "Days",
          ylab = "Progression Free Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))

ggsave(paste0(names(prad.surv.status[i]),".PF.png"), plot = print(res4), width = 8, height = 8, dpi = 500)
}
