rm(list=ls())

library(rio)

x=import("results/sensitivity.tsv")
x=x[,c(1,2,5,12,7:9,11,14,21,16:18,20,23,30,25:27,29,32)]
names(x)[3]="outcomes"

x$outcomes=ifelse(x$outcomes=="sbp_m_all","Mean 24-h systolic blood pressure",
                 ifelse(x$outcomes=="dbp_m_all","Mean 24-h diastolic blood pressure",
                        ifelse(x$outcomes=="syssd_all","Mean SD of systolic blood pressure",
                               ifelse(x$outcomes=="diasd_all","Mean SD of diastolic blood pressure",
                                      ifelse(x$outcomes=="bbps","Mean office systolic blood pressure","Mean office diastolic blood pressure")))))



x[,c(5:7,11:13,17:19)]=format(round(x[,c(5:7,11:13,17:19)],2),nsmall=2)


x=x[order(x$antib_p.value),]

for (i in grep("value",names(x),value=T))
{
  x[,i]=ifelse(x[,i]<0.001,formatC(x[,i],format="e",digits=2),round(x[,i],digit=3))
}

write.table(x,file="results/tables/Supplementary_table5.tsv",col.names=T,row.names=F,sep="\t")
