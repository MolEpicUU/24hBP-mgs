print(Sys.time())

rm(list=ls())

library(rio)

x=import("results/genera.tsv")
x$var=ifelse(x$var=="sbp_m_all","Mean 24-h systolic blood pressure",
                 ifelse(x$var=="dbp_m_all","Mean 24-h diastolic blood pressure",
                        ifelse(x$var=="syssd_all","Mean SD of systolic blood pressure",
                               ifelse(x$var=="diasd_all","Mean SD of diastolic blood pressure",
                                      ifelse(x$var=="bbps","Mean office systolic blood pressure","Mean office diastolic blood pressure")))))

x=x[,c(1:6,8:11)]

for (i in grep("estimate",names(x),value=T))
{
  x[,i]=format(round(x[,i],2),nsmall=2)
}

x=x[order(x$positive_p.value),]


for (i in grep("value",names(x),value=T))
{
  x[,i]=ifelse(x[,i]<0.001,formatC(x[,i],format="e",digits=2),round(x[,i],digit=3))
}

names(x)[1:2]=c("genus","traits")

write.table(x,file="results/tables/Supplementary_table4.tsv",col.names=T,row.names=F,sep="\t")
print(Sys.time())


