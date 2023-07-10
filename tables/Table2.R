print(Sys.time())

rm(list=ls())

library(rio)

res_val=import("results/alpha.lm_validation.csv")
res=import("results/alpha.lm.csv")
res=cbind(cbind(res[which(res$model%in%"crude"),],res[which(res$model%in%"multiple"),]),res[which(res$model%in%"BMI"),])
res$estimate_CI.1=paste(format(round(res[,3],2),nsmall=2)," (",format(round(res[,4],2),nsmall=2),", ",format(round(res[,5],2),nsmall=2),")",sep="")
res$estimate_CI.2=paste(format(round(res[,13],2),nsmall=2)," (",format(round(res[,14],2),nsmall=2),", ",format(round(res[,15],2),nsmall=2),")",sep="")
res$estimate_CI.3=paste(format(round(res[,23],2),nsmall=2)," (",format(round(res[,24],2),nsmall=2),", ",format(round(res[,25],2),nsmall=2),")",sep="")
res$cohort="ABPM cohort"
res=res[,c(34,1,2,8,31,18,32,28,33)]
names(res)[c(4,6,8)]=c("n.1","n.2","n.3")


res_val=cbind(cbind(res_val[which(res_val$model%in%"crude"),],res_val[which(res_val$model%in%"multiple"),]),res_val[which(res_val$model%in%"BMI"),])
res_val$estimate_CI.1=paste(format(round(res_val[,3],2),nsmall=2)," (",format(round(res_val[,4],2),nsmall=2),", ",format(round(res_val[,5],2),nsmall=2),")",sep="")
res_val$estimate_CI.2=paste(format(round(res_val[,13],2),nsmall=2)," (",format(round(res_val[,14],2),nsmall=2),", ",format(round(res_val[,15],2),nsmall=2),")",sep="")
res_val$estimate_CI.3=paste(format(round(res_val[,23],2),nsmall=2)," (",format(round(res_val[,24],2),nsmall=2),", ",format(round(res_val[,25],2),nsmall=2),")",sep="")
res_val$cohort="Non-ABPM cohort"
res_val=res_val[,c(34,1,2,8,31,18,32,28,33)]
res_val$var.y=c("bbps","bbpd")
names(res_val)[c(4,6,8)]=c("n.1","n.2","n.3")

res=as.data.frame(rbind(res,res_val),stringAsFactors=F)

res$var.y=ifelse(res$var.y=="sbp_m_all","Mean 24-h systolic blood pressure",
                 ifelse(res$var.y=="dbp_m_all","Mean 24-h diastolic blood pressure",
                        ifelse(res$var.y=="syssd_all","Mean SD of systolic blood pressure",
                               ifelse(res$var.y=="diasd_all","Mean SD of diastolic blood pressure",
                                      ifelse(res$var.y=="bbps","Mean office systolic blood pressure","Mean office diastolic blood pressure")))))

write.table(res,file="results/tables/Table2.tsv",col.names=T,row.names=F,sep="\t")

print(Sys.time())
