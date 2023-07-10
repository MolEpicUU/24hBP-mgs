print(Sys.time())

rm(list=ls())

library(rio)

info <- import("raw/HG3.A.7_tax.xlsx")
colnames(info)[1:3] <- c("id", "name", "level")

sbp=import("results/res_fastlm.24SBP_crude_model.tsv")
dbp=import("results/res_fastlm.24DBP_crude_model.tsv")
sbp_sd=import("results/res_fastlm.sysSD_crude_model.tsv")
dbp_sd=import("results/res_fastlm.diaSD_crude_model.tsv")


sbp$prev=format(round(sbp$prev,2),nsmall=2)
sbp$id=substring(sbp$var.x,nchar(sbp$var.x)-8, nchar(sbp$var.x))
sbp$estimate=format(round(sbp[,"estimate"],2),nsmall=2)
sbp$lower=format(round(sbp[,"lower"],2),nsmall=2)
sbp$upper=format(round(sbp[,"upper"],2),nsmall=2)


dbp$id=substring(dbp$var.x,nchar(dbp$var.x)-8, nchar(dbp$var.x))
dbp$estimate=format(round(dbp[,"estimate"],2),nsmall=2)
dbp$lower=format(round(dbp[,"lower"],2),nsmall=2)
dbp$upper=format(round(dbp[,"upper"],2),nsmall=2)

sbp_sd$id=substring(sbp_sd$var.x,nchar(sbp_sd$var.x)-8, nchar(sbp_sd$var.x))
sbp_sd$estimate=format(round(sbp_sd[,"estimate"],2),nsmall=2)
sbp_sd$lower=format(round(sbp_sd[,"lower"],2),nsmall=2)
sbp_sd$upper=format(round(sbp_sd[,"upper"],2),nsmall=2)

dbp_sd$id=substring(dbp_sd$var.x,nchar(dbp_sd$var.x)-8, nchar(dbp_sd$var.x))
dbp_sd$estimate=format(round(dbp_sd[,"estimate"],2),nsmall=2)
dbp_sd$lower=format(round(dbp_sd[,"lower"],2),nsmall=2)
dbp_sd$upper=format(round(dbp_sd[,"upper"],2),nsmall=2)

basic=merge(info[,c("name","id")],sbp[,c("id","prev","n","estimate","lower","upper","p.value","q.value")],by="id")
basic=merge(basic,dbp[,c("id","n","estimate","lower","upper","p.value","q.value")],by="id")
basic=merge(basic,sbp_sd[,c("id","n","estimate","lower","upper","p.value","q.value")],by="id")
basic=merge(basic,dbp_sd[,c("id","n","estimate","lower","upper","p.value","q.value")],by="id")

basic=basic[order(basic[,8]),]

names(basic)[names(basic)%in%grep("value",names(basic),value=T)]=paste(grep("value",names(basic),value=T),1:length(names(basic)[names(basic)%in%grep("value",names(basic),value=T)]),sep="")

for (i in grep("value",names(basic),value=T))
{
  basic[,i]=ifelse(basic[,i]<0.001,formatC(basic[,i],format="e",digits=2),round(basic[,i],digit=3))
}

names(basic)=c("id","name","prevalence","n.sbp","estimate.sbp","lower.sbp","upper.sbp","p.value.sbp","q.value.sbp",
"n.dbp","estimate.dbp","lower.dbp","upper.dbp","p.value.dbp","q.value.dbp",
"n.sbp_sd","estimate.sbp_sd","lower.sbp_sd","upper.sbp_sd","p.value.dbp_sd","q.value.dbp",
"n.dbp_sd","estimate.dbp_sd","lower.dbp_sd","upper.dbp_sd","p.value.dbp_sd","q.value.dbp")

write.table(basic,file="results/tables/Supplementary_table1.tsv",col.names=T,row.names=F,sep="\t")
print(Sys.time())

