rm(list=ls())
library(rio)
library(compareGroups)
library(rmarkdown)

data=import("0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")
data2=import("0_data/dades.pheno_malmo_2866_uppsala_395_validation.csv")

var=c("siteid","age","gender","q005a","cur_smoke","BMI","Fibrer","Energi_kcal","diab","diab_treat","HC_treat","ppi","antibio","UC","sbp_m_all","dbp_m_all","syssd_all","diasd_all","bbps","bbpd","sodium","crea","sodium_kawa")

data=data[,var]
data$siteid=paste(data$siteid,"_ABPM",sep="")
data2=data2[,var]
data2$siteid=paste(data2$siteid,"_non-ABPM",sep="")
data=data.frame(rbind(data,data2),stringsAsFactors=F)

data$cur_smoke=factor(data$cur_smoke,levels=c(0,1))
data$diab=factor(data$diab,levels=c(0,1))
data$diab_treat=factor(data$diab_treat,levels=c(0,1))
data$HC_treat=factor(data$HC_treat,levels=c(0,1))
data$ppi=factor(data$ppi,levels=c(0,1))
data$antibio=factor(data$antibio,levels=c(0,1))
data$UC=factor(data$UC,levels=c(0,1))

data$siteid=factor(data$siteid,levels=c("Uppsala_ABPM","Malmö_ABPM","Uppsala_non-ABPM","Malmö_non-ABPM"))
data$q005a=factor(data$q005a,levels=c("scandinavia","europe","asia","other"))

names(data)[which(names(data)%in%"sbp_m_all"==T)]="SBP_24h"
names(data)[which(names(data)%in%"dbp_m_all"==T)]="DBP_24h"
names(data)[which(names(data)%in%"syssd_all"==T)]="SBP_SD"
names(data)[which(names(data)%in%"diasd_all"==T)]="DBP_SD"
names(data)[which(names(data)%in%"bbps"==T)]="SBP_office"
names(data)[which(names(data)%in%"bbpd"==T)]="DBP_office"
names(data)[which(names(data)%in%"q005a"==T)]="country_birth"
names(data)[which(names(data)%in%"crea"==T)]="creatinine"

taula1=createTable(compareGroups(siteid~age+gender+country_birth+cur_smoke+BMI+Fibrer+Energi_kcal+diab+diab_treat+HC_treat+ppi+antibio+UC+SBP_24h+DBP_24h+SBP_SD+DBP_SD+SBP_office+DBP_office+sodium+creatinine+sodium_kawa,method=c(Fibrer=2,Energi_kcal=2,sodium=2,creatinine=2,sodium_kawa=2),data=data),show.p.overall=F)

export2csv(taula1, file='results/tables/table1.csv')

