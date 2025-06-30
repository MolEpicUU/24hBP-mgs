print(Sys.time())

rm(list=ls())
#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)
library(BiocParallel)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(sandwich)
library(lmtest)
library(nnet)

cores=16

##Folder to save the outputs
output.folder="./results/"


# load data
dades=import("0_data/dades.pheno_malmo_758_uppsala_2937_new.csv")


# regression functions
fastlm.fun <- function(x,y, z,x.var,y.var) {
  
  tryCatch({
    
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    prev=(table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))
    
    df<-dfbetas(fit)
    nb=which.is.max(abs(df[,2]))
    fit2 <- lm(y ~ ., data[-nb,])
    
    coef2 <- summary(fit2)$coefficients
    ci2 <- confint(fit2)
    
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data), estimate.infl = coef2["x", 1], lower.infl = ci2["x"
                                                                                                                            , 1], upper.infl = ci2["x", 2], se.infl = coef2["x", 2], p.value.infl = coef2["x", 4], n.infl = nrow(data[-nb,]),message = NA)
    
  }, warning = function(w) {
    
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n= NA,estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n= NA,estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA, message = paste("Error:", e$message))
    
  })
  
}


#variables to assess
noms=grep("^hMGS.",names(dades),value=T)

#outcome
yi="sbp_m_all"
#covariates
covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI")

#remove individuals with NAs in the outcome
dades=dades[which(is.na(dades[,"q005a"])==F),]

#run the models
print("##### time #####")
system.time(res.lm1<-bplapply(noms,function(x){
  res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm1)
}, BPPARAM = MulticoreParam(cores)))
res.lm1 <- do.call(rbind, res.lm1)
res.lm1 <- data.frame(res.lm1,stringsAsFactors=F)
res.lm1$q.value=p.adjust(res.lm1$p.value,"BH")
res.lm1$q.value.infl<-p.adjust(res.lm1$p.value.infl,method="BH")
res.lm1=res.lm1[order(res.lm1$p.value),]

#outcome
yi="dbp_m_all"

#run the models
system.time(res.lm2<-bplapply(noms,function(x){
  res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm2)
}, BPPARAM = MulticoreParam(cores)))
res.lm2 <- do.call(rbind, res.lm2)
res.lm2 <- data.frame(res.lm2,stringsAsFactors=F)
res.lm2$q.value=p.adjust(res.lm2$p.value,"BH")
res.lm2$q.value.infl<-p.adjust(res.lm2$p.value.infl,method="BH")
res.lm2=res.lm2[order(res.lm2$p.value),]



#outcome
yi="syssd_all"

#run the models
system.time(res.lm5<-bplapply(noms,function(x){
  res.lm5 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm5)
}, BPPARAM = MulticoreParam(cores)))
res.lm5 <- do.call(rbind, res.lm5)
res.lm5 <- data.frame(res.lm5,stringsAsFactors=F)
res.lm5$q.value=p.adjust(res.lm5$p.value,"BH")
res.lm5$q.value.infl<-p.adjust(res.lm5$p.value.infl,method="BH")
res.lm5=res.lm5[order(res.lm5$p.value),]



#outcome
yi="diasd_all"

#run the models
system.time(res.lm6<-bplapply(noms,function(x){
  res.lm6 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm6)
}, BPPARAM = MulticoreParam(cores)))
res.lm6 <- do.call(rbind, res.lm6)
res.lm6 <- data.frame(res.lm6,stringsAsFactors=F)
res.lm6$q.value=p.adjust(res.lm6$p.value,"BH")
res.lm6$q.value.infl<-p.adjust(res.lm6$p.value.infl,method="BH")
res.lm6=res.lm6[order(res.lm6$p.value),]


#save the result
write.table(res.lm1,file=paste(output.folder,"/res_fastlm.24SBP_multiBMI_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(res.lm2,file=paste(output.folder,"/res_fastlm.24DBP_multiBMI_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(res.lm5,file=paste(output.folder,"/res_fastlm.sysSD_multiBMI_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(res.lm6,file=paste(output.folder,"/res_fastlm.diaSD_multiBMI_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())


