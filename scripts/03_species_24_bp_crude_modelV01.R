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

cores=16

##Folder to save the outputs
output.folder="./results/"


# load data
dades=import("0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")

# regression functions
fastlm.fun <- function(x,y, z,x.var,y.var) {

  tryCatch({

    y <-y
    x <-scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    prev=(table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data),message = NA)

  }, warning = function(w) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA,message=paste("Warning:", w$message))

  }, error = function(e) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))

  })

}


#variables to assess
noms=grep("___",names(dades),value=T)

#outcome
yi="sbp_m_all"
#covariates
covari=c("age","gender","q005a","siteid_plate")

#remove individuals with NAs in the outcome
dades=dades[which(is.na(dades[,"q005a"])==F),]

#run the models
print("##### time #####")
system.time(res.lm1<-bplapply(noms,function(x){
  res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm1)
}, BPPARAM = MulticoreParam(cores)))
res.lm1 <- do.call(rbind, res.lm1)

#outcome
yi="dbp_m_all"

#run the models
system.time(res.lm2<-bplapply(noms,function(x){
  res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm2)
}, BPPARAM = MulticoreParam(cores)))
res.lm2 <- do.call(rbind, res.lm2)

#outcome
yi="syssd_all"

#run the models
system.time(res.lm5<-bplapply(noms,function(x){
  res.lm5 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm5)
}, BPPARAM = MulticoreParam(cores)))
res.lm5 <- do.call(rbind, res.lm5)

#outcome
yi="diasd_all"

#run the models
system.time(res.lm6<-bplapply(noms,function(x){
  res.lm6 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm6)
}, BPPARAM = MulticoreParam(cores)))
res.lm6 <- do.call(rbind, res.lm6)

#joint  phenotype analyses
res=data.frame(rbind(res.lm1,res.lm2,res.lm5,res.lm6),stringsAsFactors=F)
#estimate the fdr value
res$q.value <- p.adjust(res$p.value, method = "BH")

sbp=res[which(res$var.y=="sbp_m_all"),]
sbp1 <- sbp[order(sbp$p.value),]

dbp=res[which(res$var.y=="dbp_m_all"),]
dbp1 <- dbp[order(dbp$p.value),]


sbp_sd=res[which(res$var.y=="syssd_all"),]
sbp_sd1 <- sbp_sd[order(sbp_sd$p.value),]

dbp_sd=res[which(res$var.y=="diasd_all"),]
dbp_sd1 <- dbp_sd[order(dbp_sd$p.value),]


#save the result
write.table(sbp1,file=paste(output.folder,"/res_fastlm.24SBP_crude_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(dbp1,file=paste(output.folder,"/res_fastlm.24DBP_crude_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(sbp_sd1,file=paste(output.folder,"/res_fastlm.sysSD_crude_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(dbp_sd1,file=paste(output.folder,"/res_fastlm.diaSD_crude_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())


