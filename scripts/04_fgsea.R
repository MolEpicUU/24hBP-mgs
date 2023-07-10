rm(list=ls())
# load libraries and set options

library(rio)
library(fgsea)

set.seed(1)

# load data

data <- import("0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")
info <- import("raw/HG3.A.7_tax.xlsx")
gmm <- import("raw/MGS_HG3A.GMMs2MGS.RData")
gmm_info <- import("raw/GMM_reference.csv")

sbp=import("results/res_fastlm.24SBP_crude_model.tsv")
dbp=import("results/res_fastlm.24DBP_crude_model.tsv")
sbp_sd=import("results/res_fastlm.sysSD_crude_model.tsv")
dbp_sd=import("results/res_fastlm.diaSD_crude_model.tsv")

sbp$id=substring(sbp$var.x,nchar(sbp$var.x)-8, nchar(sbp$var.x))
dbp$id=substring(dbp$var.x,nchar(dbp$var.x)-8, nchar(dbp$var.x))
sbp_sd$id=substring(sbp_sd$var.x,nchar(sbp_sd$var.x)-8, nchar(sbp_sd$var.x))
dbp_sd$id=substring(dbp_sd$var.x,nchar(dbp_sd$var.x)-8, nchar(dbp_sd$var.x))

basic=data.frame(rbind(sbp,dbp,sbp_sd,dbp_sd),stringsAsFactors = F)

colnames(info)[1:3] <- c("id", "name", "level")
colnames(gmm_info) <- c("id", "name", "hl1", "hl2")

# enrichment function
gsea.fun <- function(xres, pathways0,pheno) {
  
  
  xxx=NULL
  for (z in pheno)
  {
    res=xres
    res=res[which(res$var.y%in%z),]
    stats <- -res$p.value
    names(stats) <- res$id
    
    pathways <- lapply(pathways0, function(pathway) {
      
      pathway[which(pathway %in% names(stats))]
      
    })
    
    gsea1 <- tryCatch({
      
      gsea1 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate >= 0)], na = "keep"), scoreType = "pos", eps = 0))
      data.frame(pathway = gsea1$pathway,var=z,sign="positive", estimate = gsea1$NES, p.value = gsea1$pval, size = gsea1$size, leading = sapply(gsea1$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)
      
    }, warning = function(w) {
      
      data.frame(pathway = names(pathways), var=z,sign="positive",estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)
      
    }, error = function(e) {
      
      data.frame(pathway = names(pathways), var=z,sign="positive",estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)
      
    })
    
    gsea2 <- tryCatch({
      
      gsea2 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate < 0)], na = "keep"), scoreType = "pos", eps = 0))
      data.frame(pathway = gsea2$pathway,var=z,sign="negative", estimate = gsea2$NES, p.value = gsea2$pval, size = gsea2$size, leading = sapply(gsea2$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)
      
    }, warning = function(w) {
      
      data.frame(pathway = names(pathways),var=z,sign="negative", estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)
      
    }, error = function(e) {
      
      data.frame(pathway = names(pathways),var=z,sign="negative", estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)
      
    })
    g=rbind(gsea1,gsea2)
    xxx=rbind(xxx,g)
  }
  
  xxx$q <- p.adjust(xxx$p.value, method = "BH", n = sum(!is.na(xxx$p.value)))
  x1=xxx[which(xxx$sign%in%"positive"==T),]
  x2=xxx[which(xxx$sign%in%"negative"==T),]
  xxx=merge(x1,x2, by=c("pathway","var"),all.x=T,all.y=T)
  xxx=xxx[,c(1,2,6,4,5,9,8,13,11,12,16,15)]
  names(xxx)=c("id","var","positive_size","positive_estimate","positive_p.value","positive_q.value","positive_message",
               "negative_size","negative_estimate","negative_p.value","negative_q.value","negative_message")
  xxx=xxx[order(xxx$positive_p.value),]
data.frame(xxx)

}

# gmm

gmm.res <- gsea.fun(basic, gmm,pheno=unique(names(table(basic$var.y))))
gmm.res <- data.frame(gmm.res, gmm_info[match(gmm.res$id, gmm_info$id), c("name", "hl1", "hl2")])
export(gmm.res, "results/gmm.tsv")

# run gsea on genera

genera <- lapply(unique(info$genus), function(genus) {
  
  info$id[which(info$genus == genus)]
  
})
names(genera) <- unique(info$genus)
genera <- genera[which(!names(genera) == "unclassified")]
genera.res <- gsea.fun(basic, genera,pheno=unique(names(table(basic$var.y))))

export(genera.res, "results/genera.tsv")


