rm(list=ls())
# load libraries and set options

library(rio)
library(fgsea)

set.seed(1)

# load data

data <- import("0_data/dades.pheno_malmo_758_uppsala_2937_new.csv")
info <- import("raw/scapis_metagenomics_mgs_annotations_v1.0.tsv")
gmm2 <- import("raw/scapis_metagenomics_gmm_modules_v1.0.tsv")
gmm_info <- gmm2[,1:4]

sbp=import("results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp=import("results/res_fastlm.24DBP_multiBMI_model.tsv")
sbp_sd=import("results/res_fastlm.sysSD_multiBMI_model.tsv")
dbp_sd=import("results/res_fastlm.diaSD_multiBMI_model.tsv")

sbp$id=sbp$var.x
dbp$id=dbp$var.x
sbp_sd$id=sbp_sd$var.x
dbp_sd$id=dbp_sd$var.x

basic=data.frame(rbind(sbp,dbp,sbp_sd,dbp_sd),stringsAsFactors = F)

colnames(info)[1:3] <- c("id", "name", "level")
colnames(gmm_info) <- c("id", "name", "hl1", "hl2")



gmm_list <- setNames(
  lapply(strsplit(gmm2$mgs_id, ","), trimws), # Split mgs_id by commas and trim whitespace
  gmm2$gmm_id # Use gmm_id as the names of the list
)


# enrichment function
gsea.fun <- function(xres, pathways0) {

    res=xres
    stats <- -res$p.value
    names(stats) <- res$id

    pathways <- lapply(pathways0, function(pathway) {

      pathway[which(pathway %in% names(stats))]
   })

    gsea1 <- tryCatch({

      gsea1 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate >= 0)], na = "keep"), scoreType = "pos", eps = 0))
      data.frame(pathway = gsea1$pathway,sign="positive", estimate = gsea1$NES, p.value = gsea1$pval, size = gsea1$size, leading = sapply(gsea1$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)

    }, warning = function(w) {

      data.frame(pathway = names(pathways), sign="positive",estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)

    }, error = function(e) {

      data.frame(pathway = names(pathways), sign="positive",estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)

    })

    gsea2 <- tryCatch({

      gsea2 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate < 0)], na = "keep"), scoreType = "pos", eps = 0))
      data.frame(pathway = gsea2$pathway,sign="negative", estimate = gsea2$NES, p.value = gsea2$pval, size = gsea2$size, leading = sapply(gsea2$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)

    }, warning = function(w) {

      data.frame(pathway = names(pathways),sign="negative", estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)

    }, error = function(e) {

      data.frame(pathway = names(pathways),sign="negative", estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)

    })
    xxx=rbind(gsea1,gsea2)

  

  xxx$q <- p.adjust(xxx$p.value, method = "BH", n = sum(!is.na(xxx$p.value)))
  x1=xxx[which(xxx$sign%in%"positive"==T),]
  x2=xxx[which(xxx$sign%in%"negative"==T),]
  xxx=merge(x1,x2, by=c("pathway"),all.x=T,all.y=T)
  xxx=xxx[,c(1,5,3,4,8,7,12,10,11,15,14)]
  names(xxx)=c("id","positive_size","positive_estimate","positive_p.value","positive_q.value","positive_message",
               "negative_size","negative_estimate","negative_p.value","negative_q.value","negative_message")
  xxx=xxx[order(xxx$positive_p.value),]
data.frame(xxx)

}

# gmm

gmm.res.sbp <- gsea.fun(sbp, gmm_list)
gmm.res.dbp <- gsea.fun(dbp, gmm_list)
gmm.res.sbp_sd <- gsea.fun(sbp_sd, gmm_list)
gmm.res.dbp_sd <- gsea.fun(dbp_sd, gmm_list)

gmm.res=merge(gmm.res.sbp,gmm.res.dbp,by="id",all.x=T,all.y=T)
gmm.res=merge(gmm.res,gmm.res.sbp_sd,by="id",all.x=T,all.y=T)
gmm.res=merge(gmm.res,gmm.res.dbp_sd,by="id",all.x=T,all.y=T)

gmm.res <- data.frame(gmm.res, gmm_info[match(gmm.res$id, gmm_info$id), c("name", "hl1", "hl2")])
gmm.res<- gmm.res[,c(1,42:44,2:41)]


names(gmm.res)<- c("gmm_id","name","hl1","hl2",
	"positive_size.sbp","positive_estimate.sbp","positive_p.value.sbp","positive_q.value.sbp","positivie_message.sbp","negative_size.sbp","negative_estimate.sbp","negative_p.value.sbp","negative_q.value.sbp","negative_message.sbp",
	"positive_size.dbp","positive_estimate.dbp","positive_p.value.dbp","positive_q.value.dbp","positivie_message.dbp","negative_size.dbp","negative_estimate.dbp","negative_p.value.dbp","negative_q.value.dbp","negative_message.dbp",
	"positive_size.sbp_sd","positive_estimate.sbp_sd","positive_p.value.sbp_sd","positive_q.value.sbp_sd","positivie_message.sbp_sd","negative_size.sbp_sd","negative_estimate.sbp_sd","negative_p.value.sbp_sd","negative_q.value.sbp_sd","negative_message.sbp_sd",
	"positive_size.dbp_sd","positive_estimate.dbp_sd","positive_p.value.dbp_sd","positive_q.value.dbp_sd","positivie_message.dbp_sd","negative_size.dbp_sd","negative_estimate.dbp_sd","negative_p.value.dbp_sd","negative_q.value.dbp_sd","negative_message.dbp_sd")

gmm.res=gmm.res[order(gmm.res$positive_p.value.sbp),]

export(gmm.res, "results/gmm.tsv")

# run gsea on genera

genera <- lapply(unique(info$genus), function(genus) {

  info$id[which(info$genus == genus)]

})
names(genera) <- unique(info$genus)
genera <- genera[which(!names(genera) == "unclassified")]

genera.res.sbp <- gsea.fun(sbp, genera)
genera.res.dbp <- gsea.fun(dbp, genera)
genera.res.sbp_sd <- gsea.fun(sbp_sd, genera)
genera.res.dbp_sd <- gsea.fun(dbp_sd, genera)

genera.res=merge(genera.res.sbp,genera.res.dbp,by="id",all.x=T,all.y=T)
genera.res=merge(genera.res,genera.res.sbp_sd,by="id",all.x=T,all.y=T)
genera.res=merge(genera.res,genera.res.dbp_sd,by="id",all.x=T,all.y=T)

names(genera.res)<- c("genera",
	"positive_size.sbp","positive_estimate.sbp","positive_p.value.sbp","positive_q.value.sbp","positivie_message.sbp","negative_size.sbp","negative_estimate.sbp","negative_p.value.sbp","negative_q.value.sbp","negative_message.sbp",
	"positive_size.dbp","positive_estimate.dbp","positive_p.value.dbp","positive_q.value.dbp","positivie_message.dbp","negative_size.dbp","negative_estimate.dbp","negative_p.value.dbp","negative_q.value.dbp","negative_message.dbp",
	"positive_size.sbp_sd","positive_estimate.sbp_sd","positive_p.value.sbp_sd","positive_q.value.sbp_sd","positivie_message.sbp_sd","negative_size.sbp_sd","negative_estimate.sbp_sd","negative_p.value.sbp_sd","negative_q.value.sbp_sd","negative_message.sbp_sd",
	"positive_size.dbp_sd","positive_estimate.dbp_sd","positive_p.value.dbp_sd","positive_q.value.dbp_sd","positivie_message.dbp_sd","negative_size.dbp_sd","negative_estimate.dbp_sd","negative_p.value.dbp_sd","negative_q.value.dbp_sd","negative_message.dbp_sd")

genera.res=genera.res[order(genera.res$positive_p.value.sbp),]

export(genera.res, "results/genera.tsv")

