# load libraries
library(rio)
library(ppcor)

# import data
key <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Raw/Keys/id_conversion.txt")
met <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metabolomics/Processed/scapis_metabolon_measurements_batchnorm_v2.1.tsv")
met_labels <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metabolomics/Processed/scapis_metabolon_labels_batchnorm_v2.1.tsv")
met_anno <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metabolomics/Processed/scapis_metabolon_annotations_batchnorm_v2.1.tsv")
met_meta <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metabolomics/Processed/scapis_metabolon_metadata_batchnorm_v2.1.tsv")
mgs <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv")
mgs_key <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_key_file_v1.0.tsv")
mgs_anno <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_annotations_v1.0.tsv")
mgs_batch <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_technical_variables_v1.0.tsv")
pheno_m <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Phenotypes/Phenotypes/scapis_data_20190215_malmo.txt")
pheno_u <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Phenotypes/Phenotypes/scapis_data_20190215_uppsala.txt")
pheno2 <- rio::import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/preliminary/pheno_MGS_shannon_bray_curtis_MGP_4839_upp_4980_malmo.tsv")

bp_1 <- rio::import("/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1_CHAMP/0_data/dades.pheno_malmo_758_uppsala_2937_new.csv")
bp_2 <- rio::import("/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1_CHAMP/0_data/dades.pheno_malmo_2575_uppsala_195_validation.csv")

# clean & match data
colnames(pheno_m)[1] <- "id"
colnames(pheno_u)[1] <- "id"
pheno_m$scapis_id <- pheno_m$id
pheno_u$scapis_id <- key$subject_id[match(pheno_u$id, key$export_id)]
pheno <- rbind(pheno_m, pheno_u)

pheno2 <- pheno2[which(!is.na(pheno2$q005a)), ]
pheno2$scapis_id <- pheno$scapis_id[match(pheno2$id, pheno$id)]

mgs_key <- mgs_key[which(mgs_key$selected), ]

id <- pheno$scapis_id[which(pheno$scapis_id %in% pheno2$scapis_id & pheno$scapis_id %in% met$scapis_id & pheno$scapis_id %in% mgs$scapis_id & pheno$scapis_id %in% c(bp_1$subject_id, bp_2$subject_id))]
pheno <- pheno[match(id, pheno$scapis_id), ]
pheno2 <- pheno2[match(id, pheno2$scapis_id), ]
met <- met[match(id, met$scapis_id), ]
met_labels <- met_labels[match(id, met_labels$scapis_id), ]
met_meta <- met_meta[match(id, met_meta$scapis_id), ]
mgs <- mgs[match(id, mgs$scapis_id), ]
mgs_batch <- mgs_batch[match(id, mgs_batch$scapis_id), ]

met <- met[, -1]
met_labels <- met_labels[, -1]
mgs <- mgs[, -1]

met_labels <- apply(met_labels, 2, function(x) {
    x[which(!x %in% c("measured", "imputed"))] <- NA
    x
})

# filter and transform data
filter <- which(colSums(met_labels == "measured", na.rm = TRUE) >= 100)
met <- met[, filter]
met_labels <- met_labels[, filter]

colnames <- colnames(met)
met <- lapply(colnames(met), function(x) {
    met <- met[, x]
    labels <- met_labels[, x]
    if (mean(labels == "measured", na.rm = TRUE) < 0.02) met <- ifelse(labels == "measured", 1, 0)
    if (x %in% met_anno$met_id[which(grepl("Drug", met_anno$SUB_PATHWAY))] & mean(labels == "measured", na.rm = TRUE) < 0.99) met <- ifelse(labels == "measured", 1, 0)
    met
})
met <- do.call(cbind, met)
colnames(met) <- colnames

clr <- function(data) {
    zeroes <- data == 0
    data <- data + min(data[data > 0])
    data <- log(data) - rowMeans(log(data))
    data[zeroes] <- -Inf
    data
}
mgs <- clr(mgs)

# subset to 4 significant species
species <- c("Dysosmobacter sp001916835", "Intestinimonas massiliensis", "Streptococcus sp001556435", "Dysosmobacter sp001916835", "ER4 sp900317525")
mgs <- mgs[, which(colnames(mgs) %in% mgs_anno$mgs_id[which(mgs_anno$species %in% species)])]

# prepare covariates
scapis_id <- pheno$scapis_id
age <- pheno$agev1
sex <- ifelse(pheno$Gender == 1, "Male", "Female")
bmi <- pheno$bmi
batch <- paste(pheno2$q005a, mgs_batch$extraction_plate, met_meta$delivery_batch, sep = "_") # ethnicity + metagenomic extraction plate + metabolomic delivery batch
pheno <- data.frame(age, sex, bmi, batch)

# partial spearman model
spearman.fun <- function(y, x, z) {
    res <- BiocParallel::bplapply(colnames(x), function(x2) {
        tryCatch(
            {
                data <- data.frame(y, x = x[, x2], z)
                data <- data[complete.cases(data), ]
                data <- data[, apply(data, 2, function(x) length(unique(x)) > 1)]
                data <- as.data.frame(model.matrix(~., data)[, -1])
                warnings <- c()
                withCallingHandlers(res <- ppcor::pcor.test(data$y, data$x, data[, which(!colnames(data) %in% c("y", "x"))], method = "spearman"), warning = function(x) warnings <<- c(warnings, x$message))
                data.frame(met_id = x2, rho = res$estimate, pval = res$p.value, n = res$n, message = paste(warnings, collapse = "; "))
            },
            error = function(e) {
                data.frame(met_id = x2, rho = NA, pval = NA, n = nrow(data), message = paste(e$message, collapse = "; "))
            }
        )
    }, BPPARAM = BiocParallel::MulticoreParam(16))
    do.call(rbind, res)
}

# all results
res <- lapply(colnames(mgs), function(x) {
    res <- spearman.fun(mgs[, x], met, pheno)
    data.frame(mgs_id = x, res)
})
res <- do.call(rbind, res)
rio::export(res, "/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily.tsv")

# all
res <- rio::import("/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily.tsv")
res$species <- mgs_anno$species[match(res$mgs_id, mgs_anno$mgs_id)]
res$metabolite <- met_anno$CHEMICAL_NAME[match(res$met_id, met_anno$met_id)]
res$class <- met_anno$SUPER_PATHWAY[match(res$met_id, met_anno$met_id)]
res$subclass <- met_anno$SUB_PATHWAY[match(res$met_id, met_anno$met_id)]
res <- res[, c("mgs_id", "species", "met_id", "metabolite", "class", "subclass", "rho", "pval", "n")]
res <- res[order(res$pval), ]
res <- lapply(unique(res$mgs_id), function(x) {
    res <- res[which(res$mgs_id == x), ]
})
res <- do.call(rbind, res)
rio::export(res, "/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily_all.tsv")

# top 3
res <- rio::import("/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily.tsv")
res$species <- mgs_anno$species[match(res$mgs_id, mgs_anno$mgs_id)]
res$metabolite <- met_anno$CHEMICAL_NAME[match(res$met_id, met_anno$met_id)]
res$class <- met_anno$SUPER_PATHWAY[match(res$met_id, met_anno$met_id)]
res$subclass <- met_anno$SUB_PATHWAY[match(res$met_id, met_anno$met_id)]
res <- res[, c("mgs_id", "species", "met_id", "metabolite", "class", "subclass", "rho", "pval", "n")]
res <- res[order(res$pval), ]
res <- lapply(unique(res$mgs_id), function(x) {
    res <- res[which(res$mgs_id == x), ][1:3, ]
})
res <- do.call(rbind, res)
rio::export(res, "/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily_top3.tsv")

# top 3 characterized
res <- rio::import("/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily.tsv")
res$species <- mgs_anno$species[match(res$mgs_id, mgs_anno$mgs_id)]
res$metabolite <- met_anno$CHEMICAL_NAME[match(res$met_id, met_anno$met_id)]
res$class <- met_anno$SUPER_PATHWAY[match(res$met_id, met_anno$met_id)]
res$subclass <- met_anno$SUB_PATHWAY[match(res$met_id, met_anno$met_id)]
res <- res[, c("mgs_id", "species", "met_id", "metabolite", "class", "subclass", "rho", "pval", "n")]
res <- res[order(res$pval), ]
res <- res[which(!grepl("X - ", res$metabolite)), ]
res <- lapply(unique(res$mgs_id), function(x) {
    res <- res[which(res$mgs_id == x), ][1:3, ]
})
res <- do.call(rbind, res)
rio::export(res, "/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/species_metabolites_emily_top3_characterized.tsv")

