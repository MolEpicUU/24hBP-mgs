# Load necessary libraries
library(rio)
source("/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1_CHAMP/beta_diversity/dmanova.r")

# Load datasets
df1 <- import("/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1_CHAMP/0_data/dades.pheno_malmo_758_uppsala_2937_new_mgs1_clr.csv")
df2 <- import("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_bray_curtis_dissimilarities_v1.0.tsv")

# Prepare data for analysis
df1 <- df1[, c("subject_id", "bbps", "bbpd", "sbp_m_all", "dbp_m_all", "age", "gender", "q005a", "siteid_plate", "cur_smoke", "Fibrer", "Energi_kcal", "diab_treat", "HC_treat", "sodium_kawa", "bmi")]
df1 <- df1[complete.cases(df1), ]

rownames(df2) <- df2[, 1]
df2 <- df2[, -1]

# Identify common subjects in both datasets
common_ids <- intersect(df1$subject_id, colnames(df2))

# Filter df1 based on common subjects
df1 <- df1[match(common_ids, df1$subject), ]

# Subset Bray-Curtis matrix to keep only common subjects (rows & columns)
df2 <- df2[common_ids, common_ids]

# Run DMANOVA adjusting for covariates
cov_total <- dmanova.overall(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi, data = df1)
sbp_office_total <- dmanova.overall(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbps, data = df1)
dbp_office_total <- dmanova.overall(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbpd, data = df1)
sbp_total <- dmanova.overall(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbps + sbp_m_all, data = df1)
dbp_total <- dmanova.overall(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbpd + dbp_m_all, data = df1)

sbp_office <- dmanova(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbps, data = df1)
dbp_office <- dmanova(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbpd, data = df1)
sbp <- dmanova(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbps + sbp_m_all, data = df1)
dbp <- dmanova(as.matrix(df2) ~ age + gender + q005a + siteid_plate + cur_smoke + Fibrer + Energi_kcal + diab_treat + HC_treat + sodium_kawa + bmi + bbpd + dbp_m_all, data = df1)

cov_total <- as.data.frame(cov_total$aov.tab)[1, 5]
sbp_office_total <- as.data.frame(sbp_office_total$aov.tab)[1, 5]
dbp_office_total <- as.data.frame(dbp_office_total$aov.tab)[1, 5]
sbp_total <- as.data.frame(sbp_total$aov.tab)[1, 5]
dbp_total <- as.data.frame(dbp_total$aov.tab)[1, 5]

sbp_office <- as.data.frame(sbp_office$aov.tab)[1, 6]
dbp_office <- as.data.frame(dbp_office$aov.tab)[1, 6]
sbp <- as.data.frame(sbp$aov.tab)[1, 6]
dbp <- as.data.frame(dbp$aov.tab)[1, 6]

# Save summary results using your format
res <- data.frame(BloodPressure = c("Systolic", "Diastolic"), Model1R2 = c(cov_total, cov_total), Model2R2 = c(sbp_office_total, dbp_office_total), Model3R2 = c(sbp_total, dbp_total), DiffModel2Model1Pvalue = c(sbp_office, dbp_office), DiffModel3Model2Pvalue = c(sbp, dbp))
export(res, "/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1_CHAMP/beta_diversity/bloodpressure_dmanova_20250306.csv")
