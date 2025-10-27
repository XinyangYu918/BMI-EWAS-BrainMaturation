library(dplyr)

# Load STRATIFY methylation data
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/ENIGMA-Epigenetics2/STRATIFY/beta_Quantile_PSC2_STRATIFY.rda")

# Extract CpG columns
beta_cpgs <- colnames(beta)[-(1:3)]   # adjust indices if necessary

# Clean CpG names (remove suffixes like "_1" or "_A" etc.)
beta_cpgs_clean <- gsub("_.*$", "", beta_cpgs)

# Map original names to cleaned names
cpg_map <- data.frame(orig = beta_cpgs, clean = beta_cpgs_clean)

# Keep only CpGs that exist in LASSO output
cpg_map <- cpg_map %>% filter(clean %in% out$CpG)

# Subset LASSO weights to these CpGs
sumstats_sub <- out %>% filter(CpG %in% cpg_map$clean)

# Ensure matching order
cpg_map <- cpg_map[match(sumstats_sub$CpG, cpg_map$clean), ]

# Compute methylation risk score
# Extract only CpG columns from beta
beta_matrix <- as.matrix(beta[, cpg_map$orig])

# Extract corresponding LASSO weights
weights <- sumstats_sub$Weight

# Compute methylation score per individual
meth_score <- beta_matrix %*% weights

# Add MRS to beta dataframe
beta$MethScore <- as.numeric(meth_score)

# Subset to STRATIFY sample (exclude IMAGEN)
beta.sub <- beta %>% filter(Cohort != "IMAGEN")

# Merge with BMI phenotype
STRATIFY_BMI <- read.csv("../STRATIFY_BMI.csv")

# Ensure common subject ID column
colnames(beta.sub)[1] <- "Subject"

# Merge methylation scores with phenotype
STRATIFY_data <- merge(STRATIFY_BMI, beta.sub[, c("Subject", "MethScore")],
                       by = "Subject", all.x = TRUE)

# Linear models
# Full model: BMI ~ MRS + covariates
full_model <- lm(BMI ~ MethScore + PC1 + PC2 + V1 + V2 + V3 + V4 +
                   sex + site1 + site2 + caseADHD + caseMDD + X + 
                   caseAUD + casePsychosis, data = STRATIFY_data)

# Reduced model: BMI ~ covariates only
reduced_model <- lm(BMI ~ PC1 + PC2 + V1 + V2 + V3 + V4 +
                      sex + site1 + site2 + caseADHD + caseMDD + X +
                      caseAUD + casePsychosis, data = STRATIFY_data)

# Variance explained by MRS
r2_full <- summary(full_model)$r.squared
r2_reduced <- summary(reduced_model)$r.squared
r2_MRS <- r2_full - r2_reduced

cat("Variance explained by MRS (R²):", r2_MRS, "\n")

# Summary
summary(full_model)
