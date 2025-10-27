# Load data
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/DNAmethylation/for_analysis/age14/Data_14_for_BMI.RData")

# Select complete rows
complete_rows <- complete.cases(Data_14_for_BMI[, c("BMI","PCA_1","PCA_2","PCA_3","PCA_4",
                                                    "V1","V2","V3","V4",
                                                    "sex","sit1","sit2","sit3","sit4","sit5","sit6","sit7",
                                                    "PC1","PC2","wave1","wave2")])

# Subset data
dat <- Data_14_for_BMI[complete_rows, ]

# Build covariate model matrix
cov <- model.matrix(~ PCA_1 + PCA_2 + PCA_3 + PCA_4 + V1 + V2 + V3 + V4 +
                      sex + sit1 + sit2 + sit3 + sit4 + sit5 + sit6 + sit7 +
                      PC1 + PC2 + wave1 + wave2,
                    data = dat)

# Phenotype
BMI <- dat$BMI

# Fit linear model
fit <- lm(BMI ~ cov - 1)  # '-1' because cov already includes intercept

# Extract residuals (mean-centered)
BMI_resid <- resid(fit)

# OSCA-ready file: use ID column for both FID and IID
out_df <- data.frame(FID = dat$ID, IID = dat$ID, BMI_resid = BMI_resid)

# Write file
write.table(out_df, 
            file = "MRS_lasso/BMI_residualised_OSCA.pheno", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

cat("Wrote", nrow(out_df), "rows to BMI_residualised_OSCA.pheno\n")
