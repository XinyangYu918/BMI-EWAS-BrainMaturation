setwd("MRS_lasso/")

#----------------------------------------
# Big LASSO regression: BMI residuals ~ DNAm
#----------------------------------------

# Libraries
library(data.table)
library(dplyr)
library(readr)
library(biglasso)

#----------------------------------------
# Load DNAm data (450K probes)
#----------------------------------------
# DNAm: rows = participants, columns = CpGs
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/DNAmethylation/for_analysis/age14/Data_14_for_BMI.RData")

# Select complete rows
complete_rows <- complete.cases(Data_14_for_BMI[, c("BMI","PCA_1","PCA_2","PCA_3","PCA_4",
                                                    "V1","V2","V3","V4",
                                                    "sex","sit1","sit2","sit3","sit4","sit5","sit6","sit7",
                                                    "PC1","PC2","wave1","wave2")])

# Subset data
dat <- Data_14_for_BMI[complete_rows, ]
# Extract DNAm beta values
DNAm <- as.data.frame(dat[, 2:372583])
# Set row names to participant IDs
rownames(DNAm) <- dat$ID

#----------------------------------------
# Load BMI residuals (from previous residualisation)
#----------------------------------------
resid_pheno <- read.table("BMI_residualised_OSCA.pheno")[,c(1,3)]
colnames(resid_pheno) <- c("ID", "BMI_resid")
rownames(resid_pheno) <- resid_pheno$ID

#----------------------------------------
# Run Big LASSO
#----------------------------------------
tmp.y <- resid_pheno[,'BMI_resid']

# Convert X to big.matrix
# Ensure all DNAm columns are numeric
DNAm <- DNAm %>% mutate(across(everything(), as.numeric))
X.bm <- as.big.matrix(DNAm[!is.na(tmp.y), ])
y.input <- tmp.y[!is.na(tmp.y)]

# Cross-validated LASSO (10-fold)
cvfit <- cv.biglasso(X.bm, y.input, seed = 1234, nfolds = 10, ncores = 1)

# Best lambda
lambda <- cvfit$lambda.min
print(paste0("Lambda selected: ", lambda))

# Extract non-zero coefficients
coefs <- coef(cvfit)
nonzero_idx <- which(coefs != 0)
out <- data.frame(CpG = rownames(coefs)[nonzero_idx], Weight = coefs[nonzero_idx])

# Save results
write.table(out, file = '/exports/eddie/scratch/s2112198/big_lasso_450K_BMI.txt',
            quote = F, row.names = F, col.names = T)
print("Big LASSO results saved!")


#### To adjust alpha ####
library(biglasso)
library(dplyr)

# X.bm = big.matrix of DNAm (samples x CpGs)
# y.input = BMI residuals (vector)

# Alpha values to test
alpha_vals <- c(0.25, 0.5, 0.75, 1)

# Store results
results <- data.frame(
  alpha = numeric(),
  lambda_min = numeric(),
  n_CpGs_selected = integer(),
  cvm_min = numeric()
)

# Loop over alpha values
for(alpha_val in alpha_vals) {
  
  cat("Running biglasso with alpha =", alpha_val, "...\n")
  
  # Cross-validated LASSO/Elastic Net
  cvfit <- cv.biglasso(X.bm, y.input, nfolds = 10, seed = 1234, ncores = 1, alpha = alpha_val)
  
  # Best lambda
  lambda_min <- cvfit$lambda.min
  
  # CV error at lambda_min
  cvm_min <- cvfit$cve[which.min(abs(cvfit$lambda - lambda_min))]
  
  # Extract coefficients
  coefs <- coef(cvfit)
  n_nonzero <- sum(coefs != 0)  # includes intercept
  
  # Store results
  results <- rbind(results, data.frame(
    alpha = alpha_val,
    lambda_min = lambda_min,
    n_CpGs_selected = n_nonzero - 1,  # remove intercept from count
    cvm_min = cvm_min
  ))
  
  # Optional: save top CpGs
  out <- data.frame(CpG = rownames(coefs)[which(coefs != 0)],
                    Weight = coefs[which(coefs != 0)])
  write.table(out, file = paste0("biglasso_BMI_alpha_", alpha_val, ".txt"),
              quote = F, row.names = F, col.names = T)
}

# Summary table
print(results)

# True outcome
y_true <- y.input
n <- length(y_true)

# Fold assignments from cv.biglasso
foldid <- cvfit$cv.ind   # same folds used by cv.biglasso
lambda_min <- cvfit$lambda.min  # penalty selected

# Initialise CV predictions
yhat_cv <- rep(NA, n)

# Loop through folds
for (k in sort(unique(foldid))) {
  cat("Processing fold:", k, "\n")
  
  train_idx <- which(foldid != k)
  test_idx  <- which(foldid == k)
  
  # Training subset as big.matrix
  X_train <- as.big.matrix(as.matrix(X.bm[train_idx, ]))
  
  # Fit model on training fold
  fit <- biglasso(X_train, y_true[train_idx], alpha = 1)
  
  # Test subset also as big.matrix
  X_test <- as.big.matrix(as.matrix(X.bm[test_idx, ]))
  
  # Predict on test fold
  yhat_cv[test_idx] <- predict(fit, X_test, type = "response", lambda = lambda_min)
}

# ---- Performance metrics ----
mse <- mean((y_true - yhat_cv)^2)
rss <- sum((y_true - yhat_cv)^2)
tss <- sum((y_true - mean(y_true))^2)
rsq <- 1 - rss/tss
cor_val <- cor(y_true, yhat_cv)

cat("Cross-validated Correlation:", cor_val, "\n")
cat("Cross-validated R²:", rsq, "\n")
cat("Cross-validated MSE:", mse, "\n")

# ---- Save cross-validated MRS ----
mrs_df <- data.frame(
  ID = rownames(DNAm),
  BMI_resid = y_true,
  MRS_cv = yhat_cv
)

write.csv(mrs_df, "bmiMRS_crossvalidated.csv", row.names = FALSE)

plot(y_true, yhat_cv, pch = 16, col = "steelblue",
     xlab = "Observed BMI residuals",
     ylab = "Predicted BMI residuals (cross-validated MRS)",
     main = paste0("bmiMRS: r = ", round(cor_val, 2),
                   ", R² = ", round(rsq, 3)))
abline(lm(yhat_cv ~ y_true), col = "red", lwd = 2)

