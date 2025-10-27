load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/DNAmethylation/for_analysis/age14/Data_14_for_BMI.RData")

complete_rows <- complete.cases(Data_14_for_BMI[, c("BMI","PCA_1","PCA_2","PCA_3","PCA_4",
                                                    "V1","V2","V3","V4",
                                                    "sex","sit1","sit2","sit3","sit4","sit5","sit6","sit7",
                                                    "PC1","PC2","wave1","wave2")])

Data_14_for_BMI <- Data_14_for_BMI[complete_rows, ]

library(sva)
Data_14_for_BMI$wave <- with(Data_14_for_BMI, 
                             ifelse(wave1 == 1, 1,
                                    ifelse(wave2 == 1, 2, 3)))

# Convert to factor
Data_14_for_BMI$wave <- factor(Data_14_for_BMI$wave, levels = c(1,2,3))
table(Data_14_for_BMI$wave)

# Define combat parameters
mod <- model.matrix(~ sex + PCA_1 + PCA_2 + PCA_3 + PCA_4 + V1+V2+V3+V4, data = Data_14_for_BMI)
beta_matrix <- as.matrix(t(Data_14_for_BMI[, 2:372583]))
dim(beta_matrix)

# Perform combat
beta_combat <- ComBat(dat = beta_matrix, batch = Data_14_for_BMI$wave, mod = mod)

# Rerun EWAS analysis
Methy  <- data.frame(t(beta_combat))
X <- model.matrix(~ BMI + PCA_1 + PCA_2 + PCA_3 + PCA_4 + V1 + V2 + V3 + V4 +
                    sex + sit1 + sit2 + sit3 + sit4 + sit5 + sit6 + sit7 +
                    PC1 + PC2 + wave1 + wave2,
                  data = Data_14_for_BMI[complete_rows, ])

Cov <- X[,-c(1,2)]
BMI <- data.frame(X[,2])

Num_Methy <- ncol(Methy)  
Num_Cov <- ncol(Cov) 
Num_BMI <- ncol(BMI) 

##Preparing the output files  
Origin_Beta <- matrix( data= NA, nrow=Num_Methy, ncol= Num_BMI, byrow=F, dimnames=NULL) 
colnames(Origin_Beta) <- colnames(BMI) 
rownames(Origin_Beta) <- colnames(Methy) 
Origin_SE <- matrix( data= NA, nrow=Num_Methy, ncol= Num_BMI, byrow=F, dimnames=NULL) 
colnames(Origin_SE) <- colnames(BMI) 
rownames(Origin_SE) <- colnames(Methy) 
Origin_P <- matrix( data= NA, nrow=Num_Methy, ncol= Num_BMI, byrow=F, dimnames=NULL) 
colnames(Origin_P) <- colnames(BMI) 
rownames(Origin_P) <- colnames(Methy) 

Covar <- matrix(unlist(Cov), ncol=ncol(Cov), byrow=F) 

for (i in 1:Num_BMI) { 
  for (j in 1:Num_Methy) { 
    print(i)
    Out <- summary(lm(BMI[,i]~Covar+Methy[,j])) 
    Origin_Beta[j,i] <- Out$coefficients[nrow(Out$coefficients),1] 
    Origin_SE[j,i] <- Out$coefficients[nrow(Out$coefficients),2] 
    Origin_P[j,i] <- Out$coefficients[nrow(Out$coefficients),4] 
  } 
} 
