library(cdcanthro)
d <- data.table(sex=c(1,2,1,2,2), age=c(141,54,217,155,52), 
                wt=c(57,25,72,72,17.7), ht=c(143,102,166,169,105) )
d[,age := age+0.5] # because age was given as the completed number of months
d
d <- cdcanthro(d, age, wt, ht) # if BMI is not given, it's based on wt (kg) and ht (cm)
round(d,2)

data <- read.csv("CAT12_results/demographics_imputed.csv")
data$sex <- ifelse(data$sex == 0, 1, 2)
d1 <- data.frame(ID = data$ID,
                sex = data$sex, 
                class = data$class,
                age = data$age3,
                bmi = data$bmi2)
setDT(d1)
d1[,age := age/30.4375]
d1 <- cdcanthro(d1, age = age, bmi = bmi)
d1[, bmi_cat := fifelse(bmip < 5, "below_5th",
                        fifelse(bmip < 85, "5th_85th",
                                fifelse(bmip < 95, "85th_95th", "above_95th")))]
bmi_counts <- d1[, .N, by = .(class, bmi_cat)]
bmi_counts

data <- read.csv("for demographic statistics.csv")
# Convert to data.table (if not already)
setDT(data)

# Categorize BMI3 into standard WHO categories
data[, bmi3_cat := fifelse(
  bmi3 < 18.5, "underweight",
  fifelse(
    bmi3 < 25, "normal",
    fifelse(
      bmi3 < 30, "overweight",
      "obese"
    )
  )
)]

# Count number of individuals per class and BMI3 category
bmi3_counts <- data[, .N, by = .(class, bmi3_cat)]

# View results
bmi3_counts
