library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(stringr)

#### load and format data ####
ct <- read.csv("CAT12_results/ROI_aparc_DK40_thickness.csv")
ct <- ct[,-c(3,4,11,12)]

list <- read.csv("CAT12_results/participant_id.csv")

list <- list %>%
  mutate(site = case_when(
    s1 == 1 ~ 1,
    s2 == 1 ~ 2,
    s3 == 1 ~ 3,
    s4 == 1 ~ 4,
    s5 == 1 ~ 5,
    s6 == 1 ~ 6,
    s7 == 1 ~ 7,
    TRUE ~ 8  # default case if none of the above is true
  ))

df <- subset(ct, ID %in% list$ID)
df <- merge(df, list, by = "ID")

#### for analysis between overewight and healthy controls ####
df$ParticipantID <- rep(1:713, each = 2)
results <- data.frame(ROI=colnames(df)[3:70])

for (c in 3:70)
{
  model <- lmer(df[,c] ~ class * time + (1|ParticipantID) + s1 + s2 + s3 + s4 + s5 + s6 + s7 + TIV + sex, REML = T, data = df)
  # Extract the summary
  summary_model <- summary(model)
  # Try to extract the t value and p-value for group:time interaction
  beta_value <- summary_model$coefficients[,"Estimate"]["class:time"]
  se_value <- summary_model$coefficients[,"Std. Error"]["class:time"]
  t_value <- summary_model$coefficients[,"t value"]["class:time"]
  p_value <- summary_model$coefficients[,"Pr(>|t|)"]["class:time"]
  # Store the t and p values
  results[c-2,2] <- beta_value
  results[c-2,3] <- se_value
  results[c-2,4] <- t_value
  results[c-2,5] <- p_value
}

names(results) <- c("ROI","beta","se","t","p")
results$fdr.p = p.adjust(results$p, method = "fdr")
results$bh.p = p.adjust(results$p, method = "bonferroni")

sum(results$fdr.p<0.05)
sum(results$bh.p<0.05)

write.csv(results, file = "CAT12_results/CT/results.csv", na = "", row.names = F)

#### for plot of group differences ####
resid <- data.frame(ID=df$ID)
df$ROI_sum <- rowSums(df[, colnames(df) %in% results$ROI[results$bh.p < 0.05]])
model <- lmer(ROI_sum ~ s1 + s2 + s3 + s4 + s5 + s6 + s7 + sex + (1|ParticipantID) + TIV, REML = T, data = df)
resid[,2] <- resid(model)

names(resid)[2] <- c("ROI_sum")

resid$time <- df$time
resid$class <- df$class

# Calculate the mean for each time and group combination
class_means <- resid %>%
  group_by(time, class) %>%
  summarise(mean_ROI_sum = mean(ROI_sum, na.rm = TRUE),
            se_ROI_sum = sd(ROI_sum, na.rm = TRUE) / sqrt(n()),
            lower_ci_ROI_sum = mean_ROI_sum - (1.96* se_ROI_sum),
            upper_ci_ROI_sum = mean_ROI_sum + (1.96* se_ROI_sum))

# Plot the interaction of time by class
class_means$class <- factor(class_means$class, 
                            levels = c(1, 2), 
                            labels = c("Persistently normal BMI", "Persistently overweight"))
class_means <- class_means %>%
  mutate(time = as.numeric(time)) %>%
  mutate(time = case_when(
    time == "1" ~ "14",
    time == "2" ~ "23"
  ))

# Plotting the mean predictions with confidence intervals
library(ggplot2)
library(stringr) # for str_wrap
library(forcats)
library(dplyr)

class_means <- class_means %>%
  mutate(class = fct_recode(class,
                            "Normal BMI" = "Persistently normal BMI",
                            "Overweight" = "Persistently overweight"
  ))

ggplot(class_means, aes(x = time, y = mean_ROI_sum, group = class)) +
  geom_line(aes(linetype = as.factor(class), color = as.factor(class)), size = 1.3) +
  geom_ribbon(aes(ymin = lower_ci_ROI_sum, ymax = upper_ci_ROI_sum, fill = as.factor(class), group = class), alpha = 0.2) +
  geom_point(aes(color = as.factor(class)), size = 3) +
  labs(title = "Differences in CT trajectories", x = "Age", y = "Sum of adjusted CT across ROIs") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.justification = "left",
    legend.box.just = "left",  
    legend.key.width = unit(1.2, "cm")
  ) +
  scale_linetype_manual(
    name = "Group",
    labels = str_wrap(c("Normal BMI", "Overweight"), width = 15),
    values = c("Normal BMI" = "solid", "Overweight" = "solid")
  ) +
  scale_shape_discrete(
    name = "Group",
    labels = str_wrap(c("Normal BMI", "Overweight"), width = 15)
  ) +
  scale_color_manual(
    name = "Group",
    labels = str_wrap(c("Normal BMI", "Overweight"), width = 15),
    values = c("Normal BMI" = "grey", "Overweight" = "orangered")
  ) +
  scale_fill_manual(
    name = "Group",
    labels = str_wrap(c("Normal BMI", "Overweight"), width = 15),
    values = c("Normal BMI" = "grey", "Overweight" = "orangered")
  )

ggsave("CAT12_results/CT/plot_new.pdf", width=4, height=4, dpi=300)
