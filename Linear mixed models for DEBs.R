library(tidyverse)
library(reshape2)
library(car)
library(lme4)
library(lmerTest)
library(ggeffects)

# Load disordered eating behaviours
data <- read.csv("for demographic statistics.csv")
data$class <- as.factor(data$class)
data$sex <- as.factor(data$sex)

# Reshape data to long format
dieting_long <- melt(data, id.vars = c("ID", "class", "sex", "s1", "s2", "s3", "s4", "s5", "s6", "s7"),
                     measure.vars = c("dieting1", "dieting2", "dieting3", "dieting4"),
                     variable.name = "timepoint", value.name = "dieting")

binge_long <- melt(data, id.vars = c("ID", "class", "sex", "s1", "s2", "s3", "s4", "s5", "s6", "s7"),
                   measure.vars = c("binge1", "binge2", "binge3", "binge4"),
                   variable.name = "timepoint", value.name = "binge")

purging_long <- melt(data, id.vars = c("ID", "class", "sex", "s1", "s2", "s3", "s4", "s5", "s6", "s7"),
                     measure.vars = c("purging1", "purging2", "purging3", "purging4"),
                     variable.name = "timepoint", value.name = "purging")

# Run linear mixed models
dieting_model <- lmer(dieting ~ class * timepoint + sex + s1 + s2 + s3 + s4 + s5 + s6 + s7 + (1 | ID), 
                      data = dieting_long, REML = FALSE)
binge_model <- lmer(binge ~ class * timepoint + sex + s1 + s2 + s3 + s4 + s5 + s6 + s7 + (1 | ID), 
                    data = binge_long, REML = FALSE)
purging_model <- lmer(purging ~ class * timepoint + sex + s1 + s2 + s3 + s4 + s5 + s6 + s7 + (1 | ID), 
                      data = purging_long, REML = FALSE)

# Model summary
dieting_preds <- ggpredict(dieting_model, terms = c("timepoint", "class"))
binge_preds <- ggpredict(binge_model, terms = c("timepoint", "class"))
purging_preds <- ggpredict(purging_model, terms = c("timepoint", "class"))

# Plot dieting trajectories
dieting_preds$group <- factor(dieting_preds$group,
                              levels = c(0, 1),
                              labels = c("Persistently normal BMI", "Persistently Overweight"))

ggplot(dieting_preds, aes(x = as.numeric(x), y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("14", "16", "19", "23")) +
  labs(
    x = "Age",
    y = "",
    title = "Dieting behaviours",
    color = "Group",
    fill = "Group"
  ) +
  scale_color_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  scale_fill_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),  # axis line width
        axis.ticks = element_line(color = "black"), # tick width
        #axis.ticks.length = unit(0.25, "cm"),
        text = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

ggsave("dieting_trajectories.pdf", width = 4, height = 3)

# Plot binge eating trajectories
binge_preds$group <- factor(binge_preds$group,
                              levels = c(0, 1),
                              labels = c("Persistently normal BMI", "Persistently Overweight"))

ggplot(binge_preds, aes(x = as.numeric(x), y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("14", "16", "19", "23")) +
  labs(
    x = "Age",
    y = "",
    title = "Binge eating behaviours",
    color = "Group",
    fill = "Group"
  ) +
  scale_color_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  scale_fill_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),  # axis line width
        axis.ticks = element_line(color = "black"), # tick width
        #axis.ticks.length = unit(0.25, "cm"),
        text = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

ggsave("binge_trajectories.pdf", width = 4, height = 3)

# Plot purging trajectories
purging_preds$group <- factor(purging_preds$group,
                              levels = c(0, 1),
                              labels = c("Persistently normal BMI", "Persistently Overweight"))

ggplot(purging_preds, aes(x = as.numeric(x), y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("14", "16", "19", "23")) +
  labs(
    x = "Age",
    y = "",
    title = "Purging behaviours",
    color = "Group",
    fill = "Group"
  ) +
  scale_color_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  scale_fill_manual(values = c("Persistently normal BMI" = "grey", "Persistently Overweight" = "#f38e6e")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),  # axis line width
        axis.ticks = element_line(color = "black"), # tick width
        #axis.ticks.length = unit(0.25, "cm"),
        text = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

ggsave("purging_trajectories.pdf", width = 4, height = 3)
