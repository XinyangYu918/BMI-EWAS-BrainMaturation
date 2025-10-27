library(tidyverse)
library(ggplot2)
library(RColorBrewer)

#### For 2-class solution plot ####

# Load the data
data <- read.table("/Users/yuxinyang/Documents/temp/BMI_LPA_2class.dat", header=FALSE)
colnames(data) <- c("bmi1","bmi2","bmi3","sex","s1","s2","s3","s4","s5","s6","s7","CPROB1","CPROB2","C","ID")

# Replace '*' with NA in BMI columns
bmi_cols <- c("bmi1","bmi2","bmi3")
data[bmi_cols] <- lapply(data[bmi_cols], function(x) {
  x[x == "*"] <- NA          # replace * with NA
  as.numeric(x)              # convert to numeric
})

# Assign class based on highest probability if C not available
if(!"C" %in% colnames(data)){
  data$Class <- apply(data[, c("CPROB1","CPROB2")], 1, which.max)
} else {
  data$Class <- data$C
}

# Select BMI columns
bmi_cols <- c("bmi1","bmi2","bmi3")

# Reshape data for plotting
df_long <- data %>%
  select(Class, all_of(bmi_cols)) %>%
  pivot_longer(cols = all_of(bmi_cols), names_to = "Timepoint", values_to = "BMI") %>%
  mutate(Age = case_when(
    Timepoint == "bmi1" ~ 14,
    Timepoint == "bmi2" ~ 19,
    Timepoint == "bmi3" ~ 23
  ))

# Compute mean and SD per class
df_summary <- df_long %>%
  group_by(Class, Age) %>%
  summarise(Mean_BMI = mean(BMI, na.rm=TRUE),
            SD_BMI = sd(BMI, na.rm=TRUE),
            .groups="drop")

# Define legend labels
class_labels <- c("1" = "Class 1 (N = 1096; 84.57%)",
                  "2" = "Class 2 (N = 200; 15.43%)")

# Choose colors from RColorBrewer
line_colors <- brewer.pal(3, "Set2")[1:2]  # pick first 2 colors
fill_colors <- alpha(line_colors, 0.2)     # semi-transparent for ribbons

# Plot BMI trajectories with SD shading
ggplot(df_summary, aes(x=Age, y=Mean_BMI, color=factor(Class), group=Class)) +
  geom_line(size=1.5) +
  geom_point(size=4) +
  geom_ribbon(aes(ymin=Mean_BMI-SD_BMI, ymax=Mean_BMI+SD_BMI, fill=factor(Class)),
              alpha=0.2, color=NA) +
  scale_color_manual(values = line_colors, labels = class_labels) +
  scale_fill_manual(values = fill_colors, labels = class_labels) +
  scale_x_continuous(breaks=c(14,19,23)) +
  labs(title="BMI Trajectories by LPA Classes 
(2-class solution)",
       x="Age (year)",
       y="BMI",
       color="Class",
       fill="Class") +
  theme_bw(base_size=12) +
  theme(
    legend.position="right",
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    plot.title = element_text(hjust=0.5, face="bold")
  )

ggsave("/Users/yuxinyang/Documents/temp/BMI_LPA_2class_plot.pdf", width=6, height=5, dpi=300)

#### For 3-class solution plot ####
data <- read.table("/Users/yuxinyang/Documents/temp/BMI_LPA_3class.dat", header=FALSE)

# Assign column names (adjust if needed)
colnames(data) <- c("bmi1","bmi2","bmi3","sex","s1","s2","s3","s4","s5","s6","s7",
                    "CPROB1","CPROB2","CPROB3","C","ID")

bmi_cols <- c("bmi1","bmi2","bmi3")
data[bmi_cols] <- lapply(data[bmi_cols], function(x) {
  x[x == "*"] <- NA
  as.numeric(x)
})

if(!"C" %in% colnames(data)){
  data$Class <- apply(data[, c("CPROB1","CPROB2","CPROB3")], 1, which.max)
} else {
  data$Class <- data$C
}

df_long <- data %>%
  select(Class, all_of(bmi_cols)) %>%
  pivot_longer(cols = all_of(bmi_cols), names_to = "Timepoint", values_to = "BMI") %>%
  mutate(Age = case_when(
    Timepoint == "bmi1" ~ 14,
    Timepoint == "bmi2" ~ 19,
    Timepoint == "bmi3" ~ 23
  ))

df_summary <- df_long %>%
  group_by(Class, Age) %>%
  summarise(
    Mean_BMI = mean(BMI, na.rm=TRUE),
    SD_BMI = sd(BMI, na.rm=TRUE),
    .groups="drop"
  )

class_labels <- c(
  "1" = "Class 1 (N = 304; 23.46%)",
  "2" = "Class 2 (N = 940; 72.53%)",
  "3" = "Class 3 (N = 52; 4.01%)"
)

line_colors <- c("1" = "#8DA0CB", 
                 "2" = "#66C2A5", 
                 "3" = "#FC8D62") 

fill_colors <- alpha(line_colors, 0.2)  # semi-transparent for ribbons

p <- ggplot(df_summary, aes(x=Age, y=Mean_BMI, color=factor(Class), group=Class)) +
  geom_line(size=1.5) +
  geom_point(size=4) +
  geom_ribbon(aes(ymin=Mean_BMI-SD_BMI, ymax=Mean_BMI+SD_BMI, fill=factor(Class)),
              color=NA, alpha=0.2) +
  scale_color_manual(values = line_colors, labels = class_labels) +
  scale_fill_manual(values = fill_colors, labels = class_labels) +
  scale_x_continuous(breaks=c(14,19,23), limits=c(14,24)) +
  labs(title="BMI Trajectories by LPA Classes
(3-class solution)",
       x="Age (years)",
       y="BMI",
       color="Class",
       fill="Class") +
  theme_bw(base_size=12) +
  theme(
    legend.position="right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=0.5, face="bold")
  )

ggsave("/Users/yuxinyang/Documents/temp/BMI_LPA_3class_plot.pdf", plot=p,
       width=6, height=5, dpi=300)

#### For 4-class solution plot ####
data <- read.table("/Users/yuxinyang/Documents/temp/BMI_LPA_4class.dat", header=FALSE)

colnames(data) <- c("bmi1","bmi2","bmi3","sex","s1","s2","s3","s4","s5","s6","s7",
                    "CPROB1","CPROB2","CPROB3","CPROB4","C","ID")

bmi_cols <- c("bmi1","bmi2","bmi3")
data[bmi_cols] <- lapply(data[bmi_cols], function(x) {
  x[x == "*"] <- NA
  as.numeric(x)
})

if(!"C" %in% colnames(data)){
  data$Class <- apply(data[, c("CPROB1","CPROB2","CPROB3","CPROB4")], 1, which.max)
} else {
  data$Class <- data$C
}

df_long <- data %>%
  select(Class, all_of(bmi_cols)) %>%
  pivot_longer(cols = all_of(bmi_cols), names_to = "Timepoint", values_to = "BMI") %>%
  mutate(Age = case_when(
    Timepoint == "bmi1" ~ 14,
    Timepoint == "bmi2" ~ 19,
    Timepoint == "bmi3" ~ 23
  ))

df_summary <- df_long %>%
  group_by(Class, Age) %>%
  summarise(
    Mean_BMI = mean(BMI, na.rm=TRUE),
    SD_BMI = sd(BMI, na.rm=TRUE),
    .groups="drop"
  )

class_labels <- c(
  "1" = "Class 1 (N = 609; 46.99%)",
  "2" = "Class 2 (N = 501; 38.66%)",
  "3" = "Class 3 (N = 37; 2.85%)",
  "4" = "Class 4 (N = 149; 11.50%)"
)

line_colors <- c("#66C2A5", "#8DA0CB", "#FC8D62", "#E78AC3")

fill_colors <- alpha(line_colors, 0.2)      # semi-transparent for ribbons

p <- ggplot(df_summary, aes(x=Age, y=Mean_BMI, color=factor(Class), group=Class)) +
  geom_line(size=1.5) +
  geom_point(size=4) +
  geom_ribbon(aes(ymin=Mean_BMI-SD_BMI, ymax=Mean_BMI+SD_BMI, fill=factor(Class)),
              color=NA, alpha=0.2) +
  scale_color_manual(values = line_colors, labels = class_labels) +
  scale_fill_manual(values = fill_colors, labels = class_labels) +
  scale_x_continuous(breaks=c(14,19,23), limits=c(14,24)) +
  labs(title="BMI Trajectories by LPA Classes 
(4-class solution)",
       x="Age (years)",
       y="BMI",
       color="Class",
       fill="Class") +
  theme_bw(base_size=12) +
  theme(
    legend.position="right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=0.5, face="bold")
  )

ggsave("/Users/yuxinyang/Documents/temp/BMI_LPA_4class_plot.pdf", plot=p,
       width=6, height=5, dpi=300)

