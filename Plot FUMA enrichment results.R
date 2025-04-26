library(ggplot2)
library(dplyr)
library(stringr)

## Load the top 100 CpGs associated with BMI at age 14
enrichment_results <- read.csv("top100_fuma.csv")

# Reorder geneset by adjP and remove "_" from the terms
enrichment_results$GeneSet <- enrichment_results$GeneSet %>%
  str_replace_all("_", " ")

enrichment_results$Category <- factor(enrichment_results$Category, 
                                      levels = c("Canonical_Pathways", "Cell_type_signature", "Hallmark_gene_sets", "Immunologic_signatures"))
enrichment_results <- enrichment_results %>%
  arrange(Category, adjP)
enrichment_results$GeneSet <- stringr::str_wrap(enrichment_results$GeneSet, width = 50)
enrichment_results$GeneSet <- factor(enrichment_results$GeneSet, levels = rev(enrichment_results$GeneSet))  # reverse for plotting top to bottom

ggplot(enrichment_results, aes(x = GeneSet, y = -log10(adjP), size = N_overlap, color = Category)) +
  geom_point() +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(x = "Gene Set", y = "-log10(Adjusted p-value)", 
       title = "Top 100 CpGs", 
       size = "Number of overlaps") +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave("top100_fuma.pdf", width = 8, height = 6)

## No significant results for the positively associated CpGs among top 100 CpGs.

## Load the negatively associated CpGs among top 100 CpGs
enrichment_results <- read.csv("top100_fuma_neg.csv")

enrichment_results$GeneSet <- enrichment_results$GeneSet %>%
  str_replace_all("_", " ")

enrichment_results$Category <- factor(enrichment_results$Category, 
                                      levels = c("Canonical_Pathways", "Cell_type_signature", "Hallmark_gene_sets", "Immunologic_signatures"))
enrichment_results <- enrichment_results %>%
  arrange(Category, adjP)
enrichment_results$GeneSet <- stringr::str_wrap(enrichment_results$GeneSet, width = 60)
enrichment_results$GeneSet <- factor(enrichment_results$GeneSet, levels = rev(enrichment_results$GeneSet))  # reverse for plotting top to bottom

ggplot(enrichment_results, aes(x = GeneSet, y = -log10(adjP), size = N_overlap, color = Category)) +
  geom_point() +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(x = "Gene Set", y = "-log10(Adjusted p-value)", 
       title = "Negatively associated CpGs among top 100 CpGs", 
       size = "Number of overlaps") +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave("top100_negative_fuma.pdf", width = 8, height = 7)
