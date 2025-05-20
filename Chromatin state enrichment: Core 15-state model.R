library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(minfiDataEPIC)

rg.epic <- get(data("RGsetEPIC"))
annotation <- as.data.frame(getAnnotation(rg.epic))

# =========================
# 1. Define input CpGs
# =========================

# Load significant CpGs
top100_sig_cpgs <- read.csv("top100.csv")
top100_sig_cpgs_pos <- subset(top100_sig_cpgs, beta > 0)
top100_sig_cpgs_neg <- subset(top100_sig_cpgs, beta < 0)

top100_sig_cpgs <- top100_sig_cpgs$cpg
top100_sig_cpgs_pos <- top100_sig_cpgs_pos$cpg
top100_sig_cpgs_neg <- top100_sig_cpgs_neg$cpg

# Load background CpGs
bg_cpgs <- read.csv("all BMI p values to plot manhattan.csv")$SNP

# Prepare significant CpGs with chr/pos
top100 <- annotation[rownames(annotation) %in% top100_sig_cpgs, c("chr", "pos")]
top100_pos <- annotation[rownames(annotation) %in% top100_sig_cpgs_pos, c("chr", "pos")]
top100_neg <- annotation[rownames(annotation) %in% top100_sig_cpgs_neg, c("chr", "pos")]

# Prepare background_df (all tested CpGs with chr/pos)
background_df <- annotation[rownames(annotation) %in% bg_cpgs, c("chr", "pos")]

# =========================
# 2. Enrichment Function
# =========================

run_chromatin_enrichment <- function(cpg_df, chrom_gr, background_df, tissue_id, output_dir = NULL) {
  cpg_gr <- GRanges(seqnames = cpg_df$chr,
                    ranges = IRanges(start = cpg_df$pos, end = cpg_df$pos + 1))
  
  hits <- findOverlaps(cpg_gr, chrom_gr)
  overlap_states <- mcols(chrom_gr[subjectHits(hits)])$name
  obs_counts <- table(overlap_states)
  obs_df <- as.data.frame(obs_counts)
  colnames(obs_df) <- c("Chromatin_State", "Observed")
  
  bg_gr <- GRanges(seqnames = background_df$chr,
                   ranges = IRanges(start = background_df$pos, end = background_df$pos + 1))
  bg_hits <- findOverlaps(bg_gr, chrom_gr)
  bg_states <- mcols(chrom_gr[subjectHits(bg_hits)])$name
  bg_counts <- table(bg_states)
  bg_df <- as.data.frame(bg_counts)
  colnames(bg_df) <- c("Chromatin_State", "Background")
  
  enrich_df <- merge(obs_df, bg_df, by = "Chromatin_State", all = TRUE)
  enrich_df[is.na(enrich_df)] <- 0
  enrich_df$Observed <- as.numeric(enrich_df$Observed)
  enrich_df$Background <- as.numeric(enrich_df$Background)
  
  enrich_df$Ratio <- (enrich_df$Observed / sum(enrich_df$Observed)) /
    (enrich_df$Background / sum(enrich_df$Background))
  
  enrich_df$pval <- apply(enrich_df, 1, function(row) {
    observed <- as.numeric(row["Observed"])
    background <- as.numeric(row["Background"])
    total_obs <- sum(enrich_df$Observed)
    total_bg <- sum(enrich_df$Background)
    
    fisher.test(matrix(c(
      observed,
      background,
      total_obs - observed,
      total_bg - background
    ), nrow = 2))$p.value
  })
  
  enrich_df$FDR <- p.adjust(enrich_df$pval, method = "fdr")
  sig_enrich_df <- enrich_df[enrich_df$FDR < 0.05, ]
  
  p <- ggplot(sig_enrich_df, aes(x = reorder(Chromatin_State, -Ratio), y = Ratio)) +
    geom_bar(stat = "identity", fill = "tomato") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    coord_flip() +
    labs(title = paste0(tissue_id, ": Chromatin State Enrichment (FDR < 0.05)"),
         x = "Chromatin State", y = "Observed / Expected Ratio") +
    theme_bw(base_size = 14)
  
  # Save results
  if (!is.null(output_dir)) {
    ggsave(filename = file.path(output_dir, paste0(tissue_id, "_enrichment_plot.png")), plot = p, width = 8, height = 6)
    write.csv(enrich_df, file = file.path(output_dir, paste0(tissue_id, "_enrichment_table.csv")), row.names = FALSE)
  }
  
  return(list(tissue = tissue_id, enrichment_table = enrich_df, significant = sig_enrich_df, plot = p))
}

# ============================
# 3. Loop over selected EIDs
# ============================
# Brain tissues
selected_eids <- c("E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074")

# Blood tissues
selected_eids <- c("E062", "E034", "E045", "E033", "E044", "E043", "E039", "E041",
                   "E042", "E040", "E037", "E048", "E038", "E047")

bed_dir <- "~/Downloads/all.mnemonics.bedFiles/"
bed_files <- list.files(bed_dir, pattern = "*.bed.gz", full.names = TRUE)
selected_files <- bed_files[substr(basename(bed_files), 1, 4) %in% selected_eids]

output_dir <- "chromatin_enrichment_results_selected"
dir.create(output_dir, showWarnings = FALSE)

# Positively associated CpGs among top 100 CpGs
top100_results_list2 <- list()

for (bed_file in selected_files) {
  tissue_id <- gsub("_15_coreMarks_mnemonics.bed.gz", "", basename(bed_file))
  message("Processing: ", tissue_id)
  
  chrom_gr <- import(bed_file)
  
  res <- run_chromatin_enrichment(
    cpg_df = top100,
    chrom_gr = chrom_gr,
    background_df = background_df,
    tissue_id = tissue_id,
    output_dir = output_dir
  )
  
  top100_results_list2[[tissue_id]] <- res
}

# Negatively associated CpGs among top 100 CpGs
top100_neg_results_list2 <- list()

for (bed_file in selected_files) {
  tissue_id <- gsub("_15_coreMarks_mnemonics.bed.gz", "", basename(bed_file))
  message("Processing: ", tissue_id)
  
  chrom_gr <- import(bed_file)
  
  res <- run_chromatin_enrichment(
    cpg_df = top100_neg,
    chrom_gr = chrom_gr,
    background_df = background_df,
    tissue_id = tissue_id,
    output_dir = output_dir
  )
  
  top100_neg_results_list2[[tissue_id]] <- res
}

top100_pos_results_list2 <- list()

for (bed_file in selected_files) {
  tissue_id <- gsub("_15_coreMarks_mnemonics.bed.gz", "", basename(bed_file))
  message("Processing: ", tissue_id)
  
  chrom_gr <- import(bed_file)
  
  res <- run_chromatin_enrichment(
    cpg_df = top100_pos,
    chrom_gr = chrom_gr,
    background_df = background_df,
    tissue_id = tissue_id,
    output_dir = output_dir
  )
  
  top100_pos_results_list2[[tissue_id]] <- res
}


# =========================
# 4. Summary Heatmap
# =========================
# Top 100 - Brain tissues
top100_summary_df <- do.call(rbind, lapply(names(top100_results_list), function(tissue_id) {
  df <- top100_results_list[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))

top100_pos_summary_df <- do.call(rbind, lapply(names(top100_pos_results_list), function(tissue_id) {
  df <- top100_pos_results_list[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))

top100_neg_summary_df <- do.call(rbind, lapply(names(top100_neg_results_list), function(tissue_id) {
  df <- top100_neg_results_list[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))

# Top 100 - Blood tissues
top100_summary_df2 <- do.call(rbind, lapply(names(top100_results_list2), function(tissue_id) {
  df <- top100_results_list2[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))

top100_pos_summary_df2 <- do.call(rbind, lapply(names(top100_pos_results_list2), function(tissue_id) {
  df <- top100_pos_results_list2[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))

top100_neg_summary_df2 <- do.call(rbind, lapply(names(top100_neg_results_list2), function(tissue_id) {
  df <- top100_neg_results_list2[[tissue_id]]$enrichment_table
  df$tissue <- tissue_id
  df
}))


# Create a mapping vector
eid_mapping <- c(
  "E071" = "Hippocampus middle",
  "E074" = "Substantia nigra",
  "E068" = "Anterior caudate",
  "E069" = "Cingulate gyrus",
  "E072" = "Inferior temporal lobe",
  "E067" = "Angular gyrus",
  "E073" = "Dorsolateral prefrontal cortex",
  "E070" = "Germinal matrix"
)

# Apply mapping to both dataframes
top100_pos_summary_df$tissue <- eid_mapping[top100_pos_summary_df$tissue]
top100_neg_summary_df$tissue <- eid_mapping[top100_neg_summary_df$tissue]

eid_to_name <- c(
  "E062" = "PBMCs",
  "E034" = "CD3+ T cells (PB)",
  "E045" = "CD4+ Tmem (PB)",
  "E033" = "CD3+ T cells (cord)",
  "E044" = "Tregs",
  "E043" = "CD4+ Th cells",
  "E039" = "Naive CD4+ T (CD45RA+)",
  "E041" = "Stim. CD4+ Th (IL17-)",
  "E042" = "Stim. Th17 (IL17+)",
  "E040" = "CD4+ Tmem (CD45RO+)",
  "E037" = "CD4+ Tmem",
  "E048" = "CD8+ Tmem",
  "E038" = "Naive CD4+ T (unsorted)",
  "E047" = "Naive CD8+ T"
)

top100_pos_summary_df2$tissue <- eid_to_name[top100_pos_summary_df2$tissue]
top100_neg_summary_df2$tissue <- eid_to_name[top100_neg_summary_df2$tissue]

## ------------------------------------------------------------
make_enrichment_heatmap <- function(summary_df, title = "") {
  
  ## 1a. Odds-ratio wide matrix
  ratio_mat <- summary_df %>%
    select(tissue, Chromatin_State, Ratio) %>%
    pivot_wider(names_from = Chromatin_State, values_from = Ratio, values_fill = 1) %>%
    tibble::column_to_rownames("tissue") %>%
    as.matrix()
  
  ## 1b. p-value wide matrix
  p_mat <- summary_df %>%
    select(tissue, Chromatin_State, pval) %>%
    pivot_wider(names_from = Chromatin_State, values_from = pval, values_fill = NA) %>%
    tibble::column_to_rownames("tissue") %>%
    as.matrix()
  
  ## 1c. Colour function
  col_fun <- colorRamp2(
    c(min(ratio_mat, na.rm = TRUE), 1, max(ratio_mat, na.rm = TRUE)),
    c("#ADD8E6", "white", "#FFA07A")
  )
  
  ## 1d. Build the heatmap
  ht <- Heatmap(
    ratio_mat,
    name = "Ratio",
    col = col_fun,
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill, col = "black"))  # box with border
      pval <- p_mat[i, j]
      if (!is.na(pval)) {
        stars <- if (pval < 0.001) {
          "***"
        } else if (pval < 0.01) {
          "**"
        } else if (pval < 0.05) {
          "*"
        } else {
          ""
        }
        if (stars != "") {
          grid::grid.text(stars, x = x, y = y, gp = grid::gpar(fontsize = 8, col = "black",fontface = "bold"))
        }
      }
    }
  )
  return(ht)
}


## ------------------------------------------------------------
##  2.  Build the two heatmaps
## ------------------------------------------------------------
ht_top100_pos<- make_enrichment_heatmap(top100_pos_summary_df, "Positively associated CpGs")
ht_top100_neg <- make_enrichment_heatmap(top100_neg_summary_df, "Negatively associated CpGs")


ht_top100_pos2<- make_enrichment_heatmap(top100_pos_summary_df2, "Positively associated CpGs")
ht_top100_neg2 <- make_enrichment_heatmap(top100_neg_summary_df2, "Negatively associated CpGs")

## ------------------------------------------------------------
##  3.  Draw or save
## ------------------------------------------------------------
# side-by-side
# draw(ht_top100 + ht_top200)

# stacked vertically
ht_list1 <- ht_top100_pos + ht_top100_neg
ht_list2 <- ht_top100_pos2 + ht_top100_neg2
draw(ht_list2)

pdf("top100_blood_enrichment_heatmap2.pdf", width = 8, height = 4)  
draw(ht_list2)                                                     
dev.off()      

