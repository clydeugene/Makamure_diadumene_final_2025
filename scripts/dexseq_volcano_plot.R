# Load the necessary libraries
library(tidyverse)
library(pheatmap)

# Read in the transcriptome from the commandline argument
transcriptome <- as.character(commandArgs(trailingOnly = TRUE))

# Set the working directory
setwd("data/processed/dexseq")

# Read in the dexseq results
results <- read.delim("diadumene_lineata_transcriptome.dexseq.results.dat", header = TRUE, sep = "\t", row.names = NULL)

# Add a column to the results dataframe that combines the groupID and featureID
results$geneExon <- paste(results$groupID, results$featureID, sep = "_")

# Filter out rows with NA values in the padj column
results <- results %>%
  filter(!is.na(padj))

# Select the padj and log2fold_rt_16 columns
filtered_results <- results %>%
  select(geneExon, padj, log2fold_rt_16) %>%
  column_to_rownames(var = "geneExon")

# Add a column to reverse the sign of the log2fold_rt_16 column so as to be consistently show 16 vs rt
filtered_results$log2fold_16_rt <- -filtered_results$log2fold_rt_16

# Remove the log2fold_rt_16 column
filtered_results <- filtered_results %>%
  select(-log2fold_rt_16)

# Add a column to indicate if the row is differentially expressed
filtered_results$diffexpressed <- "no"

# if log2fold_rt_16 > 0.5 and p-adj < 0.01, set as "up"
filtered_results$diffexpressed[filtered_results$log2fold_16_rt > 0.5 & filtered_results$padj < 0.01] <- "up"

# if log2fold_rt_16 < -0.5 and p-adj < 0.01, set as "down"
filtered_results$diffexpressed[filtered_results$log2fold_16_rt < -0.5 & filtered_results$padj < 0.01] <- "down"

# Count the number of genes in each category
up_exons <- sum(filtered_results$diffexpressed == "up")
down_exons <- sum(filtered_results$diffexpressed == "down")
ns_exons <- sum(filtered_results$diffexpressed == "no")

# Define the theme
theme_set(theme_classic(base_size = 12) +
  theme(
        axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
        plot.title = element_text(hjust = 0.5, margin = margin(0,0,20,0), face = "bold", size = rel(1.1), color = 'black')
        ))

# Create the volcano plot
volcano_plot <- ggplot(filtered_results, aes(x = log2fold_16_rt, y = -log10(padj), color = diffexpressed)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') + # plotting these first so they are behind the points
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#377EB8", "gray", "#E41A1C"), # to set the colors of our variable
                    labels = c(paste0("Downregulated at 16\u00B0C (n = ", down_exons, ")"), paste0("Not significant (n = ", ns_exons, ")"), paste0("Upregulated at 16\u00B0C (n = ", up_exons, ")"))) + # to set the labels in case we want to overwrite the categories from the data frame (up, down, no) ) # nolint
  coord_cartesian(ylim = c(0, 60), xlim = c(-15, 15)) + # since some genes can have -log10padj of inf, we set these limits
  labs(color = '', #legend_title, 
       x = expression("log"[2]*"Fold Change"), y = expression("-log"[10]*"(padj)")) + 
  scale_x_continuous(breaks = seq(-15, 15, 3)) + # to customize the breaks in the x axis
  ggtitle("Cold (16\u00B0C) vs Room temperature differential exon usage (DEU)")

volcano_plot

# Save the plot as png and pdf files
ggsave(paste0("dexseq_", transcriptome, "_volcano_plot.png"), plot = volcano_plot, width = 2500, height = 2500, units = "px", dpi = 300)
ggsave(paste0("dexseq_", transcriptome, "_volcano_plot.pdf"), plot = volcano_plot, width = 2500, height = 2500, units = "px", dpi = 300)