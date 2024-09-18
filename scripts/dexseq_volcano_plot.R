# Load the necessary libraries
library(tidyverse)
library(pheatmap)
library(scales)
library(ggrepel)

# Read in the transcriptome from the commandline argument
transcriptome <- as.character(commandArgs(trailingOnly = TRUE))

# Set the working directory
setwd("data/processed/dexseq")

# Read in the dexseq results
results <- read.delim("diadumene_lineata_transcriptome.dexseq.results.dat", header = TRUE, sep = "\t", row.names = NULL)

# Read in the blastp results
blastp <- read.delim("../blastp/extract_orfs_longest_orfs.blastp.outfmt6", header = FALSE, sep = "\t", row.names = NULL)

# Filter the blastp results to only include the trinity_id and gene columns
blastp_filtered <- blastp |>
  select(trinity_id = V1, gene = V2) |>
  mutate(trinity_id = sub("_i.*", "", trinity_id)) |>
  separate(gene, into = c("db", "accession", "gene_symbol"), sep = "\\|") |>
  separate(gene_symbol, into = c("gene", "species"), sep = "\\_") |>
  select("trinity_id", "gene") |>
  distinct() |>
  group_by(trinity_id) |>
  summarize(gene = paste(unique(gene), collapse = ", "))

# Add a column to the results dataframe that combines the groupID and featureID
results$geneExon <- paste(results$groupID, results$featureID, sep = "_")

# Filter out rows with NA values in the padj column
results <- results |>
  filter(!is.na(padj))

# Select the padj and log2fold_rt_16 columns
filtered_results <- results |>
  select(geneExon, padj, log2fold_rt_16) |>
  column_to_rownames(var = "geneExon")

# Add a column to reverse the sign of the log2fold_rt_16 column so as to be consistently show 16 vs rt
filtered_results$log2fold_16_rt <- -filtered_results$log2fold_rt_16

# Remove the log2fold_rt_16 column
filtered_results <- filtered_results |>
  select(-log2fold_rt_16)

# Add a column to indicate if the row is differentially expressed
filtered_results$diffexpressed <- "no"

# if log2fold_rt_16 > 0.5 and p-adj < 0.01, set as "up"
filtered_results$diffexpressed[filtered_results$log2fold_16_rt > 0.5 & filtered_results$padj < 0.01] <- "up"

# if log2fold_rt_16 < -0.5 and p-adj < 0.01, set as "down"
filtered_results$diffexpressed[filtered_results$log2fold_16_rt < -0.5 & filtered_results$padj < 0.01] <- "down"

# Count the number of genes in each category
up_exons <- comma(sum(filtered_results$diffexpressed == "up"))
down_exons <- comma(sum(filtered_results$diffexpressed == "down"))
ns_exons <- comma(sum(filtered_results$diffexpressed == "no"))

# Add gene symbols to the filtered_results df
merged_df <- filtered_results |>
  rownames_to_column(var = "full_trinity_id") |>  # Temporarily store rownames as a column
  mutate(trinity_id = sub("_E.*", "", full_trinity_id)) |>  # Extract trinity_id by removing the exon part
  left_join(blastp_filtered, by = "trinity_id") |>  # Join based on trinity_id
  select(-trinity_id) |>  # Remove the temporary trinity_id column after the join
  column_to_rownames(var = "full_trinity_id")  # Convert back to the original row names

# Top blastp  genes
top_blastp_genes <- merged_df |>
  filter(!is.na(gene)) |>
  top_n(5, wt = -padj) |>
  arrange(padj) |>
  select(gene, log2fold_16_rt, padj, diffexpressed)

# Filter top blastp genes
top_blastp_genes <- merged_df |>
  filter(!is.na(gene), diffexpressed != "no", !str_detect(gene, ",")) |>  # Filter out NA genes and genes with multiple symbols and "no" in diffexpressed
  arrange(padj) |>  # Sort by padj
  top_n(-5, padj) |>  # Select top 5 rows with the smallest padj
  select(gene, log2fold_16_rt, padj, diffexpressed)  # Select the desired columns

# Make a volcano plot of those geneExons with abs(log2fold_rt_16 > 0.5) and padj < 0.01

# Define the theme
theme_set(theme_classic(base_size = 12) +
  theme(
        axis.title.y = element_text(face = "bold", margin = margin(0,10,0,0), size = 40, color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(10,0,0,0), size = 40, color = 'black'),
        plot.title = element_text(hjust = 0.5, margin = margin(0,0,10,0), face = "bold", size = 40, color = 'black'),
        axis.text = element_text(size = 36, color = 'black'),
        legend.text = element_text(size = 32),
        legend.title = element_text(size = 32, face = "bold")
        ))

# Create the volcano plot
volcano_plot <- ggplot(filtered_results, aes(x = log2fold_16_rt, y = -log10(padj), color = diffexpressed)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') + # plotting these first so they are behind the points
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "gray", "red"), # to set the colors of our variable
                    labels = c(paste0("Downregulated (n = ", down_exons, ")"), paste0("Not significant (n = ", ns_exons, ")"), paste0("Upregulated (n = ", up_exons, ")"))) + # to set the labels in case we want to overwrite the categories from the data frame (up, down, no) )
  coord_cartesian(ylim = c(0, 80), xlim = c(-12, 12)) + # since some genes can have -log10padj of inf, we set these limits
  labs(color = 'DEU at 16\u00B0C', #legend_title, 
       x = expression("log"[2]*"Fold Change"), y = expression("-log"[10]*"(padj)")) + 
  scale_x_continuous(breaks = seq(-12, 12, 6), ) + # to customize the breaks in the x axis
  # ggtitle("DEU")
  geom_label_repel(data = top_blastp_genes, aes(label = gene), color = "black", label.size = NA, nudge_y = 25, box.padding = unit(1, "lines"), segment.color = "black", segment.curvature = 0, max.overlaps = Inf, fill = NA, direction = "x", size = 10) + # To show 5 labels for top genes
  theme(legend.position = c(0.72, 0.88), legend.key.spacing.y = unit(0.2, "cm"))

# volcano_plot

# Save the plot as png and pdf files
ggsave(paste0("dexseq_", transcriptome, "_volcano_plot.png"), plot = volcano_plot, width = 3600, height = 2700, units = "px", dpi = 300)
ggsave(paste0("dexseq_", transcriptome, "_volcano_plot.pdf"), plot = volcano_plot, width = 3600, height = 2700, units = "px", dpi = 300)