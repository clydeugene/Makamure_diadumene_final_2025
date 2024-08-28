# Load the necessary libraries
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(ggrepel) # for nice annotations
library(enrichR) # for connecting to Enrichr

# Import DGE results
original_df <- read.delim(snakemake@input[[1]], header = TRUE, sep = "\t", row.names = 1)


# Import significant genes from blastx against SwissProt Metazoa
sg <- read_delim(snakemake@input[[2]], col_names = TRUE, delim =  "\t", show_col_types = FALSE) %>%
  distinct()


# Import trinity gene list
gl <- scan(snakemake@input[[3]], "", skip = 1) %>%
  length()


# Select the relevant columns
df <- original_df %>%
  select("log2FoldChange", "padj")

# Add a column to the data frame to specify if they are up- or down- regulated
# (log2FoldChange respectively positive or negative)
df$diffexpressed <- "no"

# if log2FoldChange > 0.5 and p-adj < 0.01, set as "up"
df$diffexpressed[df$log2FoldChange > 0.5 & df$padj < 0.01] <- "up"

# if log2FoldChange < -0.5 and p-adj < 0.01, set as "down"
df$diffexpressed[df$log2FoldChange < -0.5 & df$padj < 0.01] <- "down"

# Count the number of genes that are up- or down- regulated
upregulated <- sum(df$diffexpressed == "up")
downregulated <- sum(df$diffexpressed == "down")
notsignificant <- sum(df$diffexpressed == "no")

merged_df <- df %>%
  rownames_to_column(var = "trinity_id") %>%
  left_join(sg, by = "trinity_id") %>%
  mutate(
    top_blastp_hit = ifelse(top_blastp_hit == ".", NA, as.character(top_blastp_hit))
  )

# Top blastp  genes
top_blastp_genes <- merged_df %>%
  filter(!is.na(top_blastp_hit)) %>%
  top_n(5, wt = -padj) %>%
  arrange(padj) %>%
  select(trinity_id, top_blastp_hit, log2FoldChange, padj, diffexpressed)

# Get a count of unique genes
unique_genes_count <- merged_df %>%
  filter(!is.na(top_blastp_hit)) %>%
  # top_n(5, wt = -padj) %>%
  arrange(padj) %>%
  select(trinity_id, log2FoldChange, padj, diffexpressed, top_blastp_hit) %>%
  pull(top_blastp_hit) %>%
  unique() %>%
  length()

# Define the theme
theme_set(theme_classic(base_size = 12) +
  theme(
        axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
        plot.title = element_text(hjust = 0.5, margin = margin(0,0,20,0), face = "bold", size = rel(1.1), color = 'black')
        ))

# Create the volcano plot
trinity_blastp_volcano <- ggplot(data = merged_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = top_blastp_hit)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') + # plotting these first so they are behind the points
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "gray", "red"), # to set the colors of our variable
                    labels = c(paste0("Downregulated at 16°C (n = ", downregulated, ")"), paste0("Not significant (n = ", notsignificant, ")"), paste0("Upregulated at 16°C (n = ", upregulated, ")"))) + # to set the labels in case we want to overwrite the categories from the data frame (up, down, no) )
  coord_cartesian(ylim = c(0, 70), xlim = c(-12, 12)) + # since some genes can have -log10padj of inf, we set these limits
  labs(color = '', #legend_title, 
       x = expression("log"[2]*"Fold Change"), y = expression("-log"[10]*"(padj)")) + 
  scale_x_continuous(breaks = seq(-12, 12, 2)) + # to customize the breaks in the x axis
  ggtitle("Cold temperature vs Room temperature (RT) DEGs") +
  geom_label_repel(data = top_blastp_genes, aes(label = top_blastp_hit), color = "black", label.size = NA, nudge_y = 14, box.padding = unit(1, "lines"), segment.color = "grey50", segment.curvature = 0, max.overlaps = Inf, fill = NA, direction = "x") # To show 5 labels for top genes

# Save the plot
png(filename = snakemake@output[[1]], width = 2500, height = 2500, units = "px", pointsize = 12, bg = "white", res = 300, type = "cairo")
trinity_blastp_volcano
invisible(dev.off())

ggsave(snakemake@output[[4]], trinity_blastp_volcano, width = 2500, height = 2500, units = "px", dpi = 300)

# Set Enrichr site
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    # listEnrichrSites()
    setEnrichrSite("Enrichr") # General Enrichr
}

# Get Enrichr libraries
if (websiteLive) dbs <- listEnrichrDbs()$libraryName

# Filter to only libraries starting with "KEGG"
if (websiteLive) dbs <- dbs[grep("^KEGG", dbs)]


blastp_genes <- merged_df$top_blastp_hit %>%
  na.omit() %>%
  as.character() %>%
  unique()

# Run Enrichr analysis
if (websiteLive) {
  enrichr_results <- enrichr(blastp_genes, dbs)
}


# Get the KEGG_2021 results

# Sort the results
sorted_results <- enrichr_results[[6]][order(enrichr_results[[6]]$P.value), ]

# Save the results to a file
write.table(sorted_results, snakemake@output[[3]], sep = "\t", quote = FALSE, row.names = FALSE)



# Plot the KEGG_2021 results
kegg_2021_human_plot <- if (websiteLive) { plotEnrich(enrichr_results[[6]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "Enrichr KEGG_2021_Human") }

png(filename = snakemake@output[[2]], width = 2500, height = 2500, units = "px", pointsize = 12, bg = "white", res = 300, type = "cairo")
kegg_2021_human_plot
invisible(dev.off())

ggsave(snakemake@output[[5]], kegg_2021_human_plot, width = 2500, height = 2500, units = "px", dpi = 300)