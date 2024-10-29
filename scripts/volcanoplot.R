# Load the necessary libraries
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(ggrepel) # for nice annotations
library(enrichR) # for connecting to Enrichr
library(scales) # for formatting the axis

# Import DGE results
original_df <- read.delim(snakemake@input[[1]], header = TRUE, sep = "\t", row.names = 1)


# Import significant genes from blastx against SwissProt Metazoa
sg <- read_delim(snakemake@input[[2]], col_names = TRUE, delim =  "\t", show_col_types = FALSE) |>
  distinct()

# Import cell aging data

cell_age <- read_delim("cellage3.tsv", delim = "\t", col_names = TRUE, show_col_types = FALSE) |>
  select("Gene symbol", "Type of senescence", "Senescence Effect")

cell_sig <- read_delim("signatures1.csv", delim = ";", col_names = TRUE) |>
  select("gene_symbol", "ovevrexp", "underexp")


# Import trinity gene list
gl <- scan(snakemake@input[[3]], "", skip = 1) |>
  length()


# Select the relevant columns
df <- original_df |>
  select("log2FoldChange", "padj")

# Add a column to the data frame to specify if they are up- or down- regulated
# (log2FoldChange respectively positive or negative)
df$diffexpressed <- "no"

# if log2FoldChange > 1 and p-adj < 0.01, set as "up"
df$diffexpressed[df$log2FoldChange > 1 & df$padj < 0.01] <- "up"

# if log2FoldChange < -1 and p-adj < 0.01, set as "down"
df$diffexpressed[df$log2FoldChange < -1 & df$padj < 0.01] <- "down"

# Count the number of genes that are up- or down- regulated
upregulated <- comma(sum(df$diffexpressed == "up"))
downregulated <- comma(sum(df$diffexpressed == "down"))
notsignificant <- comma(sum(df$diffexpressed == "no"))

merged_df <- df |>
  rownames_to_column(var = "trinity_id") |>
  left_join(sg, by = "trinity_id") |>
  mutate(
    top_blastp_hit = ifelse(top_blastp_hit == ".", NA, as.character(top_blastp_hit))
  )

# Select the top_blastp_hits that are not NA and where diffexpressed is not "no" and join with cell_age
diff_expr_genes <- merged_df |>
  filter(!is.na(top_blastp_hit) & diffexpressed != "no") |>
  select(top_blastp_hit, log2FoldChange, padj, diffexpressed) |>
  unique() |>  
  left_join(cell_age, by = c("top_blastp_hit" = "Gene symbol")) |>
  left_join(cell_sig, by = c("top_blastp_hit" = "gene_symbol")) |>
  filter(rowSums(across(everything(), is.na)) <= 3) |>
  unique() |>
  arrange(desc(abs(log2FoldChange)), diffexpressed, ovevrexp, top_blastp_hit) |>
  distinct(top_blastp_hit, .keep_all = TRUE) |>
  select(top_blastp_hit,
    padj,
    log2FoldChange,
    diffexpressed,
    `Senescence Effect`,
    ovevrexp,
    underexp
  ) |>
  # Create a flag column to identify rows to remove
  mutate(remove_cellage = ifelse(
    (diffexpressed == "up" & `Senescence Effect` == "Inhibits") | 
    (diffexpressed == "down" & `Senescence Effect` == "Induces"),
    TRUE, FALSE)) |>
  mutate(remove_expr = ifelse(
    (diffexpressed == "down" & ovevrexp == 1) | 
    (diffexpressed == "up" & underexp == 1),
    TRUE, FALSE)) |>
  filter(!(is.na(remove_cellage) & is.na(remove_expr)) & !(remove_cellage == TRUE & remove_expr == TRUE)) |>
  select(-remove_cellage, -remove_expr)
    
    

# Save the results to a file
write.table(diff_expr_genes,
  snakemake@output[[6]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Select the top 5 differentially expressed genes
top_5_diff_genes <- diff_expr_genes |>
  top_n(5, wt = -padj)


# Top blastp  genes
top_blastp_genes <- merged_df |>
  filter(!is.na(top_blastp_hit)) |>
  top_n(5, wt = -padj) |>
  arrange(padj) |>
  select(trinity_id, top_blastp_hit, log2FoldChange, padj, diffexpressed)

# Get a count of unique genes
unique_genes_count <- merged_df |>
  filter(!is.na(top_blastp_hit)) |>
  # top_n(5, wt = -padj) |>
  arrange(padj) |>
  select(trinity_id, log2FoldChange, padj, diffexpressed, top_blastp_hit) |>
  pull(top_blastp_hit) |>
  unique() |>
  length()

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

# Set label limits
y_limits <- c(40, NA)

# Create the volcano plot
trinity_blastp_volcano <- ggplot(data = merged_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = top_blastp_hit)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') + # plotting these first so they are behind the points
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "gray", "red"), # to set the colors of our variable
                    labels = c(paste0("Down at 16\u00B0C (n = ", downregulated, ")"), paste0("Not significant (n = ", notsignificant, ")"), paste0("Up at 16\u00B0C (n = ", upregulated, ")"))) + # to set the labels in case we want to overwrite the categories from the data frame (up, down, no) )
  coord_cartesian(ylim = c(0, 70), xlim = c(-12, 12)) + # since some genes can have -log10padj of inf, we set these limits
  labs(color = 'DEGs', #legend_title, 
       x = expression("log"[2]*" Fold Change"), y = expression("-log"[10]*"(Adjusted P-value)")) + 
  scale_x_continuous(breaks = seq(-12, 12, 6)) + # to customize the breaks in the x axis
  # ggtitle("DEGs") +
  geom_label_repel(data = top_5_diff_genes, aes(label = top_blastp_hit), color = "black", label.size = NA, nudge_y = 25, box.padding = unit(1, "lines"), segment.color = "black", segment.curvature = 0, max.overlaps = Inf, fill = NA, direction = "x", size = 10, ylim = y_limits, hjust = 0) + # To show 5 labels for top genes
  theme(legend.position = c(0.72, 0.88), legend.key.spacing.y = unit(0.2, "cm"))

# # Save the plot
# png(filename = snakemake@output[[1]], width = 2500, height = 2500, units = "px", pointsize = 12, bg = "white", res = 300, type = "cairo")
# trinity_blastp_volcano
# invisible(dev.off())

ggsave(snakemake@output[[1]], trinity_blastp_volcano, width = 3600, height = 2700, units = "px", dpi = 300)
ggsave(snakemake@output[[4]], trinity_blastp_volcano, width = 3600, height = 2700, units = "px", dpi = 300)

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


blastp_genes <- merged_df$top_blastp_hit |>
  na.omit() |>
  as.character() |>
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
kegg_2021_human_plot <- if (websiteLive) { plotEnrich(enrichr_results[[6]], showTerms = 10, numChar = 52, y = "Count", orderBy = "P.value", title = "Enrichr KEGG_2021_Human") }

kegg_plot <- kegg_2021_human_plot + 
  theme(
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # Remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # Remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # Remove vertical minor grid lines
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5, margin = margin(20,0,0,0), size = 20, color = "black"),
    plot.title = element_text(hjust = 0.5, margin = margin(0,0,20,0), face = "bold", size = 24, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 15, color = "black", face = "bold")
  ) +
  scale_fill_gradientn(
    labels = label_number(accuracy = 0.01),
    colours = c("red", "blue")
  ) +
  ggtitle("Enrichr KEGG 2021 Human")

# png(filename = snakemake@output[[2]], width = 2500, height = 2500, units = "px", pointsize = 12, bg = "white", res = 300, type = "cairo")
# kegg_2021_human_plot
# invisible(dev.off())

ggsave(snakemake@output[[2]], kegg_plot, width = 3600, height = 1200, units = "px", dpi = 300)
ggsave(snakemake@output[[5]], kegg_plot, width = 3600, height = 1200, units = "px", dpi = 300)