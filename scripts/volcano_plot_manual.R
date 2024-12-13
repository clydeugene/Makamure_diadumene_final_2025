# This script was ran manually as an adaptation of the automated pipeline. This was more efficient as it allowed quick changes to be made during the final manuscript preparation.

# Set wd to the top level directory with the Snakefile
# setwd("")

# Load the necessary libraries
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(ggrepel) # for nice annotations
library(enrichR) # for connecting to Enrichr
library(scales) # for formatting the axis
library(ComplexHeatmap) # for creating heatmaps
library(circlize) # for creating heatmaps
library(gridtext) # for formatting text in heatmaps
library(patchwork) # for combining plots
library(reshape2) # for melting data frames

# Import DGE results
original_df <- read.delim("data/processed/DESeq2/trinity.gene.counts.matrix.16_vs_rt.DESeq2.DE_results", header = TRUE, sep = "\t", row.names = 1)


# Import significant genes from blastx against SwissProt Metazoa
sg <- read_delim("data/processed/DESeq2/DESeq2_gene/differentially_expressed_genes_with_symbols.txt", col_names = TRUE, delim =  "\t", show_col_types = FALSE) |>
  distinct()


# Import trinity gene list
gl <- scan("data/processed/DESeq2/DESeq2_gene/differentially_expressed_genes.txt", "", skip = 1) |>
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
upregulated <- sum(df$diffexpressed == "up")
downregulated <- sum(df$diffexpressed == "down")
notsignificant <- sum(df$diffexpressed == "no")

merged_df <- df |>
  rownames_to_column(var = "trinity_id") |>
  left_join(sg, by = "trinity_id") |>
  mutate(
    top_blastp_hit = ifelse(top_blastp_hit == ".", NA, as.character(top_blastp_hit))
  )

# Select top_blasp_hit entries that are not NA and where diffexpressed is not "no"
diffexpressed_genes <- merged_df |>
  filter(!is.na(top_blastp_hit) & diffexpressed != "no") |>
  select(top_blastp_hit, trinity_id, log2FoldChange, padj, diffexpressed)

# Save the results to a file
write.table(diffexpressed_genes, "data/processed/DESeq2/DESeq2_gene/differentially_expressed_genes_with_log2FC.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Create a vector diffexpressed_genes
diffexpressed_genes_vector <- diffexpressed_genes |>
  pull(top_blastp_hit) |>
  unique()

# Show the vector as a list that can be copied and pasted online
diffexpressed_genes_vector |>
  as.list() |>
  str_c(collapse = ", ") |>
  writeLines(con = stdout())

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

# Create the volcano plot
trinity_blastp_volcano <- ggplot(data = merged_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = top_blastp_hit)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') + # plotting these first so they are behind the points
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "gray", "red"), # to set the colors of our variable
                     labels = c(paste0("Downregulated (n = ", downregulated, ")"), paste0("Not significant (n = ", notsignificant, ")"), paste0("Upregulated (n = ", upregulated, ")"))) + # to set the labels in case we want to overwrite the categories from the data frame (up, down, no) )
  coord_cartesian(ylim = c(0, 70), xlim = c(-12, 12)) + # since some genes can have -log10padj of inf, we set these limits
  labs(color = 'DEGs at 16\u00B0C', #legend_title, 
       x = expression("log"[2]*"Fold Change"), y = expression("-log"[10]*"(padj)")) + 
  scale_x_continuous(breaks = seq(-12, 12, 6)) + # to customize the breaks in the x axis
  # ggtitle("DEGs") +
  geom_label_repel(data = top_blastp_genes, aes(label = top_blastp_hit), color = "black", label.size = NA, nudge_y = 25, box.padding = unit(1, "lines"), segment.color = "black", segment.curvature = 0, max.overlaps = Inf, fill = NA, direction = "x", size = 10) + # To show 5 labels for top genes
  theme(legend.position = c(0.72, 0.88), legend.key.spacing.y = unit(0.2, "cm"))


# Set Enrichr site
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  # listEnrichrSites()
  setEnrichrSite("Enrichr") # General Enrichr
}

# Get Enrichr libraries
if (websiteLive) dbs <- listEnrichrDbs()$libraryName

# Filter to only libraries starting with "GO"
if (websiteLive) dbs_go <- dbs[grep("^GO", dbs)]

# Filter to only libraries starting with "KEGG"
if (websiteLive) dbs <- dbs[grep("^KEGG", dbs)]



blastp_genes <- merged_df$top_blastp_hit |>
  na.omit() |>
  as.character() |>
  unique()


# Run Enrichr analysis for GO and KEGG
# Run Enrichr analysis
if (websiteLive) {
  enrichr_results <- enrichr(blastp_genes, dbs_go)
}

if (websiteLive) {
  enrichr_results_kegg <- enrichr(blastp_genes, dbs)
}

# Filter out the GO_Biological_Process_2023 results from the list of data frames
sorted_results <- enrichr_results$GO_Biological_Process_2023 |>
  arrange(P.value)

go_2023 <- enrichr_results$GO_Biological_Process_2023 |>
  arrange(P.value)

kegg_2021 <- enrichr_results_kegg$KEGG_2021_Human |>
  arrange(P.value)

# Select only Mesenchymal Cell Differentiation (GO:0048762) and Positive regulation of meiotic nuclear division from the GO_Biological_Process_2023 results
go_2023 <- go_2023 |>
  filter(Term %in% c("Mesenchymal Cell Differentiation (GO:0048762)", "Positive Regulation Of Meiotic Nuclear Division (GO:0045836)"))

# Select only Cell cycle and Protein digestion and absorption from the KEGG_2021 results
kegg_2021 <- kegg_2021 |>
  filter(Term %in% c("Cell cycle", "Protein digestion and absorption")) |>
  mutate(Term = paste0(Term, " (KEGG 2021)"))

# Merge the filtered results
kegg_and_go <- rbind(go_2023, kegg_2021) |>
  arrange(P.value)

# Plot the merged results
kegg_and_go_plot <- plotEnrich(kegg_and_go, showTerms = 20, numChar = 70, y = "Count", orderBy = "P.value", title = "Enrichr GO Biological Process 2023 and KEGG 2021 Human") +
  theme(
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # Remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # Remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # Remove vertical minor grid lines
    axis.text = element_text(size = 15, color = 'black'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = 20, color = 'black'),
    plot.title = element_text(hjust = 0.85, margin = margin(0,0,20,0), face = "bold", size = 24, color = 'black'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 15, color = "black", face = "bold")
  ) +
  scale_fill_gradientn(
    colours = c("purple", "white"),
    limits = c(0, 0.06),
    breaks = c(0, 0.02, 0.04, 0.06),
  ) +
  ggtitle("Functional Enrichment Analysis")

# Save the kegg_and_go_plot to a file
ggsave("data/processed/DESeq2/DESeq2_gene/gene_enrichr_kegg_and_go_plot.pdf", kegg_and_go_plot, width = 3600, height = 1200, units = "px", dpi = 300)

# Create heatmaps for the differentially expressed genes in each term of kegg_and_go

# Extract the genes for each term from the genes column and add them to a list of gene vectors. Name each vector with the first two words of the term separated by an underscore
kegg_and_go_genes <- kegg_and_go |>
  mutate(genes = strsplit(Genes, ";")) |>
  pull(genes) |>
  set_names(kegg_and_go$Term) |>
  map(~ unlist(.x))

# For each item in the list, create another list extract the genes from the merged_df data frame and create a new list with data frames of matching genes and their log2FoldChange and padj values
kegg_and_go_genes_dfs <- kegg_and_go_genes |>
  map(~ merged_df |>
        filter(top_blastp_hit %in% .x) |>
        select(top_blastp_hit, log2FoldChange, padj))


# Using the complexheatmap package, create a heatmap for each term in kegg_and_go_genes_dfs
# Function to create heatmaps for each data frame in kegg_and_go_genes_dfs
create_heatmaps <- function(gene_dfs) {
  for (term in names(gene_dfs)) {
    # Extract the data frame for the current term
    df <- gene_dfs[[term]]
    
    # Check if the data frame is not empty
    if (nrow(df) > 0) {
      # Create a matrix with log2FoldChange and padj values
      gene_matrix <- as.matrix(df[, c("log2FoldChange", "padj")])
      
      # Get the padj and log2FoldChange values
      padj_values <- gene_matrix[, "padj"]
      log2fold_values <- gene_matrix[, "log2FoldChange"]
      
      # Create a heatmap of log2FoldChange values
      log2fold_heatmap <- Heatmap(log2fold_values, 
                                  name = "log\u2082 Fold Change", 
                                  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
                                  row_labels = df$top_blastp_hit,
                                  column_labels = gt_render("<span style='font-size:30pt'>Diff. Expr.</span>"),
                                  cluster_rows = FALSE, 
                                  show_row_dend = FALSE, 
                                  cluster_columns = FALSE, 
                                  width = unit(5, "cm"))
      
      # Create a heatmap of padj values
      padj_heatmap <- Heatmap(padj_values, 
                              name = "Adjusted p-value",
                              col = colorRamp2(c(0, 0.05), c("blue", "white")),
                              column_labels = c(""),
                              cluster_rows = FALSE, 
                              show_row_dend = FALSE, 
                              cluster_columns = FALSE, 
                              width = unit(2, "cm"))
      
      # Draw the combined heatmap
      draw(
        padj_heatmap + log2fold_heatmap,
        heatmap_legend_side = "left",
        column_title = term
      )
    }
  }
}

# Call the function to create heatmaps
create_heatmaps(kegg_and_go_genes_dfs)

# Save the heatmaps to separate files



# Read in the hydra tables
tbl_s4 <- read.table("data/processed/DESeq2/DESeq2_gene/table_s4_sasp.txt", header = TRUE, sep = "\t")
tbl_s5 <- read.table("data/processed/DESeq2/DESeq2_gene/table_s5_ros.txt", header = TRUE, sep = "\t")
tbl_s6 <- read.table("data/processed/DESeq2/DESeq2_gene/table_s6_genome_maintenance.txt", header = TRUE, sep = "\t")


# Create a list of the tables
tbls <- list(tbl_s4, tbl_s5, tbl_s6)

# Retain only the human_gene_name and log2FC_sexual_asexual columns
tbls <- map(tbls, ~ .x %>% select(human_gene_name, log2FC_sexual_asexual)) %>% # Remove asterisks from gene names
  map(~ .x %>% mutate(human_gene_name = gsub("\\*", "", human_gene_name))) %>% # Rename the log2FC_sexual_asexual column to log2FoldChange
  map(~ .x %>% rename(log2FoldChange = log2FC_sexual_asexual))

# Combine the tables into a single data frame
tbls_combined <- bind_rows(tbls, .id = NULL) %>% # Remove duplicates
  distinct() %>% # Remove rows with missing values
  filter(!is.na(log2FoldChange)) %>% # Rename log2FoldChange to log2FoldChange_hydra
  rename(log2FoldChange_hydra = log2FoldChange) %>% # Change log2FoldChange_hydra to dbl
  mutate(log2FoldChange_hydra = as.numeric(log2FoldChange_hydra))

# Read in the file with approved gene symbols from hugo.
# The data wrangling was done on my local machine and will need to be combined with this script later
diff_approv_symbols <- read.table("data/processed/DESeq2/DESeq2_gene/diff_gene_approv_symbols_with_log2FC.txt", header = TRUE, sep = " ")

# # Filter tbls_combined to only include genes in diff_approv_symbols
# tbls_combined <- tbls_combined %>%
#   filter(human_gene_name %in% diff_approv_symbols$Approved.symbol) %>% # Create a new column with log2FoldChange from diff_approv_symbols
#   left_join(diff_approv_symbols, by = c("human_gene_name" = "Approved.symbol")) %>% # Rename log2FoldChange to log2FoldChange_diadumene
#   rename(log2FoldChange_diadumene = log2FoldChange)

# Filter tbls_combined to only include genes in merged_df
tbls_combined <- tbls_combined %>%
  filter(human_gene_name %in% merged_df$top_blastp_hit) %>% # Create a new column with log2FoldChange from diff_approv_symbols
  left_join(merged_df, by = c("human_gene_name" = "top_blastp_hit")) %>% # Rename log2FoldChange to log2FoldChange_diadumene
  rename(log2FoldChange_diadumene = log2FoldChange)

# Create a scatter plot of absolute log2FoldChange_diadumene vs absolute log2FoldChange_hydra
scatter_plot <- ggplot(tbls_combined, aes(x = abs(log2FoldChange_diadumene), y = abs(log2FoldChange_hydra))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "| log\u2082 Fold Change | Diadumene", y = "| log\u2082 Fold Change | Hydra") + #Set the axis limits
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  theme_classic()


scatter_plot




# Change the approach for getting combining the tables and getting the scatter plot
# First, get the candidate_human_gene_names column and log2FoldChange column from each table
tbls_combined_candidate <- tbls |> 
  map(~ .x %>% select(candidate_human_gene_names, log2FC_sexual_asexual)) 
  
# Combine the tables into a single data frame
tbls_combined_candidate_df <- bind_rows(tbls_combined_candidate, .id = NULL) |>
  distinct() |>
  filter(!is.na(log2FC_sexual_asexual)) |>
  rename(log2FoldChange_hydra = log2FC_sexual_asexual) |>
  mutate(log2FoldChange_hydra = as.numeric(log2FoldChange_hydra))

# Separate comma-separated values in the candidate_human_gene_names column into separate rows
tbls_combined_candidate_df <- tbls_combined_candidate_df |>
  separate_rows(candidate_human_gene_names, sep = ",") |>
  filter(!is.na(candidate_human_gene_names)) |>
  distinct()

# Filter tbls_combined to only include genes in merged_df
tbls_combined_candidate_df <- tbls_combined_candidate_df %>%
  filter(candidate_human_gene_names %in% merged_df$top_blastp_hit) %>% # Create a new column with log2FoldChange from diff_approv_symbols
  left_join(merged_df, by = c("candidate_human_gene_names" = "top_blastp_hit")) %>% # Rename log2FoldChange to log2FoldChange_diadumene
  rename(log2FoldChange_diadumene = log2FoldChange)

# Create a scatter plot of absolute log2FoldChange_diadumene vs absolute log2FoldChange_hydra
scatter_plot <- ggplot(tbls_combined_candidate_df, aes(x = log2FoldChange_diadumene, y = log2FoldChange_hydra)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + # Add horizontal grid line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + # Add vertical grid line at x = 0
  labs(x = "log\u2082 Fold Change Diadumene", y = "log\u2082 Fold Change Hydra") + #Set the axis limits
  # coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  theme_classic()


scatter_plot














# Save the results to a file
write.table(sorted_results, "data/processed/DESeq2/DESeq2_gene/enrichr_results_GO_Biological_Process_2017b.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot the GO_Biological_Process_2023 results
go_biological_process_2023_plot <- plotEnrich(enrichr_results$GO_Biological_Process_2023, showTerms = 63, numChar = 40, y = "Count", orderBy = "P.value", title = "Enrichr GO_Biological_Process_2023")

go_plot <- go_biological_process_2023_plot + 
  theme(
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # Remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # Remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # Remove vertical minor grid lines
    axis.text = element_text(size = 5, color = 'black'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = 20, color = 'black'),
    plot.title = element_text(hjust = 0.5, margin = margin(0,0,20,0), face = "bold", size = 24, color = 'black'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 15, color = "black", face = "bold")
  ) +
  scale_fill_gradientn(
    labels = label_number(accuracy = 0.01),
    colours = c("red", "blue")
  ) +
  ggtitle("Enrichr GO Biological Process 2023")

go_plot

# Run Enrichr analysis
if (websiteLive) {
  enrichr_results <- enrichr(blastp_genes, dbs)
}



# Get the KEGG_2021 results

# Sort the results
sorted_results <- enrichr_results[[6]][order(enrichr_results[[6]]$P.value), ]

# Save the results to a file
# write.table(sorted_results, snakemake@output[[3]], sep = "\t", quote = FALSE, row.names = FALSE)



# Plot the KEGG_2021 results
kegg_2021_human_plot <- if (websiteLive) { plotEnrich(enrichr_results[[6]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "Enrichr KEGG_2021_Human") }

kegg_plot <- kegg_2021_human_plot + 
  theme(
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # Remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # Remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # Remove vertical minor grid lines
    axis.text = element_text(size = 15, color = 'black'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = 20, color = 'black'),
    plot.title = element_text(hjust = 0.5, margin = margin(0,0,20,0), face = "bold", size = 24, color = 'black'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 15, color = "black", face = "bold")
  ) +
  scale_fill_gradientn(
    labels = label_number(accuracy = 0.01),
    colours = c("red", "blue")
  ) +
  ggtitle("Enrichr KEGG 2021 Human")

kegg_plot

