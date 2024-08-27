# Load the necessary libraries
library(ComplexHeatmap) # for creating heatmaps
library(circlize) # for circular visualisation
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.

# Import the data
enrichr_table <- read.delim(snakemake@input[[1]], header = TRUE, sep = "\t", row.names = 1)
enrichr_term <- as.character(snakemake@input[[2]])
log2_mat <- as.matrix(read.table(snakemake@input[[3]], header = TRUE, row.names = 1, check.names = FALSE)) # diffExpr.P${fdr}_C${foldchange}.matrix
log2_comp <- read.delim(snakemake@input[[5]], header = FALSE, sep = "\t", row.names = 1) # sorted_logfolds.txt
diff_genes <- read.delim(snakemake@input[[4]], header = TRUE, sep = "\t")

# Select the genes for the term of interest and create a vector
enrichr_genes <- enrichr_table %>%
  filter(Term == enrichr_term) %>% # Select the term of interest
  pull(Genes) %>% # Extract the Genes column as a character vector
  str_split(";") %>% # Split the strings where there is a semicolon
  unlist() # Unlist the resulting list to create a vector

# Filter the diff_genes data frame to only include the genes that are in the enrichr_genes vector
diff_genes_2 <- diff_genes %>%
  filter(top_blastp_hit %in% enrichr_genes)

# Convert the log2_mat matrix to a data frame
log2_df <- as.data.frame(log2_mat)

# Filter log2_df by the trinity_id column of diff_genes_2
filtered_log2_df <- log2_df %>%
  filter(row.names(.) %in% diff_genes_2$trinity_id)

# Reset row names of filtered_log2_df to create a proper column
filtered_log2_df <- filtered_log2_df %>%
  rownames_to_column("trinity_id")

# Perform left join between filtered_log2_df and diff_genes_2
merged_df <- left_join(filtered_log2_df, diff_genes_2, by = "trinity_id")

# Remove the trinity_id column
merged_df <- merged_df %>%
  select(-trinity_id)

# Add two columns to the data frame with the average log2 fold change for rt_reps and 16_reps
merged_df$mean_rt_reps <- rowMeans(merged_df[,grep("rt_rep", colnames(merged_df))])
merged_df$mean_16_reps <- rowMeans(merged_df[,grep("16_rep", colnames(merged_df))])

# Keep only the top_blasp_hit, mean_rt_reps, and mean_16_reps columns
merged_df <- merged_df %>%
  select(top_blastp_hit, mean_rt_reps, mean_16_reps)

# Select the appropriate column input for the mean log2 heatmap
hclust_rep_rt <- merged_df$mean_rt_reps
hclust_rep_16 <- merged_df$mean_16_reps

# Create heatmaps for mean log2 values
heatmap_hclust_rep_rt <- Heatmap(hclust_rep_rt, name = "rt_reps", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16 <- Heatmap(hclust_rep_16, name = "16_reps", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = TRUE, cluster_rows = FALSE)

# Draw the combined heatmap in one line
png(paste("Enrichr_heatmap_", "Cell cycle", ".png", sep = ""), height = 2000, width = 2000, res = 300); draw(heatmap_hclust_rep_rt + heatmap_hclust_rep_16); dev.off()

# Save the heatmap
