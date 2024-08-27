#!/bin/bash

# Input file paths
original_df_file="data/processed/DESeq2/trinity.gene.counts.matrix.16_vs_rt.DESeq2.DE_results"
sg_file="data/processed/DESeq2/DESeq2_gene/differentially_expressed_genes_with_symbols.txt"
gl_file="data/processed/DESeq2/DESeq2_gene/differentially_expressed_genes.txt"

# Import DGE results
awk 'NR>1 {print $0}' "$original_df_file" > temp_original_df.txt

# Import significant genes from blastx against SwissProt Metazoa
awk '!seen[$0]++' "$sg_file" > temp_sg.txt

# Import trinity gene list
gl_length=$(wc -l < "$gl_file")

# Select relevant columns and perform filtering
awk -v OFS="\t" -v log2fc_thresh=0.5 -v padj_thresh=0.01 'NR==1 {print "log2FoldChange", "padj", "diffexpressed"} 
    NR>1 {
        diffexpressed = "no"
        if ($1 > log2fc_thresh && $2 < padj_thresh) {
            diffexpressed = "up"
        } else if ($1 < -log2fc_thresh && $2 < padj_thresh) {
            diffexpressed = "down"
        }
        print $1, $2, diffexpressed
    }' temp_original_df.txt > filtered_df.txt

# Merge data frames and clean up
join -t $'\t' -1 1 -2 1 <(sort -k1,1 filtered_df.txt) <(sort -k1,1 temp_sg.txt) \
    | awk -v OFS="\t" '{print $1, $2, $3, ($4 == "." ? "NA" : $4)}' > merged_df.txt

# Extract unique blastp genes
cut -f4 merged_df.txt | grep -v "NA" | sort -u > blastp_genes.txt

# Clean up temporary files
rm temp_original_df.txt temp_sg.txt filtered_df.txt merged_df.txt

echo "Blastp genes extracted and saved in blastp_genes.txt"
