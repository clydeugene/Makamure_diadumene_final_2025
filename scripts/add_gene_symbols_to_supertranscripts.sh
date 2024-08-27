#!/bin/bash

# This script replaces trinity IDs in the supertranscripts fasta file with gene symbols from the blast output file

# Start of script
echo -e "\nThe script add_gene_symbols_to_supertranscripts.sh has started in $(pwd) at $(date)\n"

# If blastout_ files or tempfile exist, remove them
if [ -f blastout_* ]; then
    rm blastout_*
fi

if [ -f tempfile ]; then
    rm tempfile
fi

# Linearize the fasta file for easier processing
awk -f code/linearizefasta.awk data/processed/dexseq/hali_transcriptome.fasta > tempfile
awk -F "\t" '{print $1}' data/processed/dexseq/hali_transcriptome.gtf > tempfile2


# Define input files
fasta_file="tempfile"
blastout_file="data/processed/DESeq2/DESeq2_gene/blastout_filtered.temp"
gtf_file="tempfile2"

# split the fasta file
split -l 500 $blastout_file blastout_
echo -e "\nThe $blastout_file file has been split into $(ls blastout_* | wc -l) files\n Starting to process each split file\n"

# Set counter
counter=1

# Define a function to process each blastout file
process_file() {
    local file=$1
    local count=$2
    while read trinity_id gene; do
        sed -i "s/$trinity_id/$trinity_id\_$gene/g" "$fasta_file"
        sed -i "s/$trinity_id/$trinity_id\_$gene/g" "$gtf_file"
    done < "$file"
}

# Iterate over blastout files and run process_file function in background
for file in blastout_*; do
    echo "Processing file $counter of $(ls blastout_* | wc -l)"
    process_file "$file" "$counter"
    let counter++
done

echo -e "\nReplacement complete\n"

# Replace the original fasta file with the updated one
tr "\t" "\n" < "$fasta_file" > data/processed/dexseq/hali_transcriptome_with_symbols.fasta

# Remove the split files and the temporary file
rm blastout_*
rm tempfile

# Add the gene symbols column to the supertranscripts gtf file
paste tempfile2 data/processed/dexseq/hali_transcriptome.gtf > tempfile3

# Replace the original gene symbols column
# awk -F "\t" -v OFS="\t" '{print $1, $3, $4, $5, $6, $7, $8, $9, $10, $1, $12, $13}' tempfile3 > data/processed/dexseq/hali_transcriptome_with_symbols.gtf
awk -v OFS="\t" '$1~"TRINITY" {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\"%s\";\t%s\t%s\n", $1, $3, $4, $5, $6, $7, $8, $9, $10, $1, $12, $13}' tempfile3 > data/processed/dexseq/hali_transcriptome_with_symbols.gtf

# # Make an array of the top 10 dexseq trinity IDs
# top10=("TRINITY_DN140401_c5_g2" "TRINITY_DN125768_c0_g1" "TRINITY_DN140401_c5_g2" "TRINITY_DN137709_c2_g1" "TRINITY_DN140401_c5_g2" "TRINITY_DN137519_c1_g1" "TRINITY_DN137235_c1_g1" "TRINITY_DN142284_c7_g2" "TRINITY_DN142417_c5_g1" "TRINITY_DN140401_c5_g2")

# For those genes with a gene symbol, add the gene symbol to the isoform name
awk -v OFS="\t" '{if ($1~"_g1_") {split($12, a, "_i"); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\%s\t%s\t\"%s\"\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $1"_i"a[2]} else {print $0}}' data/processed/dexseq/hali_transcriptome_with_symbols.gtf > data/processed/dexseq/hali_transcriptome_with_symbols_isoform.gtf


# End of script
echo -e "\nThe script add_gene_symbols_to_blast_output.sh in $(pwd) has ended at $(date)\n"