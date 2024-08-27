#!/bin/bash

# This script replaces trinity IDs in the fasta file with gene symbols from the blast output file

# Start of script
echo -e "\nThe script add_gene_symbols_to_blast_output.sh has started in $(pwd) at $(date)\n"

# If blastout_ files or tempfile exist, remove them
if [ -f blastout_* ]; then
    rm blastout_*
fi

if [ -f tempfile ]; then
    rm tempfile
fi

# Linearize the fasta file for easier processing
awk -f code/linearizefasta.awk data/processed/hali_transcriptome.nr.fasta > tempfile

# Define input files
fasta_file="tempfile"
blastout_file="data/processed/DESeq2/DESeq2_gene/blastout_filtered.temp"

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
        sed -i "s/$trinity_id/$trinity_id $gene/g" "$fasta_file"
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
tr "\t" "\n" < "$fasta_file" > data/processed/hali_transcriptome.nr_with_symbols.fasta

# Remove the split files and the temporary file
rm blastout_*
rm tempfile

# End of script
echo -e "\nThe script add_gene_symbols_to_blast_output.sh in $(pwd) has ended at $(date)\n"