#!/usr/bin/env bash

# This script filters a fasta file by contig length then selects only the longest isoforms. It takes three arguments:
# 1. The input fasta file
# 2. The name of the output file
# 3. The minimum contig length to keep

# Define the input arguments
fasta=$1
filtered=$2
contiglength=$3

# Select only those contigs that are longer than the minimum length
echo -e "\nFiltering $fasta by minimum contig length of $contiglength"

awk -f code/linearizefasta.awk $fasta | awk '{gsub(/n=/, " "); print $0}' | awk -v contiglength=$contiglength '$3>=contiglength' | sort -k1,1 | awk '{gsub(/le /, "len="); print $0}' | tr "\t" "\n" > data/raw/tempfile.txt
wait

mv data/raw/tempfile.txt $filtered
wait

# Script end
echo -e "\nScript complete"