#!/usr/bin/env bash

# This script selects only those contigs that had hits on the nr database indicating protein coding potential

# Define the input arguments
hits=$1
fasta=$2
output=$2

# Extract the contigs that had hits on the nr database
echo -e "\nExtracting the contigs that had hits on the nr database"

awk '{print $1}' $hits | sort | uniq > temp1.txt

# Linearize the fasta file
echo -e "\nLinearizing the fasta file"

awk -f code/linearizefasta.awk $fasta > temp2.txt

# Extract the contigs that had hits on the nr database
echo -e "\nExtracting the contigs that had hits on the nr database"

grep -f temp1.txt temp2.txt > temp3.txt

tr "\t" "\n" < temp3.txt > $output
wait

# Remove the temporary files
echo -e "\nRemoving the temporary files"

# rm temp1.txt temp2.txt

# End of script
echo -e "\nExtracting the contigs that had hits on the nr database finished\n"