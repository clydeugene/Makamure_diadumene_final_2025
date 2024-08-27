#!/bin/bash

# Starting script
echo -e "\nThe script verticalHeatmap.sh has started in $(pwd) at $(date)\n"

# Read matrix from command line
matrix=$1

# Initialize the last_file variable
last_file=""

# Check if the name of the current directory has the words "gene" or "isoform" in it and process accordingly
if [[ $(pwd) == *"isoform"* ]]
then
	# Extract differentially expressed isoforms from the DESeq2 results
	awk -F "\t" -v OFS="\t" 'NR>1 {print $1}' $matrix | sort | uniq | awk -v header="isoform" 'BEGIN {print header} {print}' > differentially_expressed_isoforms.txt
	wait
	
	N=1
	for file in $(ls ../../../../data/processed/salmon/rt*/quant.sf)
	do
		echo Room_$N > trans_filtered_room_TPM_$N.txt
		awk -F "\t" '{if(NR>1){print $4}}' $file >> trans_filtered_room_TPM_$N.txt
		N=$((N+1))
	done

	N=1
	for file in $(ls ../../../../data/processed/salmon/16*/quant.sf)
	do
		echo Cold_$N > trans_filtered_cold_TPM_$N.txt
		awk -F "\t" '{if(NR>1){print $4}}' $file >> trans_filtered_cold_TPM_$N.txt
		N=$((N+1))
		last_file="$file"
	done

	awk -F "\t" '{if(NR==1){print "isoforms"}else{print $1}}' "$last_file" > Vertical_TPM_Labels.txt

	rm differentially_expressed*with_symbols.txt

	paste *_room_TPM_*.txt > room_TPM_ALL.txt
	paste *_cold_TPM_*.txt > cold_TPM_ALL.txt

	paste room_TPM_ALL.txt cold_TPM_ALL.txt > vertical_TPM_ALL.txt
	paste Vertical_TPM_Labels.txt vertical_TPM_ALL.txt > vertical_TPM_ALL+L_unspaced.txt

	awk -F "\t" -v OFS="\t" '{print $0}' vertical_TPM_ALL+L_unspaced.txt > vertical_TPM_ALL+L.txt

	# Sort the vertical_TPM_ALL+L.txt file by the first column
		header=$(head -n 1 vertical_TPM_ALL+L.txt)
		tail -n +2 vertical_TPM_ALL+L.txt | sort -k1,1b > vertical_TPM_ALL+L_sorted.txt
		echo -e "$header\n$(cat vertical_TPM_ALL+L_sorted.txt)" > vertical_TPM_ALL+L_sorted.txt

	join -t $'\t' --header differentially_expressed*.txt vertical_TPM_ALL+L_sorted.txt > vertical_TPM_DE.txt

	Rscript ../../../../scripts/verticalHeatmap.R vertical_TPM_DE.txt DE_TPM+LogFold_isoform

	echo -e "\nThe script verticalHeatmap.sh in $(pwd) has finished at $(date)\n"

elif [[ $(pwd) == *"gene"* ]]
then
	# Extract differentially expressed genes from the DESeq2 results
	awk -F "\t" -v OFS="\t" 'NR>1 {print $1}' $matrix | sort | uniq | awk -v header="gene" 'BEGIN {print header} {print}' > differentially_expressed_genes.txt
	wait
	
	N=1
	for file in $(ls ../../../../data/processed/salmon/rt*/quant.sf.genes) # Figure out a way to get use the isoform quant.sf files since this is based on the filtered fasta i.e., the longest isoform fasta
	do
		echo Room_$N > gene_filtered_room_TPM_$N.txt
		awk -F "\t" '{if(NR>1){print $4}}' $file >> gene_filtered_room_TPM_$N.txt
		N=$((N+1))
	done

	N=1
	for file in $(ls ../../../../data/processed/salmon/16*/quant.sf.genes)
	do
		echo Cold_$N > gene_filtered_cold_TPM_$N.txt
		awk -F "\t" '{if(NR>1){print $4}}' $file >> gene_filtered_cold_TPM_$N.txt
		N=$((N+1))
		last_file="$file"
	done

	awk -F "\t" '{if(NR==1){print "genes"}else{print $1}}' "$last_file" > Vertical_TPM_Labels.txt

	rm differentially_expressed*with_symbols.txt

	paste *_room_TPM_*.txt > room_TPM_ALL.txt
	paste *_cold_TPM_*.txt > cold_TPM_ALL.txt

	paste room_TPM_ALL.txt cold_TPM_ALL.txt > vertical_TPM_ALL.txt
	paste Vertical_TPM_Labels.txt vertical_TPM_ALL.txt > vertical_TPM_ALL+L_unspaced.txt

	awk -F "\t" -v OFS="\t" '{print $0}' vertical_TPM_ALL+L_unspaced.txt > vertical_TPM_ALL+L.txt

	# Sort the vertical_TPM_ALL+L.txt file by the first column
		header=$(head -n 1 vertical_TPM_ALL+L.txt)
		tail -n +2 vertical_TPM_ALL+L.txt | sort -k1,1b > vertical_TPM_ALL+L_sorted.txt
		echo -e "$header\n$(cat vertical_TPM_ALL+L_sorted.txt)" > vertical_TPM_ALL+L_sorted.txt

	join -t $'\t' --header differentially_expressed*.txt vertical_TPM_ALL+L_sorted.txt > vertical_TPM_DE.txt

	Rscript ../../../../scripts/verticalHeatmap.R vertical_TPM_DE.txt DE_TPM+LogFold_gene

	echo -e "\nThe script verticalHeatmap.sh in $(pwd) has finished at $(date)\n"

fi

# rm differentially_expressed*with_symbols.txt

# paste *_room_TPM_*.txt > room_TPM_ALL.txt
# paste *_cold_TPM_*.txt > cold_TPM_ALL.txt

# paste room_TPM_ALL.txt cold_TPM_ALL.txt > vertical_TPM_ALL.txt
# paste Vertical_TPM_Labels.txt vertical_TPM_ALL.txt > vertical_TPM_ALL+L_unspaced.txt

# awk -F "\t" -v OFS="\t" '{print $0}' vertical_TPM_ALL+L_unspaced.txt > vertical_TPM_ALL+L.txt

# # Sort the vertical_TPM_ALL+L.txt file by the first column
# 	header=$(head -n 1 vertical_TPM_ALL+L.txt)
# 	tail -n +2 vertical_TPM_ALL+L.txt | sort -k1,1b > vertical_TPM_ALL+L_sorted.txt
# 	echo -e "$header\n$(cat vertical_TPM_ALL+L_sorted.txt)" > vertical_TPM_ALL+L_sorted.txt

# join -t $'\t' --header differentially_expressed*.txt vertical_TPM_ALL+L_sorted.txt > vertical_TPM_DE.txt

# Rscript ../../../scripts/verticalHeatmap.R vertical_TPM_DE.txt DE_TPM+LogFold

# echo -e "\nThe script verticalHeatmap.sh in $(pwd) has finished at $(date)\n"
