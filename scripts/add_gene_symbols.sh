#!/bin/bash

# Start of script
      echo -e "\nThe script add_gene_symbols.sh has started in $(pwd) at $(date)\n"

# Set variables
        # diffexpr="$1"
        level="$1"
        blastout="$2"

if [[ $(pwd) == *"isoform"* ]]
then
	# Prepare the blast output file. Sorting order: qseqid, evalue, bitscore, pident. Lastly extract gene symbol from the top ranking qseqid
        sort -k1,1 -k11,11n -k12,12nr -k3,3nr ../../../../$blastout | awk -F "\t" -v OFS="\t" '{split($1, a, ".p"); split($2, b, "|"); split(b[3], c, "_HU"); print c[1], a[1]}' | uniq -f1 | awk -v OFS="\t" '{print $2,$1}' > blastout_filtered.temp

	# Prepare the differentially expressed genes file
        tail -n +2 differentially_expressed_${level}s.txt | sort -k1b,1 > diffexpr_filtered.temp

        # Add gene symbols to the differentially expressed genes
        join -t $'\t' diffexpr_filtered.temp blastout_filtered.temp | awk -v OFS="\t" 'BEGIN{print "trinity_id","top_blastp_hit"} {print $0}' > differentially_expressed_isoforms_with_symbols.txt

elif [[ $(pwd) == *"gene"* ]]
then
	# Prepare the blast output file. Sorting order: qseqid, evalue, bitscore, pident. Lastly extract gene symbol from the top ranking qseqid
        sort -k1,1 -k11,11n -k12,12nr -k3,3nr ../../../../$blastout | awk -F "\t" -v OFS="\t" '{split($1, a, "_i"); split($2, b, "|"); split(b[3], c, "_HU"); print c[1], a[1]}' | uniq -f1 | awk -v OFS="\t" '{print $2,$1}' > blastout_filtered.temp

        # Prepare the differentially expressed genes file
        tail -n +2 differentially_expressed_${level}s.txt | sort -k1b,1 > diffexpr_filtered.temp
        
        # Add gene symbols to the differentially expressed genes
        join -t $'\t' diffexpr_filtered.temp blastout_filtered.temp | awk -v OFS="\t" 'BEGIN{print "trinity_id","top_blastp_hit"} {print $0}' > differentially_expressed_genes_with_symbols.txt

fi

# End of script
      echo -e "\nThe script add_gene_symbols.sh in $(pwd) has ended at $(date)\n"        
