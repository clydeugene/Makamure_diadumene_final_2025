#!/bin/bash

# Starting makeNormalizedHeatmap.sh
echo -e "\nThe script makeNormalizedHeatmap.sh has started in $(pwd) at $(date)\n"

fdr=$1
foldchange=$2

# Initialize the heatmap name variable
heatmap_name=$(echo "$(basename $(pwd))_heatmap_combined")
heatmap_TMM_name=$(echo "$(basename $(pwd))_heatmap_TMM_combined")

# Check if the name of the current directory has the word "gene" or "isoform" in it and process accordingly
if [[ $(pwd) == *"isoform"* ]]
then
        # Extract the significant isoforms from the DESeq2 results
        join <(awk 'NR>1 {print $1}' differentially_expressed_isoforms.txt) <(awk -F "\t" -v OFS="\t" 'NR>1 {print $0}' trinity*.DE_results | sort -k1b,1) -o 1.1,2.7 > significant_DE_results.txt
        wait

        # Exit if significant_DE_results.txt is not created
        if [[ ! -s significant_DE_results.txt ]]
        then
                echo -e "\nERROR: significant_DE_results.txt is empty or not created. Exiting...\n"
                exit 1
        fi

        awk -F "\t" -v OFS="\t" '{print $1}' significant_DE_results.txt > significant_isoforms.txt

        split -l 500 significant_isoforms.txt significant_isoformsa

        process() {
        for name in $(awk '{print $1}' $File)
        do
                awk -F "\t" -v OFS="\t" -v name="$name" '{if($1==name){print NR, $0}}' Hclust_TPMTable_2.txt >> Hclust_TPMTable_2_sigs_$File
        done
        }

        for File in $(ls significant_isoformsa*)
        do
                process &
        done

        wait

        cat Hclust_TPMTable_2_sigs_significant_isoformsa* > Hclust_TPMTable_2_sigs.txt
        rm Hclust_TPMTable_2_sigs_significant_isoformsa*
        rm significant_isoformsa*

        sort -k1n,1 Hclust_TPMTable_2_sigs.txt | awk -F "\t" -v OFS="\t" '{print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' > Hclust_TPMTable_3_sigs_notop.txt

        head -n 1 Hclust_TPMTable_2.txt > Hclust_TPMTable_3_sigs.txt
        cat Hclust_TPMTable_3_sigs_notop.txt >> Hclust_TPMTable_3_sigs.txt

        ###

        awk  -F "\t" -v OFS="\t" 'NR>1{print$1}' Hclust_TPMTable_3_sigs.txt | sort -k1b,1 > Isoforms.txt

        awk  -F "\t" -v OFS="\t" 'BEGIN{print "Isoforms", "RT_vs_16"}{print}' significant_DE_results.txt > Isoform_logfolds.txt

        head -n 1 Isoform_logfolds.txt | awk -F "\t" -v OFS="\t" '{gsub(" ","\t");print}' > Isoform_logfolds.sorted.mod2

        join Isoforms.txt Isoform_logfolds.txt | awk -F "\t" -v OFS="\t" '{gsub(" ","\t");print}' >> Isoform_logfolds.sorted.mod2

        awk -F "\t" -v OFS="\t" '{print $1, $2, $2}' Isoform_logfolds.sorted.mod2 > Isoform_logfolds.sorted.mod

        ###

        awk  -F "\t" -v OFS="\t" 'NR>1{print$1}' Hclust_TPMTable_3_sigs.txt > Isoforms.txt

        split -l 500 Isoforms.txt Isoformsa

        process() {
        for name in $(awk '{print $1}' $File)
        do
                awk -F "\t" -v OFS="\t" -v name="$name" '{if($1==name){print}}' vertical_TPM_DE.txt >> vertical_TPM_DE_$File
        done
        }

        for File in $(ls Isoformsa*)
        do
                process &
        done

        wait

        cat vertical_TPM_DE_Isoformsa* > vertical_TPM_DE_sigs.txt
        rm vertical_TPM_DE_Isoformsa*
        rm Isoformsa*

        head -n 1 Hclust_TPMTable_2.txt > vertical_TPM_DE_sigs2.txt
        cat vertical_TPM_DE_sigs.txt >> vertical_TPM_DE_sigs2.txt

        #Sorted by log2fold
        head -n 1 Isoform_logfolds.sorted.mod > Isoform_logfolds.sorted.mod2
        tail -n +2 Isoform_logfolds.sorted.mod | sort -k2nr,2 > Isoform_logfolds.sorted.mod2
        sort -k2,2nr significant_DE_results.txt | awk -v OFS="\t" '{print $1,$2,$2}' > sorted_logfolds.txt
        Rscript ../../../../scripts/makeNormalizedHeatmap.R Isoform_logfolds.sorted.mod2 vertical_TPM_DE_sigs2.txt $heatmap_name
        Rscript ../../../../scripts/makeNormalizedHeatmapTMM.R sorted_logfolds.txt diffExpr.P${fdr}_C${foldchange}.matrix $heatmap_TMM_name

elif [[ $(pwd) == *"gene"* ]]
then    
        # Extract the significant genes from the DESeq2 results
        join <(awk 'NR>1 {print $1}' differentially_expressed_genes.txt) <(awk -F "\t" -v OFS="\t" 'NR>1 {print $0}' trinity*.DE_results | sort -k1b,1) -o 1.1,2.7 > significant_DE_results.txt
        wait

        # Exit if significant_DE_results.txt is not created
        if [[ ! -s significant_DE_results.txt ]]
        then
                echo -e "\nERROR: significant_DE_results.txt is empty or not created. Exiting...\n"
                exit 1
        fi
        
        awk -F "\t" -v OFS="\t" '{print $1}' significant_DE_results.txt > significant_genes.txt

        split -l 500 significant_genes.txt significant_genesa

        process() {
        for name in $(awk '{print $1}' $File)
        do
                awk -F "\t" -v OFS="\t" -v name="$name" '{if($1==name){print NR, $0}}' Hclust_TPMTable_2.txt >> Hclust_TPMTable_2_sigs_$File
        done
        }

        for File in $(ls significant_genesa*)
        do
                process &
        done

        wait

        cat Hclust_TPMTable_2_sigs_significant_genesa* > Hclust_TPMTable_2_sigs.txt
        rm Hclust_TPMTable_2_sigs_significant_genesa*
        rm significant_genesa*

        sort -k1n,1 Hclust_TPMTable_2_sigs.txt | awk -F "\t" -v OFS="\t" '{print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' > Hclust_TPMTable_3_sigs_notop.txt

        head -n 1 Hclust_TPMTable_2.txt > Hclust_TPMTable_3_sigs.txt
        cat Hclust_TPMTable_3_sigs_notop.txt >> Hclust_TPMTable_3_sigs.txt

        ###

        awk  -F "\t" -v OFS="\t" 'NR>1{print$1}' Hclust_TPMTable_3_sigs.txt | sort -k1b,1 > Genes.txt

        awk  -F "\t" -v OFS="\t" 'BEGIN{print "Genes", "RT_vs_16"}{print}' significant_DE_results.txt > Gene_logfolds.txt

        head -n 1 Gene_logfolds.txt | awk -F "\t" -v OFS="\t" '{gsub(" ","\t");print}' > Gene_logfolds.sorted.mod2

        join Genes.txt Gene_logfolds.txt | awk -F "\t" -v OFS="\t" '{gsub(" ","\t");print}' >> Gene_logfolds.sorted.mod2

        awk -F "\t" -v OFS="\t" '{print $1, $2, $2}' Gene_logfolds.sorted.mod2 > Gene_logfolds.sorted.mod

        ###

        awk  -F "\t" -v OFS="\t" 'NR>1{print$1}' Hclust_TPMTable_3_sigs.txt > Genes.txt

        split -l 500 Genes.txt Genesa

        process() {
        for name in $(awk '{print $1}' $File)
        do
                awk -F "\t" -v OFS="\t" -v name="$name" '{if($1==name){print}}' vertical_TPM_DE.txt >> vertical_TPM_DE_$File
        done
        }

        for File in $(ls Genesa*)
        do
                process &
        done

        wait

        cat vertical_TPM_DE_Genesa* > vertical_TPM_DE_sigs.txt
        rm vertical_TPM_DE_Genesa*
        rm Genesa*

        head -n 1 Hclust_TPMTable_2.txt > vertical_TPM_DE_sigs2.txt
        cat vertical_TPM_DE_sigs.txt >> vertical_TPM_DE_sigs2.txt

        #Sorted by log2fold
        head -n 1 Gene_logfolds.sorted.mod > Gene_logfolds.sorted.mod2
        tail -n +2 Gene_logfolds.sorted.mod | sort -k2nr,2 > Gene_logfolds.sorted.mod2
        sort -k2,2nr significant_DE_results.txt | awk -v OFS="\t" '{print $1,$2,$2}' > sorted_logfolds.txt
        Rscript ../../../../scripts/makeNormalizedHeatmap.R Gene_logfolds.sorted.mod2 vertical_TPM_DE_sigs2.txt $heatmap_name
        Rscript ../../../../scripts/makeNormalizedHeatmapTMM.R sorted_logfolds.txt diffExpr.P${fdr}_C${foldchange}.matrix $heatmap_TMM_name

fi
wait

# End of script
echo -e "\nThe script makeNormalizedHeatmap.sh in $(pwd) has finished at $(date)\n"