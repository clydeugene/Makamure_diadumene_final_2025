#!/bin/bash

# Start of script
echo -e "\nThe script correlateReplicates.sh has started at $(date)\n"

# Define samples file
    samples=$1

# Declare an empty array to store the sample types
    sample_types=()

# Fill the array with unique values from the first column of the samples file
    while read condition rep path; do
        # Check if the value is not already in the array
        if [[ ! " ${sample_types[@]} " =~ " $condition " ]]; then
            # Add the value to the array
            sample_types+=("$condition")
        fi
    done < ../../../../$samples

# Define variables to hold the two sample types
    sample_type1=${sample_types[0]}
    sample_type2=${sample_types[1]}

if [[ $(pwd) == *"isoform"* ]]
then

    for type in ${sample_types[@]}
    do
        rm ${type}_samples.txt
        
        # Read conditions from sample file
        while read condition rep path
        do
            if [[ $condition == $type ]]
            then
                echo "${rep}" >> ${type}_samples.txt
            fi
        done < ../../../../$samples

        # Prepare the quant.sf files
        for rep in $(cat ${type}_samples.txt)
        do
            quants=../../salmon/${rep}/quant.sf
            awk 'NR > 1' $quants | cut -f1,4 | sort -k1,1 | cut -f2 > heatmap_DE.$rep.temp
        done
        wait
        
    done

    # Get the ids from a representative of each file to use as the row names
    quant_file1=../../salmon/${sample_type1}_rep_1/quant.sf
    awk 'NR > 1' $quant_file1 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.${sample_type1}_ids.temp

    quant_file2=../../salmon/${sample_type2}_rep_1/quant.sf
    awk 'NR > 1' $quant_file2 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.${sample_type2}_ids.temp

    # Combine the temp files
    paste $(ls heatmap_DE.*${sample_type1}*.temp) > heatmap_DE_${sample_type1}.temp2
    rm heatmap_DE.*${sample_type1}*.temp
    awk -v OFS='\t' -v type=${sample_type1} 'BEGIN{print "",type"_rep_1",type"_rep_2",type"_rep_3",type"_rep_4",type"_rep_5"} {print}' heatmap_DE_${sample_type1}.temp2 > heatmap_DE_${sample_type1}.temp3

    paste $(ls heatmap_DE.*${sample_type2}*.temp) > heatmap_DE_${sample_type2}.temp2
    rm heatmap_DE.*${sample_type2}*.temp
    awk -v OFS='\t' -v type=${sample_type2} 'BEGIN{print "",type"_rep_1",type"_rep_2",type"_rep_3",type"_rep_4",type"_rep_5"} {print}' heatmap_DE_${sample_type2}.temp2 > heatmap_DE_${sample_type2}.temp3

    # Create the correlation heatmaps
    Rscript ../../../../scripts/correlateReplicates.R heatmap_DE_${sample_type1}.temp3 ${sample_type1}
    wait
    Rscript ../../../../scripts/correlateReplicates.R heatmap_DE_${sample_type2}.temp3 ${sample_type2}
    wait

    Rscript ../../../../scripts/correlateReplicatesCombined.R heatmap_DE_${sample_type1}.temp3 heatmap_DE_${sample_type2}.temp3 isoform

elif [[ $(pwd) == *"gene"* ]]
then
    
    for type in ${sample_types[@]}
    do
        rm ${type}_samples.txt
        
        # Read conditions from sample file
        while read condition rep path
        do
            if [[ $condition == $type ]]
            then
                echo "${rep}" >> ${type}_samples.txt
            fi
        done < ../../../../$samples

        # Prepare the quant.sf files
        for rep in $(cat ${type}_samples.txt)
        do
            quants=../../salmon/${rep}/quant.sf.genes
            awk 'NR > 1' $quants | cut -f1,4 | sort -k1,1 | cut -f2 > heatmap_DE.$rep.temp
        done
        wait
        
    done

    # Get the ids from a representative of each file to use as the row names
    quant_file1=../../salmon/${sample_type1}_rep_1/quant.sf.genes
    awk 'NR > 1' $quant_file1 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.${sample_type1}_ids.temp

    quant_file2=../../salmon/${sample_type2}_rep_1/quant.sf.genes
    awk 'NR > 1' $quant_file2 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.${sample_type2}_ids.temp

    # Combine the temp files
    paste $(ls heatmap_DE.*${sample_type1}*.temp) > heatmap_DE_${sample_type1}.temp2
    rm heatmap_DE.*${sample_type1}*.temp
    awk -v OFS='\t' -v type=${sample_type1} 'BEGIN{print "",type"_rep_1",type"_rep_2",type"_rep_3",type"_rep_4",type"_rep_5"} {print}' heatmap_DE_${sample_type1}.temp2 > heatmap_DE_${sample_type1}.temp3

    paste $(ls heatmap_DE.*${sample_type2}*.temp) > heatmap_DE_${sample_type2}.temp2
    rm heatmap_DE.*${sample_type2}*.temp
    awk -v OFS='\t' -v type=${sample_type2} 'BEGIN{print "",type"_rep_1",type"_rep_2",type"_rep_3",type"_rep_4",type"_rep_5"} {print}' heatmap_DE_${sample_type2}.temp2 > heatmap_DE_${sample_type2}.temp3

    # Create the correlation heatmaps
    Rscript ../../../../scripts/correlateReplicates.R heatmap_DE_${sample_type1}.temp3 ${sample_type1}
    wait
    Rscript ../../../../scripts/correlateReplicates.R heatmap_DE_${sample_type2}.temp3 ${sample_type2}
    wait

    Rscript ../../../../scripts/correlateReplicatesCombined.R heatmap_DE_${sample_type1}.temp3 heatmap_DE_${sample_type2}.temp3 gene

    
    
    # rm rt_samples.txt 16_samples.txt

    # # Read conditions from samples.txt
    # while read line
    # do
    #     condition=$(echo $line | awk '{print $1}')
    #     if [[ $condition == "rt" ]]
    #     then
    #         echo $line | awk '{print $2}' >> rt_samples.txt
    #     elif [[ $condition == "16" ]]
    #     then
    #         echo $line | awk '{print $2}' >> 16_samples.txt
    #     fi
    # done < ../../raw/samples.txt

    # # Prepare the quant.sf files
    # for type in $(cat rt_samples.txt)
    # do
    #     File1=../salmon/${type}/quant.sf.genes
    #     awk 'NR > 1' $File1 | cut -f1,4 | sort -k1,1 | cut -f2 > heatmap_DE.$type.temp2
    # done
    # wait

    # for type in $(cat 16_samples.txt)
    # do
    #     File1=../salmon/${type}/quant.sf.genes
    #     awk 'NR > 1' $File1 | cut -f1,4 | sort -k1,1 | cut -f2 > heatmap_DE.$type.temp2
    # done
    # wait

    # File01=../salmon/Hal_1_rt.merged.SE_trimmed.fq/quant.sf.genes
    # awk 'NR > 1' $File01 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.Hal_0_rt.merged.SE_trimmed.fq.temp2

    # File02=../salmon/Hal_2_16.merged.SE_trimmed.fq/quant.sf.genes
    # awk 'NR > 1' $File02 | cut -f1,4 | sort -k1,1 | cut -f1 > heatmap_DE.Hal_0_16.merged.SE_trimmed.fq.temp2

    # paste $(ls heatmap_DE.*rt*.temp2) > heatmap_DE_rt.temp2
    # rm heatmap_DE.*rt*.temp2
    # awk -v OFS='\t' 'BEGIN{print "","Rep_RT_1","Rep_RT_3","Rep_RT_5","Rep_RT_6","Rep_RT_7"}{print}' heatmap_DE_rt.temp2 > heatmap_DE_rt.temp3
    # Rscript ../../../scripts/correlateReplicates.R heatmap_DE_rt.temp3 rt
    # mv heatmap_DE_rt.temp3 heatmap_DE_room.temp3

    # paste $(ls heatmap_DE.*16*.temp2) > heatmap_DE_16.temp2
    # rm heatmap_DE.*16*.temp2
    # awk -v OFS='\t' 'BEGIN{print "","Rep_16_2","Rep_16_3","Rep_16_4","Rep_16_5","Rep_16_6"}{print}' heatmap_DE_16.temp2 > heatmap_DE_16.temp3
    # Rscript ../../../scripts/correlateReplicates.R heatmap_DE_16.temp3 16
    # mv heatmap_DE_16.temp3 heatmap_DE_cold.temp3
    
    # Rscript ../../../scripts/correlateReplicatesCombined.R heatmap_DE_cold.temp3 heatmap_DE_room.temp3 gene

fi

# End of script
echo -e "\nThe script correlateReplicates.sh has finished at $(date)\n"
