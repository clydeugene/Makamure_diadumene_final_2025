#!/bin/bash

# This script will pool all the bam files in the current directory and create a single bam file for IGV visualization

# Starting script
echo -e "\nThe script poolForIGV.sh has started at $(date)\n"

# Merge the bam files
echo -e "\nMerging the bam files\n"

for i in 16 rt
do
    samtools merge ${i}_pooled.bam ${i}_rep_1.nSorted.star.hali_transcriptome.fasta.bam ${i}_rep_2.nSorted.star.hali_transcriptome.fasta.bam ${i}_rep_3.nSorted.star.hali_transcriptome.fasta.bam ${i}_rep_4.nSorted.star.hali_transcriptome.fasta.bam ${i}_rep_5.nSorted.star.hali_transcriptome.fasta.bam &
done
wait

# Convert the bam files to bed files
for i in 16 rt
do
    bamToBed -split -i ${i}_pooled.bam > ${i}_pooled.bed &
done
wait


# Make a genome size file using the method at https://www.biostars.org/p/173963/#174150
echo -e "\nMaking a genome size file\n"
samtools faidx hali_transcriptome.fasta
cut -f1,2 hali_transcriptome.fasta.fai > hali_transcriptome.genome


# Make bedGraph files
for i in 16 rt
do
    scale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c ${i}_pooled.bam)")
    cut -f 1-3 ${i}_pooled.bed | sort -k1,1 | bedtools genomecov -scale $scale -i - -g hali_transcriptome.genome -bg > ${i}_pooled_reads.bedGraph &
done
wait

# Convert the bedGraph files to bigWig files
for i in 16 rt
do
    bedGraphToBigWig ${i}_pooled_reads.bedGraph hali_transcriptome.genome ${i}_pooled_reads.bw &
done
wait

# Remove the intermediate files
rm *_pooled_reads.bedGraph

# Compress the bed files
gzip *_pooled.bed

# End of script
echo -e "\nThe script poolForIGV.sh has ended at $(date)\n"


# $ TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c rt_pooled.bam)") && echo $TmpScale