#!/bin/bash

# Start of script
echo -e "\n*******************************\nThe script doAlignments_to_Diadumene_Genome.sh has started on $(date)\n*******************************\n"

# Read in the command line arguments
refGenome=$1

# Make the analysis directory if it doesn't exist
mkdir -p data/processed/gatk/genome
cd data/processed/gatk/genome

# Make a logs directory if it doesn't exist
mkdir -p logs

# Starting step 1
echo -e "\n**************************************************************\nStep 1 Started - Aligning trimmed reads to the reference genome\n**************************************************************\n"

# Starting genome indexing
echo -e "\n*******************************\nIndexing the reference genome has started on $(date)\n*******************************\n"

# Index the reference genome
bwa index ../../../../$refGenome
samtools faidx ../../../../$refGenome
gatk CreateSequenceDictionary -R ../../../../$refGenome

# Starting alignment
echo -e "\n*******************************\nAligning the trimmed reads to the reference genome\n*******************************\n"

# Create an alignment function
align_reads() {
    local read=$1

    # Get the base name of the file
    local base=$(basename $read .merged.SE_trimmed.fq.gz)

    echo -e "*******************************\nAligning $base started - $(date)\n*******************************\n"

    # Align the reads to the Diadumene lineata genome
    bwa mem -t 5 ../../../../$refGenome $read | samtools view -@ 5 -Sb -o $base.bam

    echo -e "\n*******************************\nDone aligning $base - $(date)\n*******************************\n"
}

# Start a counter
i=1

# Run the following command to align the trimmed reads to the Diadumene lineata genome using BWA
for read in $(ls ../../../raw/*merged.SE_trimmed.fq.gz)
do
    # Check if the counter is divisible by 40 to prevent overloading the system
    if [ $(echo $i | awk '{print $1%40}') -eq 0 ]; then
        wait
    fi

    # Define a log file
    logfile="logs/$(basename $read .merged.SE_trimmed.fq.gz).align_reads.log"

    # Align the reads to the Diadumene lineata genome
    align_reads $read > $logfile 2>&1 &
    i=$(expr $i + 1)
done
wait

# End of step 1
echo -e "\n**************************************************************\nStep 1 Completed - Done aligning trimmed reads to the reference genome\n**************************************************************\n"

# Starting step 2
echo -e "\n**************************************************************\nStep 2 Started - Variant calling\n**************************************************************\n"

# Report current step
echo -e "\n*******************************\nRemoving soft clipped bases from bam files - $(date)\n*******************************\n"

# Create a bam clipping function
remove_soft_clip_bam() {
    local bam=$1
    local base=$(basename $bam .bam)

    echo -e "\n*******************************\nStarting removing soft clipped bases for $base - $(date)\n*******************************\n"

    # Sort the bam file by coordinates
    samtools sort -@ 5 -o $base.sorted.bam $bam

    # Clean and convert bam file with Picard tools
    picard AddOrReplaceReadGroups I=$base.sorted.bam O=$base.sorted.rg.bam VALIDATION_STRINGENCY=SILENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

    # Remove soft-clipped bases
    jvarkit biostar84452 -o $base.cliprmvd.bam --samoutputformat BAM $base.sorted.rg.bam

    # Index the cleaned bam file
    samtools index $base.cliprmvd.bam


    # Echo done
    echo -e "\n*******************************\nDone removing soft clipped bases for $base - $(date)\n*******************************\n"
}

# Restart the counter
i=1

# Run the following pipeline to remove soft-clipped bases from the bam files
for bam in $(ls *.bam)
do
    # Check if the counter is divisible by 40 to prevent overloading the system
    if [ $(echo $i | awk '{print $1%40}') -eq 0 ]; then
        wait
    fi

    # Define a log file
    logfile="logs/$(basename $bam .bam).variant_calling.log"

    # Call variants
    remove_soft_clip_bam $bam >> $logfile 2>&1 &
    i=$(expr $i + 1)
    
done
wait

# Done removing soft-clipped bases
echo -e "\n*******************************\nDone removing soft-clipped bases from bam files - $(date)\n*******************************\n"


# If there are any cliprmvd bam files in the current directory, move them to the appropriate directories
if ls *cliprmvd.bam 1> /dev/null 2>&1; then
    # Moving bam files
    echo -e "\n*******************************\nMoving bam files to appropriate directories - $(date)\n*******************************\n"

    # Make directories for the bam files if they don't exist
    mkdir -p gatk_rt gatk_16

    # Move bam files to the appropriate directories
    mv *rt*cliprmvd.ba* gatk_rt
    mv *16*cliprmvd.ba* gatk_16

    # Remove intermediate files
    rm *sorted.*ba*
    
    # Report done
    echo -e "\n*******************************\nDone moving bam files to appropriate directories - $(date)\n*******************************\n"
fi

# Report current step
echo -e "\n*******************************\nDownsampling the bam files with soft-clipped bases removed - $(date)\n*******************************\n"

# Remove reads_clipped.tsv file if it exists
if [ -f reads_clipped.tsv ]; then
    rm reads_clipped.tsv
fi

# Create a reads count file
for bam in $(ls gatk_*/*.cliprmvd.bam)
do
    sample=$(basename $bam .cliprmvd.bam)
    reads=$(samtools view -c $bam)

    echo -e "$sample\t$reads" >> reads_clipped.tsv
done

# Get the smallest read count from the reads count file
smallest_read_count=$(awk 'NR==1 || $2 < min {min = $2} END {print min}' reads_clipped.tsv)

# Report the smallest read count
echo -e "\n*******************************\nDownsampling the bam files to the smallest read count: $smallest_read_count - $(date)\n*******************************\n"

# Create a downsampling function
downsample_bam() {
    local bam=$1
    local smallest_read_count=$2
    local base=$(basename $bam .cliprmvd.bam)
    local out_dir=$(dirname $bam)

    echo -e "\n*******************************\nStarting downsampling for $base - $(date)\n*******************************\n"

    # Downsample the bam file using http://lindenb.github.io/jvarkit/Biostar145820.html
    jvarkit biostar145820 -n $smallest_read_count -o $out_dir/$base.downsampled.bam $bam
    wait

    echo -e "\n*******************************\nDone downsampling for $base - now sorting bam file - $(date)\n*******************************\n"

    # Sort the downsampled bam file
    samtools sort -@ 5 -o $out_dir/$base.downsampled.sorted.bam $out_dir/$base.downsampled.bam

    echo -e "\n*******************************\nDone sorting bam for $base - now separating mapped and unmapped reads - $(date)\n*******************************\n"

    # Separate mapped reads from unmapped reads
    samtools view -@ 5 -b -F 4 $out_dir/$base.downsampled.sorted.bam > $out_dir/$base.downsampled.mapped.bam
    samtools view -@ 5 -b -f 4 $out_dir/$base.downsampled.sorted.bam > $out_dir/$base.downsampled.unmapped.bam

    echo -e "\n*******************************\nDone separating mapped and unmapped reads for $base - now merging the mapped and unmapped reads - $(date)\n*******************************\n"

    # Merge the mapped and unmapped reads
    samtools merge -@ 5 -f -c $out_dir/$base.downsampled.bam $out_dir/$base.downsampled.mapped.bam $out_dir/$base.downsampled.unmapped.bam

    echo -e "\n*******************************\nDone merging mapped and unmapped reads for $base - now indexing the downsampled bam file - $(date)\n*******************************\n"

    # Index the downsampled bam file
    samtools index $out_dir/$base.downsampled.bam

    # Echo done
    echo -e "\n*******************************\nDone downsampling steps for $base - $(date)\n*******************************\n"
}

# Restart the counter
i=1

# Run the following pipeline to downsample the bam files
for bam in $(ls gatk_*/*.cliprmvd.bam)
do
    # Check if the counter is divisible by 40 to prevent overloading the system
    if [ $(echo $i | awk '{print $1%40}') -eq 0 ]; then
        wait
    fi

    # Define a log file
    logfile="logs/$(basename $bam .cliprmvd.bam).downsampling.log"

    # Downsample the bam files
    downsample_bam $bam $smallest_read_count > $logfile 2>&1 &
    i=$(expr $i + 1)
    
done
wait

# End of downsampling
echo -e "\n*******************************\nDone downsampling the bam files - $(date)\n*******************************\n"

# Report current step
echo -e "\n*******************************\nStarting somatic variant calling - $(date)\n*******************************\n"

# Create a somatic variant calling function
run_mutect2() {
    local bam=$1
    local base=$(basename $bam .downsampled.bam)
    local out_dir=$(dirname $bam)

    echo -e "\n*******************************\nStarting somatic variant calling for $base - $(date)\n*******************************\n"

    # # # Run splitNCigarReads
    # # gatk SplitNCigarReads -R ../../../../$refGenome -I $out_dir/$base.sorted.rg.bam -O $out_dir/$base.sorted.rg.split.bam --read-validation-stringency LENIENT

    # # Run the HaplotypeCaller tool to call variants
    # gatk HaplotypeCaller -R ../../../../$refGenome -I $out_dir/$base.cliprmvd.bam -O $out_dir/$base.vcf

    # # Run the VariantFiltration tool to filter the variants
    # gatk VariantFiltration -R ../../../../$refGenome -V $out_dir/$base.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O $out_dir/$base.filtered.vcf

    # Run Mutect2 to call somatic variants
    gatk Mutect2 -R ../../../../$refGenome -I $bam -O $out_dir/$base.somatic.vcf
    wait

    # Run the FilterMutectCalls tool to filter the somatic variants
    gatk FilterMutectCalls -R ../../../../$refGenome -V $out_dir/$base.somatic.vcf -O $out_dir/$base.somatic.filtered.vcf

    # Echo done
    echo -e "\n*******************************\nDone variant calling for $base - $(date)\n*******************************\n"
}

# Restart the counter
i=1

# Run the following pipeline to call somatic variants
for bam in $(ls gatk_*/*.downsampled.bam)
do
    # Check if the counter is divisible by 40 to prevent overloading the system
    if [ $(echo $i | awk '{print $1%40}') -eq 0 ]; then
        wait
    fi

    # Define a log file
    logfile="logs/$(basename $bam .downsampled.bam).somatic_variant_calling.log"

    # Call somatic variants
    run_mutect2 $bam > $logfile 2>&1 &
    i=$(expr $i + 1)
    
done
wait

# Done with somatic variant calling
echo -e "\n*******************************\nDone with somatic variant calling - $(date)\n*******************************\n"

# Separate SNPs and INDELs
echo -e "\n*******************************\nSeparating SNPs and INDELs - $(date)\n*******************************\n"

# Create a separate_snps_indels function
separate_snps_indels() {
    local vcf=$1
    local base=$(basename $vcf .vcf)
    local out_dir=$(dirname $vcf)

    echo -e "\n*******************************\nStarting separating SNPs and INDELs for $base - $(date)\n*******************************\n"

    # Separate SNPs and INDELs
    picard SplitVcfs I=$vcf SNP_OUTPUT=$out_dir/$base.snps.vcf INDEL_OUTPUT=$out_dir/$base.indels.vcf STRICT=false

    # Echo done
    echo -e "\n*******************************\nDone separating SNPs and INDELs for $base - $(date)\n*******************************\n"
}

# Run the following loop to separate SNPs and INDELs
for vcf in $(ls gatk_*/*.somatic.filtered.vcf)
do
    # Define a log file
    logfile="logs/$(basename $vcf .somatic.filtered.vcf).separate_snps_indels.log"

    # Separate SNPs and INDELs
    separate_snps_indels $vcf > $logfile 2>&1
done
wait

# Done separating SNPs and INDELs
echo -e "\n*******************************\nDone separating SNPs and INDELs - $(date)\n*******************************\n"

# # End of step 2
# echo -e "\n**************************************************************\nStep 2 Completed - Done variant calling\n**************************************************************\n"

# Starting step 3
echo -e "\n**************************************************************\nStep 3 Started - snp statistics\n**************************************************************\n"


# Report current step
echo -e "\n*******************************\nStarting snp count for the vcf files - $(date)\n*******************************\n"

# Remove count files if they exist
if [ -f snp_and_indel_counts.tsv ]; then
    rm snp_and_indel_counts.tsv
fi

if [ -f snp_counts.tsv ]; then
    rm snp_counts.tsv
fi

if [ -f indel_counts.tsv ]; then
    rm indel_counts.tsv
fi

rm snp_and_indel_counts_no_pass.tsv snp_counts_no_pass.tsv indel_counts_no_pass.tsv

# Create a snp count function
snp_counts() {
    local vcf=$1
    local base=$(basename $vcf .somatic.filtered.vcf)
    local dir_name=$(dirname $vcf)

    echo -e "\n*******************************\nStarting snp count for $base - $(date)\n*******************************\n"

    # Count the number of snps in the vcf file and append the count to snp_and_indel_counts.tsv
    snp_indel_count=$(grep -v "^#" $vcf | grep "PASS" | wc -l)
    echo -e "$base\t$snp_indel_count" >> snp_and_indel_counts.tsv
    snp_indel_count2=$(grep -v "^#" $vcf | wc -l)
    echo -e "$base\t$snp_indel_count2" >> snp_and_indel_counts_no_pass.tsv

    # Count the number of snps in the vcf file and append the count to snp_counts.tsv
    snp_count=$(grep -v "^#" $dir_name/$base.somatic.filtered.snps.vcf | grep "PASS" | wc -l)
    echo -e "$base\t$snp_count" >> snp_counts.tsv
    snp_count2=$(grep -v "^#" $dir_name/$base.somatic.filtered.snps.vcf | wc -l)
    echo -e "$base\t$snp_count2" >> snp_counts_no_pass.tsv

    # Count the number of indels in the vcf file and append the count to indel_counts.tsv
    indel_count=$(grep -v "^#" $dir_name/$base.somatic.filtered.indels.vcf | grep "PASS" | wc -l)
    echo -e "$base\t$indel_count" >> indel_counts.tsv
    indel_count2=$(grep -v "^#" $dir_name/$base.somatic.filtered.indels.vcf | wc -l)
    echo -e "$base\t$indel_count2" >> indel_counts_no_pass.tsv

    # Echo done
    echo -e "\n*******************************\nDone snp count for $base - $(date)\n*******************************\n"
}

# Run the following loop to count the number of snps in the vcf files
for vcf in $(ls gatk_*/*.somatic.filtered.vcf)
do
    # Define a log file
    logfile="logs/all_snp_counts.log"

    # Count the number of snps in the vcf file
    snp_counts $vcf >> $logfile 2>&1
done
wait

# Report done
echo -e "\n*******************************\nDone snp count for the vcf files - $(date)\n*******************************\n"

# End of step 3
echo -e "\n**************************************************************\nStep 3 Completed - Done snp statistics\n**************************************************************\n"

# End of script
echo -e "\n*******************************\nThe script doAlignments_to_Diadumene_Genome.sh has ended on $(date)\n*******************************\n"

