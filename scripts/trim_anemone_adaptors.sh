#Stephen Justinen
echo trim_anemone_adaptors.sh start
process ()
{
	trim_galore $file
}
for file in $(ls Hal*.merged.SE.fastq.gz)
do
	process &
done
wait
echo trim_anemone_adaptors.sh done


# trim_galore version 0.6.6

# Add trim_galore citation using the example below from https://doi.org/10.1186/s13578-023-01012-8 
# Raw fastq files were trimmed to remove adaptors and low quality bases using trim_galore [84] ver- sion 0.6.6., a wrapper for cutadapt [85] (version 2.8)

# 84. Krueger F. Trim Galore. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. Accessed 1 Apr 2021.
# 85. Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet journal. 2011;17:10–2.