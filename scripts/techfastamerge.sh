# Modified from a script written by Stephen Justinen
echo techfastamerge.sh start

# Merge technical replicates for each sample

parallel 'zcat /data3/sjustinen/Hessinger_Data/technical_fastqs{}*.fastq.gz | gzip > zenodo/{}.merged.SE.fastq.gz' ::: Hal_1_rt Hal_2_16 Hal_3_16 Hal_3_rt Hal_4_16 Hal_5_16 Hal_5_rt Hal_6_16 Hal_6_rt Hal_7_rt

# End of script!
echo techfastamerge.sh end
