#Stephen Justinen
echo techfastamerge.sh start
process ()
{
	zcat ${file}*.fastq.gz > ${file}.merged.SE.fastq
        gzip ${file}.merged.SE.fastq
}
for file in Hal_1_rt Hal_2_16 Hal_3_16 Hal_3_rt Hal_4_16 Hal_5_16 Hal_5_rt Hal_6_16 Hal_6_rt Hal_7_rt
do
	process &
done
wait
echo techfastamerge.sh end
