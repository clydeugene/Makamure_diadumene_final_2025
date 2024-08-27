#!/bin/bash

# Define input/output files

input_file=$1
cdhit_file=$2
output_file=$3
barplot_file=$4

# Print start time and notify user of the start of the data extraction process

echo -e "\n$(date): Starting to extract data from input file: $input_file\n"

# Extract data from input file

awk '{if($2!~"Database:"){print}}' $input_file | 
awk '{if ($3!~"hits"){print}}' | 
awk '/Fields/{print a}{a=$0} f{print;f=0} /Fields/{f=1}' | 
awk -v OFS="\t" '{if($2~"Query"){print $0, $3}else{print}}' | 
awk  -v OFS="\t" -F '\t' '{if($1~"Query"){a=$2}else{print a, $1, $2, $3, $4, $5, $6, $7, $8 , $9, $10, $11}}' | 
sort -k1b,1 | 
grep -v BLASTN > $output_file

# Print end time and notify user of the number of lines in the output file

echo -e "\n$(date): Data extraction from input file: $input_file is complete\n\nThe output file has $(wc -l $output_file | awk '{print $1}') lines\n"

# Make bar plot of the extracted results

A=$(grep -v Diadumene $output_file | grep anemone | wc -l)
B=$(grep Diadumene $output_file | wc -l)
C=$(grep -v Diadumene $output_file | grep coral | wc -l)
D=$(cat $output_file | wc -l)
E=$(grep TRINITY $cdhit_file | wc -l)

Rscript scripts/make_anemone_barplot.R $A $B $C $D $E $barplot_file