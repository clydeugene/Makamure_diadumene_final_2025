#!/usr/bin/awk -f

# This script is used to convert a multi-line fasta file to a single line file
# Usage: awk -f linearizefasta.awk < input.fasta > output.fasta
# Source: https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345

/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}
     {printf("%s",$0);}
END  {printf("\n");}