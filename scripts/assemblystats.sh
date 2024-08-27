#!/usr/bin/env bash

# This script is used to run assembly statistics on the assembly files

# Define input files
initialassembly=$1
filteredassembly=$2
nrblast=$3

# Define output files
initialassemblystats=$4
filteredassemblystats=$5
nrblaststats=$6

# Run assembly statistics
echo -e "\nRunning assembly statistics on initial assembly"

$TRINITY_HOME/util/TrinityStats.pl $initialassembly > $initialassemblystats

echo -e "\nRunning assembly statistics on filtered assembly"

$TRINITY_HOME/util/TrinityStats.pl $filteredassembly > $filteredassemblystats

echo -e "\nRunning assembly statistics on nr blast assembly"

$TRINITY_HOME/util/TrinityStats.pl $nrblast > $nrblaststats

# End of script
echo -e "\nAssembly statistics completed\n"
