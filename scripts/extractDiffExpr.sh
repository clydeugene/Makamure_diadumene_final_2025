#!/bin/bash

# Start of script
      echo -e "\nThe script extractDiffExpr.sh has started in $(pwd) at $(date)\n"

# Set variables
      matrix="$1"
      samples="$2"
      log2fold="$3"
      q_value="$4"
      diffExprlogfile="$5"
      clusterslog="$6"
      analysis_type="$7"
      # log2fold=$(printf "%.1f" "$(echo "$2/2" | bc -l)")
      # q_value=$(echo $(printf "%.0e" "$3" | sed 's/-0/-/'))      

# Extract those differentially expressed (DE) genes that are at least 2-fold differentially expressed at a q-value of <= 0.05 in any of the pairwise sample comparisons
      $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
      --matrix ../../../$matrix \
      --samples ../../../$samples \
      -P $log2fold \
      -C $q_value \
       > ../../../$logfile 2>&1
            # The above generates several output files with a prefix diffExpr.P5e-2_C1', indicating the parameters chosen for filtering, where P (FDR actually) is set to 0.05, and fold change (C) is set to 2^(1) or 2-fold.

# Extract the differentially expressed genes or isoforms into a separate file
      awk -F "\t" -v OFS="\t" 'NR>1 {print $1}' diffExpr.P"$q_value"_C"$log2fold".matrix | sort | uniq | awk -v header="$analysis_type" 'BEGIN{print header} {print}' > differentially_expressed_"$analysis_type".txt


# Report the number of differentially expressed genes or isoforms
      echo -e "\nThe number of differentially expressed $analysis_type is $(wc -l differentially_expressed_"$analysis_type".txt | awk '{print $1}')"

# Extract transcript clusters by expression profile by cutting the dendrogram
      # $TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
      # --Ptree 60 \
      # -R diffExpr.P"$q_value"_C"$log2fold".matrix.RData \
      # > ../../../$clusterslog 2>&1

# End of script
      echo -e "\nThe script extractDiffExpr.sh in $(pwd) has ended at $(date)\n"
