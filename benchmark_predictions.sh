#!/bin/bash
# benchmark_predictions.sh
# This script generated benchmark plots for TARGET-SLs predictions on the CCLE cell lines

# change working dir to script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"


mkdir -p plots/benchmark

############################################################
# Read Config File                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"
source all_settings.cfg

############################################################
# Benchmark Driver Predictions                             #
############################################################

echo -e "\n\n---------------------------"
echo -e "Running driver gene benchmark"
echo -e "---------------------------\n\n"

#Rscript --vanilla "" -m $mode -n $network_choice -a $

############################################################
# Benchmark Essential Gene Predictions                     #
############################################################

echo -e "\n\n---------------------------"
echo -e "Running essential gene benchmark"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/benchmark_essential_gene_predictions.R" -n $network_choice -a $benchmark_algorithms -t $threads -c "all" -N $n_benchmark_predictions



