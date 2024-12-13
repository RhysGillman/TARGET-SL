#!/bin/bash
# TARGET_SL.sh
# This script predicts essential genes and drug sensitivity based on driver gene predictions

# change working dir to script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"


############################################################
# Read Config File                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"
source all_settings.cfg

if [[ "$mode" == "benchmark" ]]
then
    echo -e "\n\n---------------------------"
    echo -e "Generating Random Predictions for Benchmarking"
    echo -e "---------------------------\n\n"
    Rscript --vanilla "scripts/random_predictions.R" -n $network_choice -c "ALL" -t $threads -c "ALL"
fi

############################################################
# Find Consensus Drivers                                   #
############################################################

echo -e "\n\n---------------------------"
echo -e "Finding Consensus Drivers"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/get_consensus_drivers.R" -m $mode -n $network_choice -a $consensus_algorithms

############################################################
# Essential Gene Predictions                               #
############################################################

echo -e "\n\n---------------------------"
echo -e "Predicting Essential Genes"
echo -e "---------------------------\n\n"


Rscript --vanilla "scripts/predict_essential_genes.R" -m $mode -s $sample_info -n $network_choice -a "ALL" -t $threads -c "ALL"

############################################################
# Drug Predictions                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Predicting Drug Sensitivity"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/predict_drug_sensitivity.R" -m $mode -s $sample_info -n $network_choice -a "ALL" -t $threads -c "ALL"

