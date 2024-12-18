#!/bin/bash
# benchmark_predictions.sh
# This script generated benchmark plots for TARGET-SLs predictions on the CCLE cell lines

# change working dir to script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"


mkdir -p plots/benchmark
mkdir -p cache

############################################################
# Read Config File                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"
source all_settings.cfg


###########################################################
# Compare Resource Utilisation                            #
###########################################################

echo -e "\n\n---------------------------"
echo -e "Comparing Resource Utilisation"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/benchmark_resource_utilisation.R" -n $network_choice -a $benchmark_algorithms -c "all"

###########################################################
# Evaluation Prediction Similarity                        #
###########################################################

echo -e "\n\n---------------------------"
echo -e "Evaluating Prediction Similarity"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/"