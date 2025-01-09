#!/bin/bash
# create_benchmark_data.sh
# This wrapper runs the scripts needed to generate the benchmark data
# It is not necessary to run this script, as the benchmark data is available online

# change working dir to script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"

mkdir -p data
mkdir -p benchmark_data

############################################################
# Read Config File                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"
source all_settings.cfg

#download data




#choose cells

############################################################
# Choose Cells                                             #
############################################################

echo -e "\n\n---------------------------"
echo -e "Choosing CCLE Cell Lines"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/choose_cells.r"


############################################################
# Prepare Gold Standards                                   #
############################################################

echo -e "\n\n---------------------------"
echo -e "Preparing gold standards"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/prepare_all_gold_standards.r"


############################################################
# Formatting Benchmark Data                                #
############################################################

echo -e "\n\n---------------------------"
echo -e "Wrangling Benchmark Data"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/prepare_specific_data.r" -n $network_choice


############################################################
# Map Genomic LOF/GOFs                                     #
############################################################

echo -e "\n\n---------------------------"
echo -e "Mapping Genome Level LOF/GOFs"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/prepare_LOFGOF_annotations.R" -t $threads



############################################################
# Annotate CCLE LOF/GOFs                                   #
############################################################

echo -e "\n\n---------------------------"
echo -e "Annotating LOF/GOFs"
echo -e "---------------------------\n\n"

Rscript --vanilla "scripts/annotate_LOF_GOFs.R" -t $threads




