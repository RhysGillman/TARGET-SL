#!/bin/bash
# prioritise_drivers.sh
# This script prioritises driver genes in selected samples using the selected algorithms

# change working dir to script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"


############################################################
# Functions                                                #
############################################################

memory_usage () {
    local pid=$1
    local maxmem=0
    while [ -d "/proc/${pid}" ]; do
        local mem=`cat /proc/${pid}/status | grep VmRSS | awk '{print $2}'`
        if [[ ${mem} -gt ${maxmem} ]]; then
            local maxmem=${mem}
        fi
        sleep 1
    done
    echo $maxmem
}

############################################################
# Read Config File                                         #
############################################################

echo -e "\n\n---------------------------"
echo -e "Reading Config File"
echo -e "---------------------------\n\n"
source all_settings.cfg

############################################################
# Setup Directories                                        #
############################################################

mkdir -p log
mkdir -p results/$mode/network_$network_choice
mkdir -p tmp

if [[ "$cancer_types" == "ALL" ]]
then
 
    dos2unix benchmark_data/network_${network_choice}/lineages.txt
    readarray -t cancer_types < benchmark_data/network_${network_choice}/lineages.txt

fi


for cancer_type in ${cancer_types[@]}; do
    echo -e "\n\n---------------------------"
    echo -e "Commencing Driver Prioritisation for cell type: $cancer_type"
    echo -e "---------------------------\n\n"
    

    ############################################################
    # Run DawnRank                                             #
    ############################################################
    
    if (($run_DawnRank==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running DawnRank for $cancer_type"
        echo -e "---------------------------\n\n"
        # Start time
        start=$(date +%s.%N)
        Rscript --vanilla "scripts/run_DawnRank.R" -m $mode -n $network_choice -c $cancer_type -s $sample_info -r $rna -d $mut > log/DawnRank_${network_choice}_${cancer_type}.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )
        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        #runtime=$(echo "$end - $start" | bc -l)
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running DawnRank"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/DawnRank_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/DawnRank_${network_choice}_${cancer_type}_stats.txt
        
        
    fi
    
    
    ############################################################
    # Run OncoImpact                                           #
    ############################################################
    
    if (($run_OncoImpact==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running OncoImpact for $cancer_type"
        echo -e "---------------------------\n\n"
        
        # Prepare data files for OncoImpact and store in tmp/
        Rscript --vanilla "scripts/prepare_OncoImpact_data.R" -m $mode -n $network_choice -c $cancer_type -t $threads -s $sample_info -r $rna -d $mut -v $cnv -w "$SCRIPT_DIR" -T $tmp
        # Start time
        start=$(date +%s.%N)
        ####### Running OncoImpact ########
        perl scripts/OncoImpact/oncoIMPACT.pl tmp/tmp_${cancer_type}_OncoImpact_config.cfg &> log/OncoImpact_${network_choice}_${cancer_type}.log & 

        # Get the process ID (PID) of the  script
        pid=$!
        # Monitor memory usage while running and wait until the process is done
        max_mem=$( memory_usage $pid )
        wait $pid
        # Get finish time
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running OncoImpact"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/OncoImpact_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/OncoImpact_${network_choice}_${cancer_type}_stats.txt
  
        # Formatting results
        Rscript --vanilla "scripts/format_OncoImpact_results.R" -n $network_choice -c $cancer_type -T $tmp
        # Remove temporary files
        rm -rf "results/network_$network_choice/OncoImpact/$cancer_type/ANALYSIS"
        rm -rf "results/network_$network_choice/OncoImpact/$cancer_type/COMPLETE_SAMPLES"
        rm -rf "results/network_$network_choice/OncoImpact/$cancer_type/INCOMPLETE_SAMPLES"
        rm tmp/tmp_${cancer_type}_OncoImpact_cnv.txt
        rm tmp/tmp_${cancer_type}_OncoImpact_config.cfg
        rm tmp/tmp_${cancer_type}_OncoImpact_EXP.txt
        rm tmp/tmp_${cancer_type}_OncoImpact_network.txt
        rm tmp/tmp_${cancer_type}_OncoImpact_SNP.txt
    fi
    
    ############################################################
    # Run PersonaDrive                                         #
    ############################################################
    
    if (($run_PersonaDrive==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PersonaDrive for $cancer_type"
        echo -e "---------------------------\n\n"
        
        Rscript --vanilla "scripts/prepare_PersonaDrive_data.R" -m $mode -n $network_choice -c $cancer_type -s $sample_info -r $rna -d $mut
        echo "1. Personalized Bipartite Networks (PBNs)..." > log/PersonaDrive_${network_choice}_${cancer_type}.log

        # Start time
        start=$(date +%s.%N)

        python3 scripts/PersonaDrive/constructing_PBNs.py -o "$SCRIPT_DIR/results/network_$network_choice/PersonaDrive/$cancer_type" >> log/PersonaDrive_${network_choice}_${cancer_type}.log &
    
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem1=$( memory_usage $pid )
        wait $pid

        echo "2 - Rank Mutated Genes ..." >> log/PersonaDrive_${network_choice}_${cancer_type}.log
        python3 scripts/PersonaDrive/PersonaDrive.py -o "$SCRIPT_DIR/results/network_$network_choice/PersonaDrive/$cancer_type" >> log/PersonaDrive_${network_choice}_${cancer_type}.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem2=$( memory_usage $pid )
        wait $pid

        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Take greatest maxmem value
        if [ $max_mem1 -gt $max_mem2 ]
        then
            max_mem=$max_mem1
        else
            max_mem=$max_mem2
        fi
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PersonaDrive"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/PersonaDrive_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PersonaDrive_${network_choice}_${cancer_type}_stats.txt

        Rscript --vanilla "scripts/format_PersonaDrive_results.R" -m $mode -n $network_choice -c $cancer_type -s $sample_info
        
        # Remove temporary files
        rm tmp/tmp_PersonaDrive_DEGs.csv
        rm tmp/tmp_PersonaDrive_MUT.csv
        rm tmp/tmp_PersonaDrive_network.tsv
    fi
    
    
    ############################################################
    # Run sysSVM2                                              #
    ############################################################
    
    if (($run_sysSVM2==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running sysSVM2 for $cancer_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p validation_data/network_$network_choice/ANNOVAR_input/$cancer_type
        mkdir -p results/network_$network_choice/sysSVM2/$cancer_type
        
        
        # Prepare input data
        
        Rscript --vanilla "scripts/create_ANNOVAR_input_files.R"  -m $mode -n $network_choice -c $cancer_type -d $maf -s $sample_info 
        
        # Start time
        start=$(date +%s.%N)
        
        # Run
        
        Rscript --vanilla "scripts/run_sysSVM.r" -m $mode -s $sample_info -n $network_choice -c $cancer_type -a $annovar_path -d $maf -v $cnv -V $cnv_CNR  > $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cancer_type}.log &
        
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running sysSVM2"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cancer_type}_stats.txt
        
        Rscript --vanilla "scripts/format_sysSVM2_results.R" -m $mode -n $network_choice -c $cancer_type
        
    
    
    fi
    
    ############################################################
    # Run PRODIGY                                              #
    ############################################################
    
    if (($run_PRODIGY==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PRODIGY for $cancer_type"
        echo -e "---------------------------\n\n"
        
        # Start time
        start=$(date +%s.%N)
        Rscript --vanilla "scripts/run_PRODIGY.R" -n $network_choice -c $cancer_type -t $threads > log/PRODIGY_${network_choice}_${cancer_type}.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )
        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PRODIGY"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/PRODIGY_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PRODIGY_${network_choice}_${cancer_type}_stats.txt
    fi
    
    ############################################################
    # Run SCS                                                  #
    ############################################################
    
    if (($run_SCS==1))
    then
        
        echo -e "\n\n---------------------------"
        echo -e "Running SCS for $cancer_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/$mode/network_$network_choice/SCS/$cancer_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_SCS_data.R" -m $mode -n $network_choice -c $cancer_type
        
        if [[ "$network_choice" == "own" ]]
        then
          cp -f data/own_networks/SCS.mat tmp/tmp_network.mat
          cd scripts
          cd SCS
        else
          cd scripts
          matlab -batch "create_matlab_network('../benchmark_data/network_$network_choice/network_directed.csv')"
          cd SCS
        fi
        

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_SCS('$network_choice', '$cancer_type', '$mode')" > $SCRIPT_DIR/log/SCS_${network_choice}_${cancer_type}.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running SCS"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/SCS_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/SCS_${network_choice}_${cancer_type}_stats.txt

        cd ..
        matlab -batch "get_SCS_results_names('$network_choice', '$cancer_type', '$mode')"
        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_SCS_results.R" -m $mode -n $network_choice -c $cancer_type
        
        # Remove temporary files
        rm scripts/SCS/CNV_internate.txt
        rm scripts/SCS/EXPR_internate.txt
        rm scripts/SCS/SNP_internate.txt
        rm tmp/tmp_SCS_cnv.txt
        rm tmp/tmp_SCS_EXP.txt
        rm tmp/tmp_SCS_SNP.txt
        rm tmp/tmp_network.mat
        
    fi
    
    ############################################################
    # Run PNC                                                  #
    ############################################################
    
    if (($run_PNC==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PNC for $cancer_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/$mode/network_$network_choice/PNC/$cancer_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_PNC_data.R" -m $mode -n $network_choice -c $cancer_type


        if [[ "$network_choice" == "own" ]]
        then
          cp -f data/own_networks/PNC.mat tmp/tmp_network.mat
          cd scripts
          cd PNC
        else
          cd scripts
          matlab -batch "create_matlab_network('../benchmark_data/network_$network_choice/network_directed.csv')"
          cd PNC
        fi

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_PNC('$network_choice', '$cancer_type', '$gurobi_path','$mode')" > $SCRIPT_DIR/log/PNC_${network_choice}_${cancer_type}.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PNC"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/PNC_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/PNC_${network_choice}_${cancer_type}_stats.txt

        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_PNC_results.R" -m $mode -n $network_choice -c $cancer_type
        
        
        # Remove temporary files
        rm tmp/tmp_PNC_pseudonormal_expression.txt
        rm tmp/tmp_PNC_tumour_expression.txt
        rm tmp/tmp_network.mat
    fi
    
    ############################################################
    # Run De Novo Methods                                      #
    ############################################################
    
    if (($run_combined_de_novo==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running Combined De Novo Methods for $cancer_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/$mode/network_$network_choice/combined_de_novo_methods/$cancer_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_combined_de_novo_methods_data.R" -n $network_choice -c $cancer_type


        if (($network_choice=="own"))
        then
          cp -f data/own_networks/de_novo_net2.mat tmp/tmp_network.mat
          cd scripts
          cd combined_de_novo_methods
        else
          cd scripts
          matlab -batch "create_matlab_network('../benchmark_data/network_$network_choice/network_directed.csv')"
          cd combined_de_novo_methods
        fi

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_Benchmark_control('$network_choice', '$cancer_type', '$gurobi_path')" > $SCRIPT_DIR/log/combined_de_novo_methods_${network_choice}_${cancer_type}.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running combined de novo methods"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/combined_de_novo_methods_${network_choice}_${cancer_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/combined_de_novo_methods_${network_choice}_${cancer_type}_stats.txt

        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_combined_de_novo_methods_results.R" -n $network_choice -c $cancer_type
        
        # Remove temporary files
        rm tmp/tmp_combined_de_novo_methods_pseudonormal_expression.txt
        rm tmp/tmp_combined_de_novo_methods_tumour_expression.txt
        rm tmp/tmp_network.mat
        
    fi
    
    ############################################################
    # Run PhenoDriverR                                         #
    ############################################################
    
    if (($run_PhenoDriverR==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PhenoDriverR for $cancer_type"
        echo -e "---------------------------\n\n"
        
        
        # Start time
        start=$(date +%s.%N)
        
        # Run
        
        Rscript --vanilla "scripts/run_PhenoDriverR.R" -n $network_choice -c $cancer_type -t $threads  > "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cancer_type}.log" &
        
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PhenoDriverR"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cancer_type}_stats.txt"
        echo -e "$runtime\t$max_mem" >> "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cancer_type}_stats.txt"
        
    
    fi
    
    ############################################################
    # Run badDriver Simulation                                 #
    ############################################################
    
    #if (($run_badDriver==1))
    #then
    #    echo -e "\n\n---------------------------"
    #    echo -e "Running badDriver Simulation for $cancer_type"
    #    echo -e "---------------------------\n\n"
        
        
        # Run
        
    #    Rscript --vanilla "scripts/badDriver.R" -n $network_choice -c $cancer_type -t $threads -r $badDriver_ref_set -s $badDriver_n_sim
        
    
    
    #fi
    
    
done
