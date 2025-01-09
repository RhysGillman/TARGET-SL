#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Downloads all of the data required for comparing driver gene prioritisation algorithms using CCLE data"
   echo
   echo "Syntax: download_data [-h|a|c|g|d|s]"
   echo "options:"
   echo "h     Print help"
   echo "a     Download all data"
   echo "c     Download CCLE data"
   echo "g     Download GDSC data"
   echo "d     Download DGIdb data"
   echo "s     Download STRINGv11 data"
   echo "l     Download LOF/GOF Prediction data"
   echo "L     Download larger LOF/GOF Prediction data (Warning Large File)"
   echo
}
############################################################

############################################################
# Options                                                  #
############################################################

# Set Defaults
CCLE_download=0
GDSC_download=0
DGI_download=0
STRINGv11_download=0
lofgof_download=0
lofgof_large_download=0

# Process Input Options
while getopts ":hacgdsl" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        a) # download all data
            CCLE_download=1
            GDSC_download=1
            DGI_download=1
            STRINGv11_download=1
            lofgof_download=1
            lofgof_large_download=0;;
        c) # download CCLE Data
            CCLE_download=1;;
        g) # download GDSC Data
            GDSC_download=1;;
        d) # download DGI Data
            DGI_download=1;;
        s) # download STRINGv11 Data
            STRINGv11_download=1;;
        l) # download small LOFGOF Data
            lofgof_download=1;;
        L) # download large LOFGOF Data
            lofgof_large_download=1;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
   esac
done

############################################################
# Main Program                                             #
############################################################

CUR_DIR=$(pwd)

if (($CCLE_download==1))
then
    #----------------------------------------------------------------------
    # DepMap CCLE Data from https://depmap.org/portal/download/all/
    # Latest Update: 22Q4
    # Pipeline details available @ https://github.com/broadinstitute/depmap_omics
    # REF: Mahmoud Ghandi, Franklin W. Huang, Judit Jané-Valbuena, Gregory V. Kryukov, ... Todd R. Golub, Levi A. Garraway & William R. Sellers. 2019. Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature 569, 503-508 (2019)

    mkdir -p "$CUR_DIR"/data/CCLE
    cd "$CUR_DIR"/data/CCLE
    # 22Q4 Cell Line information
    curl -J -O -L "https://figshare.com/ndownloader/files/38466923"
    curl -J -O -L "https://figshare.com/ndownloader/files/38357447"
    # 22Q4 TPM CCLE Expression data (Log2(TPM+1), Protein-Coding Genes) (RSEM)
    curl -J -O -L "https://figshare.com/ndownloader/files/38357462"
    # 22Q4 Expected Counts (RSEM)
    curl -J -O -L "https://figshare.com/ndownloader/files/38357459"
    # 22Q4 Somatic Mutations (Mutec2)
    curl -J -O -L "https://figshare.com/ndownloader/files/38357492"
    # 22Q4 Gene-Level Copy Number
    curl -J -O -L "https://figshare.com/ndownloader/files/38357438"
    # 22Q4 Segment-Level Copy Number
    curl -J -O -L "https://figshare.com/ndownloader/files/38357441"
    # 22Q4 Gene Dependency
    curl -J -O -L "https://figshare.com/ndownloader/files/38357387"
    # 22Q4 Gene Effect
    curl -J -O -L "https://figshare.com/ndownloader/files/38357390"
    # 22Q4 README
    curl -J -O -L "https://figshare.com/ndownloader/files/38357510"
    # 23Q2 PRISM Drug Repurposing
    ## Matrix data
    curl -J -O -L "https://figshare.com/ndownloader/files/41408700"
    ## Matrix metadata
    curl -J -O -L "https://figshare.com/ndownloader/files/41408664"
    ## Long data
    curl -J -O -L "https://figshare.com/ndownloader/files/41408655"
    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi

if (($GDSC_download==1))
then
    #-------------------------------------------------------------------------
    # GDSC Drug Screening

    mkdir -p "$CUR_DIR"/data/GDSC
    cd "$CUR_DIR"/data/GDSC
    # 8.4 Cell Line Information
    #curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.4/Cell_Lines_Details.xlsx"
    # 8.4 Drug Information
    #curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.4/screened_compounds_rel_8.4.csv"
    # GDSC1 Release IC50s (24/7/22)
    #curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.4/GDSC1_fitted_dose_response_24Jul22.xlsx"
    # GDSC2 Release IC50s (24/7/22)
    #curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.4/GDSC2_fitted_dose_response_24Jul22.xlsx"
    
    curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx"
    curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx"
    curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/screened_compounds_rel_8.5.csv"
    curl -J -O -L "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/Cell_Lines_Details.xlsx"
    
    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi

if (($DGI_download==1))
then
    #-------------------------------------------------------------------------
    # DGI Drug Gene Interaction Database
    # Release Feb 2022
    mkdir -p "$CUR_DIR"/data/DGIdb
    cd "$CUR_DIR"/data/DGIdb
    #curl -J -O -L "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv"
    curl -J -O -L "https://www.dgidb.org/data/2022-Feb/interactions.tsv"
    curl -J -O -L "https://www.dgidb.org/data/2022-Feb/genes.tsv"
    curl -J -O -L "https://www.dgidb.org/data/2022-Feb/drugs.tsv"
    curl -J -O -L "https://www.dgidb.org/data/2022-Feb/interactions.tsv"
    
    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi



if (($STRINGv11_download==1))
then
    #-------------------------------------------------------------------------
    # STRINGv11
    mkdir -p "$CUR_DIR"/data/STRINGv11
    cd "$CUR_DIR"/data/STRINGv11/
    #v11 directionality
    curl -J -O -L "https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz"
    #v11 aliases
    curl -J -O -L "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz"
    # unzipping files
    gunzip 9606.protein.actions.v11.0.txt.gz
    gunzip 9606.protein.aliases.v11.0.txt.gz
    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi


if (($lofgof_download==1))
then
    #-------------------------------------------------------------------------
    # LOF GOF Mutation Predictions https://itanlab.shinyapps.io/goflof/ 
    mkdir -p "$CUR_DIR"/data/LOFGOF
    # https://www.sciencedirect.com/science/article/pii/S0002929721003840
    
    cd "$CUR_DIR"/data/LOFGOF/
    #HGMD-based Annotations
    curl -J -O -L "https://itanlab.shinyapps.io/goflof/_w_719a1fcc/session/e4e040aa8501ad86590a389bae9f6816/download/downloadFeats2?w=719a1fcc"
    #ClinVar-based annotations
    curl -J -O -L "https://itanlab.shinyapps.io/goflof/_w_719a1fcc/session/e4e040aa8501ad86590a389bae9f6816/download/downloadDataClinvar2?w=719a1fcc"

    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi


if ((lofgof_large_download==1))
then
    #-------------------------------------------------------------------------
    # LOF GOF Mutation Predictions https://itanlab.shinyapps.io/goflof/ 
    mkdir -p "$CUR_DIR"/data/LOFGOF
    # https://www.sciencedirect.com/science/article/pii/S0002929721003840
    
    cd "$CUR_DIR"/data/LOFGOF/
    
    # https://www.biorxiv.org/content/10.1101/2022.06.08.495288v2
    git clone "https://gitlab.com/itan-lab/logofunc-predictions.git"
    
    cd "logofunc-predictions"
    
    gunzip LoGoFuncVotingEnsemble_preds_final.csv.gz

    cd "$CUR_DIR"
    #-------------------------------------------------------------------------
fi
