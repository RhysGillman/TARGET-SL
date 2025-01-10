# TARGET-SL

**TARGET-SL** (**T**umour-specific **A**lgorithm for **R**anking **GE**netic **T**argets via **S**ynthetic **L**ethality) is a framework for essential gene and drug sensitivity prediction based on personalised driver prioritisation. As such, it combines a wrapper function to simplify the execution of personlised driver prioritisation algorithms using the same input data with novel functions to predict essential genes from their output based on variant effect prediction and synthetic lethality.

This repository contains code from several personalised driver prioritisation algorithms that have been slightly modified to simplify their execution on our benchmarking dataset. We thank the original authors for the use of their code in our work. The original sources of this code are listed below:

### De Novo Network Methods (CSN_NCUA From Manuscript)

https://github.com/WilfongGuo/Benchmark_control

### DawnRank

https://github.com/MartinFXP/DawnRank

### OncoImpact

https://github.com/CSB5/OncoIMPACT

### PersonaDrive

https://github.com/abu-compbio/PersonaDrive

### PhenoDriverR

https://github.com/NWPU-903PR/PhenoDriverR

### PNC

https://github.com/NWPU-903PR/PNC

### PRODIGY

https://github.com/Shamir-Lab/PRODIGY

### SCS

Kindly provided by the original authors

### sysSVM2

https://github.com/ciccalab/sysSVM2


## Dependencies 

Some functions in TARGET-SL require specific dependencies. 

### R

We have supplied a file called install_R_dependencies.R in the scripts/ directory to simplify installing all of the R dependencies

```
Rscript --vanilla scripts/install_R_dependencies.R
```

### Python

PersonaDrive depends on several python packages that can be installed using the file in the scripts/PersonaDrive/requirements.txt file.

```
pip scripts/PersonaDrive/requirements.txt 
```

### MATLAB

MATLAB is required to run SCS, PNC and the de novo Network methods

### Gurobi

If users wish to run PNC or the de novo Network methods, they will need to setup gurobi (https://www.gurobi.com/) with a license

### ANNOVAR

An ANNOVAR (https://annovar.openbioinformatics.org/en/latest/) installation is required to run sysSVM2

## System Requirements

This software was tested in an Ubuntu linux system on a consumer-grade PC. Reproduction of manuscript figures should be possible on any system running R.

## Additional Data Files

Before running this code, you will need to acquire the benchmark data and pre-run results. These are available at https://doi.org/10.5281/zenodo.14625324. Once downloaded, unzip the contents and place the results/ and benchmark_data/ directories within the main TARGET-SL directory.

## How to Run TARGET-SL

This repository contains several wrapper functions to simplify running TARGET-SL. All options are controlled are controlled by making changes 
to the all_settings.cfg file.

TARGET-SL can be run in two different modes:
- **Benchmark Mode**: Benchmarks the performance of personalised driver prioritisation algorithms using cell line data from CCLE
- **Predict Mode**: Runs the user-specified driver prioritisation algorithms on user input data and predicts essential genes and drug sensitivity
  - Must be human hg38 data using HGNC gene symbols
  - Currently only DawnRank, OncoImpact, sysSVM2, and PersonaDrive are supported

### Benchmark Mode

1. Run personalised driver prioritisation algorithms (optional)

This is the longest running step of TARGET-SL. It can be skipped if the results directory is already populated from https://doi.org/10.5281/zenodo.14625324

```
bash prioritise_drivers.sh
```

2. Compare driver algorithms

This step generates plots to make initial comparisons of the output of the algorithms (upset plots) and their resource utilisation (memory and time)

```
bash compare_algorithms.sh
```

3. Predict essential genes and drug sensitivity

```
bash TARGET_SL.sh
```

4. Benchmark predictions

This step uses the output of the previous step and compares against ground truth gene essentiality and drug sensitivity data from the CCLE to benchmark performance

```
bash benchmark_predictions.sh
```

### Predict Mode

TARGET-SL predict mode is currently in a very early iteration as I have not had access to novel datasets to run. Expect this section to be updated soon. If you have any problems, please submit an issue and I will do my best to update the software.

Required data formattin is described in the relevant sections of the all_settings.cfg file. Once data is prepared as described, provide the file paths and change the mode to "predict"

```
bash prioritise_drivers.sh
bash TARGET_SL.sh
```
