# PGRP_Snakemake
PGRP_Snakemake

## Introduction

## Overview of pipeline





## Installation
***PLACE HOLDER***

### Dependencies
***PLACE HOLDER***

### Conda based

```bash
conda create -n PGRP_Snakemake -c bioconda -c conda-forge python=3.7 mamba

conda activate PGRP_Snakemake

mamba install -c bioconda -c conda-forge trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 subread=2.0.1 multiqc=1.11 snakemake=6.15.0 parallel-fastq-dump=0.6.7 bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2  bioconductor-qvalue

```
### Set-up
```bash
# Set up root dierctory for SRA files
# Enter '4' in interactive editor. Then enter your path: [absolute_path]/PGRP_Snakemake/output
vdb-config  -i --interactive-mode textual

# Set timeout to 100 seconds
vdb-config -s /http/timeout/read=100000
```
## Running pipeline without scheduler
***PLACE HOLDER***

## Running pipeline with SLURM scheduler
***PLACE HOLDER***
```bash
sbatch --mem=16g -c 4 --time=13-11:00:00 -o snakemake.out -e snakemake.err --wrap="./run.sh"
```

## Citations

####tximport

Soneson C, Love MI, Robinson MD (2015). “Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.” F1000Research, 4. doi: 10.12688/f1000research.7563.1.
