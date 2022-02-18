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

mamba install -c bioconda -c conda-forge -c anaconda trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 subread=2.0.1 multiqc=1.11 snakemake=6.15.0 parallel-fastq-dump=0.6.7 bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2 bioconductor-qvalue sambamba graphviz gffread tpmcalculator lxml

```
### Set-up

Configure SRA Toolkit.

```bash
# Set up the root dierctory for SRA files
# Enter '4' in interactive editor. Then enter your path: [absolute_path]/PGRP_Snakemake/output
vdb-config  -i --interactive-mode textual

# Set timeout to 100 seconds
vdb-config -s /http/timeout/read=100000
```
The pipeline requires an indexed reference genome and GTF file as input. By default, the pipeline expects the reference genome and gtf file to be located in the ref/ directory. You can either create a directory called 'ref' to store the reference genome and GTF file or update config.json to use a different location for the reference materials.

```bash
# Add the ref directory to you working directory
mkdir ref
cd ref

# Download a reference genome fasta file using curl, wget, ...etc.
# Download the reference GTF or GFF file. If downloading a GFF, convert it to GTF using the following (make sure gff is unzipped)
gffread [GFF_file] -T -F --keep-exon-attrs -o [output_name].gtf

# Update config.json with the path to the GTF file
vim config.json
# update "ref/[GTF_name]" under GTFName

# Index the reference genome with STAR (make sure genome fasta is unzipped)
# Helps to run on a compute cluster (computationally expensive)
STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles [genome.fa] --sjdbGTFfile [reference.gtf] --sjdbOverhang 99   --genomeSAindexNbases 12
```

## Running pipeline without scheduler
```bash
# Check the pipeline prior to run
snakemake -np

# Visualize the pipeline as a DAG
snakemake --dag [output] | dot -Tpdf -Gnodesep=0.75 -Granksep=0.75 > dag.pdf

# Run the pipeline 
snakemake --cores [available cores]
```

## Running pipeline with SLURM scheduler

```bash
sbatch --mem=16g -c 4 --time=13-11:00:00 -o snakemake.out -e snakemake.err --wrap="./run.sh"
```

## Citations

