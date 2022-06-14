# PGRP_Snakemake
## Introduction

The PGRP_Snakemake pipeline provides an efficient and modular workflow for processing RNA sequencing data from multiple sources. The workflow allows users to go from downloading raw data to differential expression analysis in a single run. Additionally, the pipeline can utilize raw sequencing reads directly from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) or from local storage. 

**Pipeline overview:**
- Download fastq from SRA (SRA Toolkit)
- Quality control on raw reads (FastQC)
- Trimming (Trim Galore)
- Quality control on trimmed reads (FastQC)
- Map reads to reference (STAR)
- Count reads (HTseq/FeatureCounts/TPMcalculator)
- Normalize counts TPM/FPKM (custom scripts)
- Summary statistics of normalized counts (custom scripts)
- Differential expression analysis (Trinity/DESeq2)

## Installation

Clone the repository:

```bash
git clone git@github.com:plantgenomicslab/PGRP_Snakemake.git
```

### Dependencies
- Trim Galore 0.6.7
- SRA Toolkit 2.11.0
- STAR
- HTseq 1.99.2
- Subread 2.0.1
- multiqc 1.11
- Snakemake 6.15.0
- parallel-fastq-dump 0.6.7
- Samtools 1.14
- ggplot2
- Trinity 2.13.2
- Graphviz
- Gffread
- TPMcalculator
- Bioconductor Qvalue

### Setting up a Conda environment 

We recommend setting up a conda environment for PGRP_Snakemake to easily install the above dependencies. The latest version of Miniconda can be installed [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links). Use the commands below to set up the PGRP_Snakemake environment.

```bash
conda create -n PGRP_Snakemake -c bioconda -c conda-forge python=3.7 mamba

conda activate PGRP_Snakemake

mamba install -c bioconda -c conda-forge -c anaconda trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 subread=2.0.1 multiqc=1.11 snakemake=7.5.0 parallel-fastq-dump=0.6.7 bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2 bioconductor-qvalue sambamba graphviz gffread tpmcalculator lxml

```
### Configuration

Add ```config.json``` to the repo. This file controls various inputs to the workflow and must be updated by the user. A template for ```config.json``` is availble in the ```example/``` directory. 
```
cd .../PGRP_Snakemake
cp example/example_config.json .
```

Configure SRA Toolkit (only necessary if using SRA). The following commands allow you to specify where large SRA files get stored and ensure that your connection doesn't time out when downloading data from NCBI's SRA database.

```bash
# Set up the root dierctory for SRA files
# Enter '4' in interactive editor. Then enter your path.
vdb-config  -i --interactive-mode textual

# Add the new path to config.json under 'sra_dir'
#vim config.json

# Set timeout to 100 seconds
vdb-config -s /http/timeout/read=100000
```

### Input Files
#### Reference genome
The pipeline requires an indexed reference genome and GTF file as input. By default, the pipeline expects the reference genome and gtf file to be located in the ```ref/``` directory. You can either create a directory called 'ref' to store the reference genome and GTF file or update config.json to use a different location for the reference materials.

To add a reference genome to the pipeline download fasta and GFF3 files from an appropriate source. Then:

```bash
# Add the ref directory to your working directory
mkdir ref && cd ref

# Make sure file is unzipped
gffread [GFF_file] -T -F --keep-exon-attrs -o [output_name].gtf

# Update config.json with the relative path to the GTF file
#vim config.json

# Index the reference genome with STAR (make sure genome fasta is unzipped)
# Helps to run on a compute cluster (computationally expensive)
cd ref
STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles [genome.fa] --sjdbGTFfile [reference.gtf] --sjdbOverhang 99   --genomeSAindexNbases 12
```

#### Workflow control file
The pipeline requires the user to create a control file called ```RunsbyExperiment.tsv```. An examples of this file is provided in ```examples/```.

```RunsbyExperiment.tsv``` provides the pipleline with information about the data that you want to process. These can be either SRA runs or locally stored data. When downloading SRA data, there may be multiple SRA runs (SRR...) for each SRA experiment (SRX...) where experiments represents the sequencing performed on a particular sample. Experiments should be given meainingful titles to aid in the interpretation of the pipeline output. Finally, the relationships between replicates and treatments should be flushed out.

```
# Example format of RunsbyExperiment.tsv
head sraRunsbyExperiment.tsv

Run	Experiment	Replicate	Treatment
SRR5210841	SRX2524297	ZT0_rep1	ZT0
SRR5210842	SRX2524297	ZT0_rep1	ZT0
SRR5210843	SRX2524297	ZT0_rep1	ZT0
SRR5210844	SRX2524298	ZT4_rep1	ZT4
SRR5210845	SRX2524298	ZT4_rep1	ZT4
SRR5210846	SRX2524298	ZT4_rep1	ZT4
SRR5210847	SRX2524299	ZT8_rep1	ZT8
SRR5210848	SRX2524299	ZT8_rep1	ZT8
SRR5210849	SRX2524299	ZT8_rep1	ZT8
```
Because these files can tedious to generate for projects with many samples, a script called joinSraRelations.py is provided. This script takes an SRA project id (SRP...) as input and fetches the associated run information over the SRA API. The final two inputs are regex expressions to parse out the human readable treatment/replicate text from the full SRA run titles.

```bash
# Example usage for SRA project SRP098160
# Full run titles for this project take the form: 'GSM2471308: ZT0_rep1; Glycine max; RNA-Seq'
# Experiment titles will be of the form 'ZT0_rep1'
./scripts/joinSraRelations.py SRP098160 "ZT\d{1,2}_rep\d" "ZT\d{1,2}"

# To not deal with regex and fix manually in text editor
./scripts/joinSraRelations.py SRP098160 ".*" ".*"
```

If **running locally**, the 'Run' and 'Replicate' fields of ```RunsbyExperiment.tsv``` my be omitted. The 'Replicate' field should represent the prefixes of all fastq files to be included in the analysis. 'Treatment' then should be the desired name of the treatment to which each replicate belongs. Make sure to update ```config.json``` with the nomenclature for paired-ends. 

#### DE control file
Another optional control file called ```contrasts.txt``` must be created if differential expression analysis is desired. If the file is missing then DE analysis will not be performed. ```contrasts.txt``` provides tab separated contrasts to be used in differential gene expression analysis using DESeq2. The first column defines treatments and the second column defines replicates associated with each treatment. Pairwise differentially expressed genes will be computed for each combination of treatments. See ```example/``` for an example:

```
# Example format of contrasts.txt
cat contrasts.txt

ZT0	ZT0_rep1
ZT0	ZT0_rep2
ZT0	ZT0_rep3
ZT4	ZT4_rep1
ZT4	ZT4_rep2
ZT4	ZT4_rep3
ZT8	ZT8_rep1
ZT8	ZT8_rep2
ZT8	ZT8_rep3
```

## Running the pipeline 

Separate snakefiles are provided for running from the SRA and from local data, namely ```Snakefile_SRA``` and ```Snakefile_local```. To run on SLURM there are also separate shell scripts for each option, ```run_SRA.sh``` and ```run_local```. Make sure to use the correct workflow for your analysis.

### Running without scheduler
```bash
# Check the pipeline prior to run
snakemake --snakefile [Snakemake_local|Snakemake_SRA] -np

# Visualize the pipeline as a DAG
snakemake --snakefile [Snakemake_local|Snakemake_SRA] \
          --dag [output] | \
          dot -Tpdf -Gnodesep=0.75 -Granksep=0.75 > dag.pdf

# Run the pipeline 
snakemake --snakefile [Snakemake_local|Snakemake_SRA] \
          --cores [available cores]
```

### Running with SLURM scheduler

```bash
# Running from sra
sbatch --mem=4g \
       -c 2 \
       --time=13-11:00:00 \
       -o snakemake.out \
       -e snakemake.err \
       --wrap="./run_SRA.sh"

# Running from local
sbatch --mem=4g \
       -c 2 \
       --time=13-11:00:00 \
       -o snakemake.out \
       -e snakemake.err \
       --wrap="./run_local.sh"
```

## Citations
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), pp. 10-12. doi:https://doi.org/10.14806/ej.17.1.200
- Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
- Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
- Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014
- Vera Alvarez R, Pongor LS, Mariño-Ramírez L, Landsman D. TPMCalculator: one-step software to quantify mRNA abundance of genomic features. Bioinformatics. 2019 Jun 1;35(11):1960-1962. doi: 10.1093/bioinformatics/bty896. PMID: 30379987; PMCID: PMC6546121.
- Anders, S., Pyl, P. T., & Huber, W. (2015). HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics (Oxford, England), 31(2), 166–169. https://doi.org/10.1093/bioinformatics/btu638
- Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
