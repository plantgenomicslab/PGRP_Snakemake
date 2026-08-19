# PGRP_Snakemake
## Introduction

The PGRP_Snakemake pipeline provides an efficient and modular workflow for processing RNA sequencing data from multiple sources. The workflow provides users with a complete and distributed pipeline from downloading raw data to differential expression analysis. The pipeline can utilize raw sequencing reads directly from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) or from local storage and is compatable with paired- or single-end libraries.

**Pipeline overview:**
- Download fastq from SRA (SRA Toolkit)
- Quality control on raw reads (FastQC)
- Trimming (fastp)
- Quality control on trimmed reads (FastQC)
- Map reads to reference (STAR)
- Count reads (RSEM/HTseq/FeatureCounts/TPMcalculator)
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
- snakemake 7.24.0
- parallel-fastq-dump 0.6.7
- Samtools 1.14
- ggplot2
- Trinity 2.13.2
- Graphviz
- Gffread
- TPMcalculator
- Bioconductor Qvalue
- RSEM
- tabulate 0.8.10
- fastp 0.23.2

### Setting up a Conda environment 

We recommend [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for fast, single-binary environment management. If you don't have it, conda (or mamba) works as a drop-in replacement — just substitute the binary name in the commands below.

#### Quick install (micromamba)

```bash
# Linux/macOS one-liner; installs to ~/.local/bin/micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

If `micromamba` is unavailable on your system, replace `micromamba` with `conda` (or `mamba`) in every command in this README — flags and arguments are compatible.

#### Set up the environment

```bash
# Create env + install all dependencies in one shot
micromamba create -n PGRP_Snakemake -c bioconda -c conda-forge -c anaconda \
  python=3.7 \
  tabulate=0.8.10 trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 \
  subread=2.0.1 multiqc=1.11 snakemake=7.24.0 parallel-fastq-dump=0.6.7 \
  bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2 \
  bioconductor-qvalue sambamba graphviz gffread tpmcalculator lxml rsem fastp=0.23.2

micromamba activate PGRP_Snakemake
```

> Conda fallback: replace `micromamba create` with `conda create` and `micromamba activate` with `conda activate`. The `-c` channel order is identical.
### Configuration

Add ```config.json``` to the repo. This file controls various inputs to the workflow and must be updated by the user. A template for ```config.json``` is availble in the ```example/``` directory. 
```
cd ../PGRP_Snakemake
cp example/example_config.json config.json
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

## Input Files
### Reference genome
The pipeline requires an indexed reference genome and GTF file as input. To add a reference genome to the pipeline download fasta and GFF3 files from an appropriate source. Then:

```bash
# Make sure file is unzipped
gffread [GFF_file] -T -F --keep-exon-attrs -o [genome].gtf

# Update config.json with the relative path to the GTF file and the reference folder
#vim config.json

# Index the reference genome with STAR (make sure genome fasta is unzipped)
# Helps to run on a compute cluster (computationally expensive)
STAR  --runThreadN 48 \
  --runMode genomeGenerate \
  --genomeDir . \
  --genomeFastaFiles [genome.fa] \
  --sjdbGTFfile [genome.gtf] \
  --sjdbOverhang 99 \
  --genomeSAindexNbases 12

# If using RSEM, prepare the reference 
rsem-prepare-reference -p 48 --gtf [genome.gtf] [genome.fa] [rsem_prep]
```

### Workflow control file
The pipeline requires the user to create a control file called ```RunsByExperiment.tsv```. An examples of this file is provided in ```examples/```.

```RunsByExperiment.tsv``` provides the pipleline with information about the data that you want to process. These can be either SRA runs or locally stored data. When downloading SRA data, there may be multiple SRA runs (SRR...) for each SRA experiment (SRX...) where experiments represents the sequencing performed on a particular sample. Experiments should be given meainingful titles to aid in the interpretation of the pipeline output. Finally, the relationships between replicates and treatments should be flushed out for count aggregation and DEG analysis.

#### Columns

| Column | Required | Meaning |
|--------|----------|---------|
| `Run` | yes | One sequencing run. For local data, the fastq file prefix. |
| `Experiment` | no | SRA experiment accession (SRX...). May be omitted for local data. |
| `Replicate` | yes | Biological replicate a run belongs to. Several runs may share one replicate. |
| `Treatment` | yes | Condition/group a replicate belongs to. Contrasts are computed between treatments. |

The treatment column is named `Treatment`. Control files written before this was
standardised used `Sample` for the same field; that name is still accepted, but
`Treatment` is what the examples and the generator scripts produce.

#### Generating the control file

Two generators are provided for local data — pick the one that matches your file
naming. Both write to the current directory.

```bash
# Paired-end reads named '<treatment>_rep<N>_R1.fastq.gz' / '..._R2.fastq.gz'.
# Also writes a pairwise sample_contrasts file covering every treatment pair.
./scripts/create_RunsbyExperiment.py /path/to/reads
#   -> RunsByExperiment_<timestamp>.tsv   (Run, Treatment, Replicate, Sample)
#   -> sample_contrasts_<timestamp>.tsv   (one treatment pair per line)
# Rename/point config.yml at these before running the pipeline.

# Any .fq/.fastq/.fq.gz/.fastq.gz naming. Strips a trailing _rep<N>/_Rep<N>
# from each file prefix to derive the treatment name.
python scripts/Experiment_name_composer.py /path/to/reads [output.tsv]
#   -> RunsByExperiment.tsv               (Run, Replicate, Treatment)
```

`Experiment_name_composer.py` puts the original file prefix in both `Run` and
`Replicate`, and the prefix with the `_rep*`/`_Rep*` suffix removed in
`Treatment`.

```
# Example format of RunsbyExperiment.tsv

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

If **running locally**, the 'Run' and 'Experiment' fields of ```RunsbyExperiment.tsv``` may be omitted. The 'Replicate' field should represent the prefixes of all fastq files to be included in the analysis. 'Treatment' then should be the desired name of the treatment to which each replicate belongs. Make sure to update ```config.json``` with the nomenclature for paired-ends. 

#### DE control files
Differential expression needs a second control file, ```replication_relationship.txt``` (its path is set by ```rep_relations``` in ```config.yml```). It is the ```--samples_file``` handed to DESeq2 via Trinity: tab separated, first column the treatment, second column a replicate belonging to it, **one line per replicate**. Pairwise differentially expressed genes are computed for each combination of treatments.

The pipeline regenerates this file from ```RunsByExperiment.tsv``` at the path given by ```rep_relations``` every time the Snakefile is parsed, so there is normally nothing to write by hand. To produce or inspect it up front — useful when setting up contrasts before committing to a full run — use:

```bash
./scripts/make_replication_relationship.py
#   reads  RunsByExperiment.tsv
#   writes replication_relationship.txt

# non-default paths
./scripts/make_replication_relationship.py -i my_runs.tsv -o deg_samples.txt
```

If you do want to maintain this file by hand — to analyse a subset of replicates, say — just edit it. Whenever its contents disagree with ```RunsByExperiment.tsv``` the pipeline treats it as yours, leaves it untouched, and prints a warning naming the differences rather than overwriting your edits. Delete the file (or pass ```--force``` to the script) to go back to a generated one.

Note that ```RunsByExperiment.tsv``` holds one row per *run* while this file holds one row per *replicate*: runs sharing a replicate are collapsed to a single line.

```
# Example format of replication_relationship.txt
cat replication_relationship.txt

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

## Optional: BBDuk pre-alignment contaminant filter

`PGRP_Snakemake` ships with an optional pre-alignment filter (BBDuk single-pass at `k=31`) that removes rRNA, chloroplast, mitochondrion, and vector reads before STAR alignment. On plant total-RNA libraries this typically lifts STAR's uniquely-mapped rate by 5–30 percentage points; the gain is smaller for polyA-selected libraries. The filter is **off by default** (`bbduk_enable: false` in `config.yml`); enable it once references are built.

### Reference set

Four references are concatenated for a single BBDuk pass:

- **rRNA k-mers** — BBTools-bundled `ribokmers.fa.gz` (broad-phylogeny rRNA)
- **Chloroplast (12 species)** — `data/bbduk_refs/cp_12sp.fa.gz`
- **Mitochondrion (12 species)** — `data/bbduk_refs/mt_12sp.fa.gz`
- **Vector** — NCBI UniVec → `data/bbduk_refs/univec.fa.gz`

The 12 species span all major land-plant clades (algae → bryophyte → lycophyte → gymnosperm → seven angiosperm clades). The canonical accession list lives at `data/bbduk_refs/SPECIES_LIST.tsv`; edit that file and re-run the build script to add or swap species.

### Setup

```bash
# 1. Install env once (micromamba preferred; replace with `conda` if needed)
micromamba env create -f envs/bbduk.yaml
micromamba activate pgrp-bbduk

# 2. Build references once (~5 min, fetches 24 GenBank records via efetch + UniVec via curl)
bash scripts/build_bbduk_refs.sh

# 3. Enable in config.yml
#    bbduk_enable: true
```

### Driver-side environment

`rules/bbduk_filter.smk` resolves `ribokmers.fa.gz` by inspecting the BBTools install relative to `$(which bbduk.sh)`. This resolution runs in the Snakemake **driver process**, not on the SLURM workers, so the `pgrp-bbduk` conda env (or at least `bbduk.sh` on PATH) must be active wherever you launch `snakemake`. If you can't activate the env on the driver — common when the driver lives on a login node with no `bbmap` install — set `bbduk_ribokmers` in `config.yml` to an absolute path to `ribokmers.fa.gz` to skip auto-resolution.

### Output

For each sample the rule writes:

- `output/{replicate}/{sample}/bbduk/{sample}{_1,_2}_clean.fq.gz` — fed to STAR
- `output/{replicate}/{sample}/bbduk/{sample}{_1,_2}_contam.fq.gz` — archived contaminant reads
- `output/{replicate}/{sample}/bbduk/{sample}_stats.txt` and `{sample}_refstats.txt` — per-sample contamination breakdown
- `output/{replicate}/{sample}/logs/{sample}_bbduk.log` — BBDuk stderr

Aggregate per-sample stats into one tidy TSV after a run:

```bash
python scripts/aggregate_bbduk_stats.py output bbduk_summary.tsv
```

> Note: for paired-end libraries the `total_reads` / `contam_reads` / `clean_reads` columns count individual reads (R1 + R2), not read pairs — that is the convention BBDuk's own `stats.txt` `#Total` line uses.

### Disabling

Set `bbduk_enable: false` in `config.yml` (the default) to bypass the filter; the DAG falls back to `trim_*` → `align_*` directly with no other changes required.

### Smoke test

```bash
export SMOKE_R1=/path/to/R1.fq.gz SMOKE_R2=/path/to/R2.fq.gz
bash tests/test_bbduk_filter.sh
```

Expected: `PASS` plus a refstats table with non-zero rRNA / cp / mt / univec read counts.

## Running the pipeline 

### Running without scheduler
```bash
# Check the pipeline prior to run
snakemake --snakefile Snakefile -np

# Visualize the pipeline as a DAG
snakemake --snakefile Snakefile \
          --dag [output] | \
          dot -Tpdf -Gnodesep=0.75 -Granksep=0.75 > dag.pdf

# Run the pipeline 
snakemake --snakefile Snakefile \
          --cores [available cores]
```

### Running with SLURM scheduler

```bash
sbatch --mem=4g \
       -c 2 \
       --time=13-11:00:00 \
       -o snakemake.out \
       -e snakemake.err \
       --wrap="./run.sh"
```
## Upgrading an existing checkout

If you already have a working directory with results in it, pulling is safe for your data but needs one manual step.

**`config.yml` is tracked by git.** Your copy almost certainly has local edits (your `genomeDir`, `GTFname`, `rawInputDir`, `RSEM_prepared_genome`), and upstream changes touch the same region, so a plain `git pull` will refuse or conflict:

```bash
git status                 # confirm what is locally modified
git stash push -m pre-pull
git pull
git stash pop              # if config.yml conflicts, keep YOUR values
```

You do **not** need to adopt the upstream `config.yml`. The only functional change there is that `ref:` became optional — leaving your existing `ref:` line in place works exactly as before.

**Nothing else needs touching.** `RunsByExperiment.tsv` and `replication_relationship.txt` are untracked, and `output/` and `.snakemake/` are gitignored, so a pull cannot disturb finished results. Completed alignments are preserved and will not re-run.

Then resume normally:

```bash
bash run.sh
```

`run.sh` already passes `--rerun-incomplete`. If the previous run was killed and you get `Directory cannot be locked`, run `snakemake --unlock --snakefile Snakefile` once first.

Two things you may notice on the next run:

- **`replication_relationship.txt` gets shorter.** It used to list one line per sequencing *run*, so replicates spanning several runs appeared more than once. It is now one line per replicate. If you have edited the file by hand, the pipeline detects that its contents disagree with `RunsByExperiment.tsv`, leaves it alone, and prints a warning naming the differences — delete the file, or run `scripts/make_replication_relationship.py --force`, to go back to a generated one.
- **A stray `refcds_length.tsv` next to your reference directory.** HTseq normalization used to write `cds_length.tsv` by gluing the filename onto `genomeDir`, landing it one level up. It now goes inside `genomeDir`. The old file is orphaned and can be deleted.

## error checking

# Slurm accounting for both jobs
`sacct -j 5396911,5396912 --format=JobID,State,ExitCode,Elapsed,NodeList,MaxRSS,MaxVMSize%20`

# Full Slurm job records (helpful if they were OOM-killed or preempted)
```
scontrol show job -dd 5396911
scontrol show job -dd 5396912
```

# Tool logs from your rules
```
tail -n +200 output/AgteqT02_rep2/AgteqT02_rep2/logs/AgteqT02_rep2_raw_fastqc.log
tail -n +200 output/AgteqT12_rep1/AgteqT12_rep1/logs/AgteqT12_rep1_trim.log
```
## Nextflow
```
nextflow run main.nf -params-file params.yaml -c resources-slurm.config -resume -with-report report.html -with-trace trace.txt -with-timeline timeline.html -with-dag flowchart.png
```

## Citations
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), pp. 10-12. doi:https://doi.org/10.14806/ej.17.1.200
- Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
- Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
- Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014
- Vera Alvarez R, Pongor LS, Mariño-Ramírez L, Landsman D. TPMCalculator: one-step software to quantify mRNA abundance of genomic features. Bioinformatics. 2019 Jun 1;35(11):1960-1962. doi: 10.1093/bioinformatics/bty896. PMID: 30379987; PMCID: PMC6546121.
- Anders, S., Pyl, P. T., & Huber, W. (2015). HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics (Oxford, England), 31(2), 166–169. https://doi.org/10.1093/bioinformatics/btu638
- Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
