import sys
import os
import pandas as pd

SAMPLES_FILE = pd.read_csv("sraRunsbyExperiment.tsv", sep="\t")
SAMPLE_LIST = list(set(SAMPLES_FILE["Run"].values.tolist()))
REPLICATE_LIST = list(set(SAMPLES_FILE["Replicate"].values.tolist()))
REPLICATE_LOOKUP = SAMPLES_FILE.groupby("Replicate")['Run'].apply(list).to_dict()

PAIR_LIST=["1", "2"]

CONTRASTS_FILE = pd.read_csv("contrasts.txt", sep="\t", header=None)

# Create sample output folders
os.makedirs("output/counts/", exist_ok=True)
os.makedirs("output/sra/", exist_ok=True)
os.makedirs("output/sra/logs/", exist_ok=True)
os.makedirs("output/logs/", exist_ok=True)
os.makedirs("output/DEG/", exist_ok=True)
os.makedirs("output/counts/featureCounts/", exist_ok=True)
os.makedirs("output/counts/htseq/", exist_ok=True)
os.makedirs("output/counts/tpmcalculator", exist_ok=True)

for path in REPLICATE_LIST:
    os.makedirs("output/" + path + "/raw/", exist_ok=True)
    os.makedirs("output/" + path + "/bam/", exist_ok=True)
    os.makedirs("output/" + path + "/logs/", exist_ok=True)
    os.makedirs("output/" + path + "/trim/", exist_ok=True)

# Load the config file as a dictionary
cf = open("config.json")
config_dict = json.load(cf)
cf.close()

# Check for an indexed reference genome based on required files in config.json
for ref in config_dict["ref"]:
    if not os.path.exists(ref):
        print("The " + ref + " file is missing. Have you downloaded and indexed a refernce genome? Terminating...")
        exit()

configfile: "config.json"

def allInput():
    inputs = ["output/counts/tpmcalculator/tpm-calculator.tpm",
              "output/counts/featureCounts/featureCounts.cnt",
              "output/counts/htseq/htseq-count.tsv",
              "output/counts/featureCounts/featureCounts.tpm.tsv",
              "output/counts/featureCounts/featureCounts.fpkm.tsv",
              "output/counts/htseq/htseq-count.tpm.tsv",
              "output/counts/htseq/htseq-count.fpkm.tsv",
              "output/counts/tpmcalculator/tpm-calculator.tpm.mean.tsv",
              "output/counts/featureCounts/featureCount_clean.cnt",
              "output/DEG/featureCount_clean.cnt." + CONTRASTS_FILE.iloc[1][0]  + "_vs_" + CONTRASTS_FILE.iloc[-1][0]  + ".DESeq2.DE_results"]
    for replicate in REPLICATE_LIST:
        inputs.append("output/" + replicate + "/bam/" + replicate + ".bam")
        for sample in REPLICATE_LOOKUP[replicate]:
            inputs.append("output/" + replicate + "/" + sample + "/bam/" + sample + ".bamAligned.sortedByCoord.out.bam")
            for pair in PAIR_LIST:
                inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + "_" + pair + ".fastq.gz")
                inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + "_" + pair + "_fastqc.zip")
                inputs.append("output/" + replicate + "/" + sample + "/trim/" + sample + "_" + pair + ".fq.gz")
    return(inputs)

rule all:
    input:
        allInput()

rule fetchSRA:
    output: "output/sra/{sample}.sra"
    message: "-----Fetching {wildcards.sample} SRA files-----"
    threads: config["threads"]["fetchSRA"]
    log: "output/sra/logs/{sample}_downloadSRA.log"
    run:
        shell("prefetch {wildcards.sample} --output-file {output} 2> {log}")
        # Check sra files to make sure they are valid
        shell("echo '--------Validating {wildcards.sample}.sra--------'")
        shell("set +e")
        shell("if ! vdb-validate {output}; then echo 'vdb-validate found errors in {output}. Check log here:{log}. Exiting......' && exit 1; fi")

rule convertSRAtoFastq:
    input: "output/sra/{sample}.sra"
    output:
        "output/{replicate}/{sample}/raw/{sample}_1.fastq.gz",
        "output/{replicate}/{sample}/raw/{sample}_2.fastq.gz"
    message: "-----Converting {wildcards.sample} SRA to Fastq files-----"
    threads: config["threads"]["convertSRAtoFastq"]
    log: "output/{replicate}/{sample}/logs/{sample}_fastqdump.log"
    run:
        shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-e --gzip \
				--outdir output/{wildcards.replicate}/{wildcards.sample}/raw \
				2> {log}")

rule fastqc_raw:
    input:
        "output/{replicate}/{sample}/raw/{sample}_1.fastq.gz",
        "output/{replicate}/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{replicate}/{sample}/raw/{sample}_1_fastqc.zip",
        "output/{replicate}/{sample}/raw/{sample}_2_fastqc.zip"
    message: "-----Running Fastqc_raw {wildcards.sample}-----"
    threads: config["threads"]["fastqc_raw"]
    log: "output/{replicate}/{sample}/logs/{sample}_raw_fastqc.log"
    run:
        shell("fastqc --threads {threads} {input} 2> {log}")

rule trim:
    input:
        fwd_fastq = "output/{replicate}/{sample}/raw/{sample}_1.fastq.gz",
        rev_fastq = "output/{replicate}/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{replicate}/{sample}/trim/{sample}_1.fq.gz",
        "output/{replicate}/{sample}/trim/{sample}_2.fq.gz"
    message: "-----Trimming {wildcards.sample}-----"
    log: "output/{replicate}/{sample}/logs/{sample}_trim.log"
    threads: config["threads"]["trim"]
    run:
        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
             --cores {threads} --max_n 40 --gzip -o output/{wildcards.replicate}/{wildcards.sample}/trim {input.fwd_fastq} {input.rev_fastq} \
             2> {log}")
        shell("for file in output/{wildcards.replicate}/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc --threads {threads} {output}")

rule align:
    input:
        fwd_fastq = "output/{replicate}/{sample}/trim/{sample}_1.fq.gz",
        rev_fastq = "output/{replicate}/{sample}/trim/{sample}_2.fq.gz"
    output:"output/{replicate}/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
    message: "-----Aligning {wildcards.sample}-----"
    log: "output/{replicate}/{sample}/logs/{sample}_align.log"
    threads: config["threads"]["align"]
    run:
        shell("STAR --runMode alignReads --runThreadN {threads} \
                --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 25000 \
		--genomeDir " + config["genomeDir"]  + " \
	        --readFilesCommand gunzip -c --readFilesIn {input.fwd_fastq} {input.rev_fastq} \
                --outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.replicate}/{wildcards.sample}/bam/{wildcards.sample}.bam  \
		2> {log}")
        # Perform quick check on output bam file to ensure it is not corrupted
        shell("echo '--------Checking {output}----------'")
        shell("set +e")
        shell("if ! samtools quickcheck {output}; then echo 'samtools quickcheck found errors in {output}. Check log here: {log}. Exiting......' && exit 1; fi")

def mergeInput(wildcards):
    inputs = expand("output/" + wildcards.replicate + "/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=REPLICATE_LOOKUP[wildcards.replicate])
    return(inputs)

rule merge:
    input: mergeInput
    output: "output/{replicate}/bam/{treatment}.bam"
    message: "-----Merging bam files by experiment-----"
    log: "output/{replicate}/bam/{treatment}_merge.log"
    threads: config["threads"]["merge"]
    run:
        shell("sambamba merge --nthreads {threads} {output} {input}")
        # Create bam index folder for HTseq
        shell("samtools index {output}")

rule TPMCalculator:
    input: expand("output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)
    output: "output/counts/tpmcalculator/tpm-calculator.tpm"
    message: "------- Calculating TPM --------"
    log: "output/counts/tpmcalculator/tpm-calculator.log"
    threads: config["threads"]["TPMCalculator"]
    run:
        shell("TPMCalculator \
               --mode intersection-strict \
               --type exon \
               --idattr gene_id \
               -g " + config_dict["GTFname"] + " " +" ".join(expand("-b output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)))

rule featureCounts:
    message: "-----Generating raw counts (featureCounts)-----"
    input: expand("output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)
    output: 
        raw = "output/counts/featureCounts/featureCounts.cnt",
        cleaned = "output/counts/featureCounts/featureCount_clean.cnt"
    log: "output/counts/featureCounts/featureCounts.log"
    threads: config["threads"]["featureCount"]
    run:
        shell("featureCounts -o {output.raw}  -T {threads} -p -a " + config_dict["GTFname"] + " {input} \
    		2> {log}")
        shell("cat {input} |  egrep -v '#' | \
                sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g; s/output\/[A-Za-z0-9_-]*\/bam\///g' \
                > {output.cleaned}")

rule HTseq:
    message: "-----Generating raw counts (HTseq)-----"
    input:expand("output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)
    output: "output/counts/htseq/htseq-count.tsv"
    log: "output/counts/htseq/HTseq.log"
    threads: config["threads"]["HTseq"]
    run:
        # Compute raw gene-wise counts
        shell("htseq-count \
                 --format bam \
                 --order pos \
                 --mode intersection-strict \
                 --type exon \
                 --idattr gene_id \
                 --nprocesses {threads} \
                 --counts_output {output} \
                 {input} " + config_dict["GTFname"] + "  \
                 2> {log}")
        # Add headers to label HTseq tsv output fields
        shell("sed -i '1 i\gene\\t" + "\\t".join(input) + "' {output} &&\
               sed -i s/'output\/[A-Za-z0-9_-]*\/bam\/'//g {output}")

rule normalizeFeatureCounts:
    input:
        "output/counts/featureCounts/featureCount_clean.cnt"
    output:
        "output/counts/featureCounts/featureCounts.tpm.tsv",
        "output/counts/featureCounts/featureCounts.fpkm.tsv"
    message: "------- Normalizing featureCounts ---------"
    log: "output/counts/normalize_featureCounts.log"
    threads: config["threads"]["normalizeFeatureCounts"]
    run:
        shell("python normalizeCounts.py featureCounts " + config["GTFname"] + " {input} output/counts/featureCounts/featureCounts")

rule normalizeHTseq:
    input:
        "output/counts/htseq/htseq-count.tsv"
    output:
        tpm = "output/counts/htseq/htseq-count.tpm.tsv",
        fpkm = "output/counts/htseq/htseq-count.fpkm.tsv",
        cleaned = "output/counts/htseq/htseq-count_clean.cnt"
    message: "------- Normalizing HTseq ---------"
    log: "output/counts/htseq/normalize_HTseq.log"
    threads: config["threads"]["normalizeHTseq"]
    run:
        shell("sed 's/\.bam//g' {input} > {output.cleaned}")
        shell("python normalizeCounts.py HTseq " + config["GTFname"] + " {output.cleaned} output/counts/htseq/htseq-count")

rule DEG:
    input: "output/counts/featureCounts/featureCount_clean.cnt"
    output: "output/DEG/featureCount_clean.cnt." + CONTRASTS_FILE.iloc[1][0]  + "_vs_" + CONTRASTS_FILE.iloc[-1][0]  + ".DESeq2.DE_results"
    message: "-------Calculating DEGs {output}---------"
    log: "output/DEG/DEG.log"
    threads: config["threads"]["DEG"]
    run:
        shell("run_DE_analysis.pl --matrix {input} --method DESeq2 --samples_file contrasts.txt --output output/DEG")

rule summarize:
    input: "output/counts/tpmcalculator/tpm-calculator.tpm"
    output: "output/counts/tpmcalculator/tpm-calculator.tpm.mean.tsv"
    message: "--------------Summarizing Normalized Counts--------------"
    threads: config["threads"]["summarize"]
    run:
        shell("./summarizeNormalizedCounts.py {input}")
