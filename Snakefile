import sys
import os
import pandas as pd
SAMPLE_FILE = pd.read_table('sraid.txt', sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()
REPLICATE_LIST=["1", "2"]

# Create sample output folders
os.makedirs("output/counts/", exist_ok=True)
os.makedirs("output/sra/", exist_ok=True)
os.makedirs("output/sra/logs/", exist_ok=True)
os.makedirs("output/logs/", exist_ok=True)

for path in SAMPLE_LIST:
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

rule all:
    input:
        expand("output/sra/{sample}.sra", sample=SAMPLE_LIST),
        expand("output/{sample}/raw/{sample}_{replicate}.fastq.gz", sample=SAMPLE_LIST, replicate=REPLICATE_LIST),
        expand("output/{sample}/raw/{sample}_{replicate}_fastqc.zip", sample=SAMPLE_LIST, replicate=REPLICATE_LIST),
        expand("output/{sample}/trim/{sample}_{replicate}_val_{replicate}.fq.gz", sample=SAMPLE_LIST, replicate=REPLICATE_LIST),
        expand("output/{sample}/trim/{sample}_{replicate}.fq.gz", sample=SAMPLE_LIST, replicate=REPLICATE_LIST),
        expand("output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST),
        "output/counts/featureCounts.cnt"

rule fetchSRA:
    output: "output/sra/{sample}.sra"
    message: "-----Fetching {wildcards.sample} SRA files-----"
    threads: config["threads"]["fetchSRA"]
    log: "output/sra/logs/{sample}_downloadSRA.log"
    run:
        shell("prefetch {wildcards.sample} --output-file output/sra/{wildcards.sample}.sra 2> {log}")

rule convertSRAtoFastq:
    input: "output/sra/{sample}.sra"
    output:
        "Output/{sample}/raw/{sample}_1.fastq.gz",
        "output/{sample}/raw/{sample}_2.fastq.gz"
    message: "-----Downloading {wildcards.sample} Fastq files-----"
    threads: config["threads"]["convertSRAtoFastq"]
    log: "output/{sample}/logs/{sample}_fastqdump.log"
    run:
        shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-3 --gzip \
				--outdir output/{wildcards.sample}/raw \
				2> {log}")

rule fastqc_raw:
    input:
        "output/{sample}/raw/{sample}_1.fastq.gz",
        "output/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{sample}/raw/{sample}_1_fastqc.zip",
        "output/{sample}/raw/{sample}_2_fastqc.zip"
    message: "-----Running Fastqc_raw {wildcards.sample}-----"
    threads: config["threads"]["fastqc_raw"]
    log: "output/{sample}/logs/{sample}_raw_fastqc.log"
    run:
        shell("fastqc --threads {threads} {output} 2> {log}")

rule trim:
    input:
        fwd_fastq = "output/{sample}/raw/{sample}_1.fastq.gz",
        rev_fastq = "output/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{sample}/trim/{sample}_1_val_1.fq.gz",
        "output/{sample}/trim/{sample}_2_val_2.fq.gz"
    message: "-----Trimming {wildcards.sample}-----"
    log: "output/{sample}/logs/{sample}_trim.log"
    threads: config["threads"]["trim"]
    run:
        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
             --cores 2 --max_n 40 --gzip -o output/{wildcards.sample}/trim {input.fwd_fastq} {input.rev_fastq} \
             2> {log}")

rule fastqc_trim:
    input:
        "output/{sample}/trim/{sample}_1_val_1.fq.gz",
        "output/{sample}/trim/{sample}_2_val_2.fq.gz"
    output:
        "output/{sample}/trim/{sample}_1.fq.gz",
        "output/{sample}/trim/{sample}_2.fq.gz"
    message: "-----Running Fastqc_trim {wildcards.sample}-----"
    threads: config["threads"]["fastqc_trim"]
    log: "output/{sample}/logs/{sample}_trim_fastqc.log"
    run:
        shell("for file in output/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc --threads {threads} {output} 2> {log}")

rule align:
    input:
        fwd_fastq = "output/{sample}/trim/{sample}_1.fq.gz",
        rev_fastq = "output/{sample}/trim/{sample}_2.fq.gz"
    output:"output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
    message: "-----Aligning {wildcards.sample}-----"
    log: "output/{sample}/logs/{sample}_align.log"
    threads: config["threads"]["align"]
    run:
        shell("STAR --runMode alignReads --runThreadN {threads} \
                --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 \
		--genomeDir " + config["genomeDir"]  + " \
	        --readFilesCommand gunzip -c --readFilesIn {input.fwd_fastq} {input.rev_fastq} \
                --outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.sample}/bam/{wildcards.sample}.bam  \
		2> {log}")

rule featureCounts:
    message: "-----Generating feature counts-----"
    input: expand("output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST)
    output: "output/counts/featureCounts.cnt"
    log: "output/counts/featureCounts.log"
    threads: config["threads"]["featureCount"]
    run:
        shell("featureCounts -o output/counts/featureCounts.cnt  -T {threads} -p -a " + config_dict["GTFname"] + " {input} \
    		2> {log}")
