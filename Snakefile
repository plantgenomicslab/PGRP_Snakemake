import sys
import os
import pandas as pd
SAMPLE_FILE = pd.read_table('sraid.txt', sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()

# Create sample output folders
os.makedirs("output/counts/", exist_ok=True)
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
        expand("output/{sample}/raw/{sample}_{replicate}.fastq.gz", sample=SAMPLE_LIST, replicate=["1", "2"]),
        expand("output/{sample}/trim/{sample}_{replicate}.fq.gz", sample=SAMPLE_LIST, replicate=["1", "2"]),
        expand("output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST),
        "output/counts/featureCounts.cnt"

rule fetchSRA:
    output: "output/sra/{sample}.sra"
    message: "-----Fetching {wildcards.sample} SRA files-----"
    threads: config["threads"]["fetchSRA"]
    run:
        shell("prefetch {wildcards.sample} --output-file output/sra/{wildcards.sample}.sra")

rule downloadSRA:
    input: "output/sra/{sample}.sra"
    output: 
        "output/{sample}/raw/{sample}_1.fastq.gz", 
        "output/{sample}/raw/{sample}_2.fastq.gz"
    message: "-----Downloading {wildcards.sample} Fastq files-----"
    threads: config["threads"]["downloadSRA"]
    log: "output/{sample}/logs/{sample}_downloadSRA.log"
    run:
        shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-3 --gzip \
				--outdir output/{wildcards.sample}/raw \
				2> {log}")

rule trim:
    input: 
        fwd_fastq = "output/{sample}/raw/{sample}_1.fastq.gz",
        rev_fastq = "output/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{sample}/trim/{sample}_1.fq.gz",
        "output/{sample}/trim/{sample}_2.fq.gz"
    message: "-----Trimming {wildcards.sample}-----"
    log: "output/{sample}/logs/{sample}_trim.log"
    threads: config["threads"]["trim"]
    run:
        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
             --cores 2 --max_n 40 --gzip -o output/{wildcards.sample}/trim {input.fwd_fastq} \
             {input.rev_fastq} 2> {log}")
        shell("for file in output/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc {output}")

rule align:
    input: 
        "output/{sample}/trim/{sample}_1.fq.gz",
        "output/{sample}/trim/{sample}_2.fq.gz"
    output:"output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
    message: "-----Aligning {wildcards.sample}-----"
    log: "output/{sample}/logs/{sample}_align.log"
    threads: config["threads"]["align"]
    run:
        shell("STAR --runMode alignReads --runThreadN 1 \
                --readFilesCommand gunzip -c --outFilterMultimapNmax 10 --alignIntronMin 25 \
                --alignIntronMax 10000 --genomeDir " + config["genomeDir"]  + " --readFilesIn \
                output/{wildcards.sample}/trim/{wildcards.sample}_1.fq.gz \
                output/{wildcards.sample}/trim/{wildcards.sample}_2.fq.gz \
			--outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.sample}/bam/{wildcards.sample}.bam 2> {log}")

rule featureCounts:
    message: "-----Generating feature counts-----"
    input: expand("output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST)
    output: "output/counts/featureCounts.cnt"
    log: "output/counts/featureCounts.log"
    shell: "featureCounts -o output/counts/featureCounts.cnt -p -a " + config_dict["GTFname"] + " {input} 2> {log}"
