import sys
import os
import pandas as pd
SAMPLE_FILE = pd.read_table('sraid.txt', sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()

# Create sample output folders
for path in SAMPLE_LIST:
    os.makedirs("output/" + path + "/stat/", exist_ok=True)
    os.makedirs("output/" + path + "/raw/", exist_ok=True)
    os.makedirs("output/" + path + "/SRA/", exist_ok=True)
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
        expand("output/{sample}/SRA/{sample}.sra", sample=SAMPLE_LIST),
        expand("output/{sample}/raw/{sample}_{replicate}.fastq.gz", sample=SAMPLE_LIST, replicate=["1", "2"]),
        expand("output/{sample}/trim/{sample}_{replicate}.fq.gz", sample=SAMPLE_LIST, replicate=["1", "2"]),
        expand("output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST)

rule downloadSRA:
    message: "-----Downloading SRA files-----"
    output: "output/{sample}/raw/{sample}_{replicate}.fastq.gz"
    threads: config["threads"]["downloadSRA"]
    log: "output/{sample}/logs/{sample}_downloadSRA.log"
    run:
        shell("prefetch {wildcards.sample} --output-file output/{wildcards.sample}/SRA")
        shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-3 --gzip \
				--outdir output/{wildcards.sample}/raw \
				2> {log}")
        #for sample in config["samples"]:
        #    shell("prefetch " + sample)
        #    shell("parallel-fastq-dump --sra-id " + sample +
        #          " --threads " + str(config["maxThreads"]) +
        #          " --split-3 --outdir ./raw_data  --gzip 2> {log}")

rule trim:
    input: "output/{sample}/raw/{sample}_{replicate}.fastq.gz"
    output:"output/{sample}/trim/{sample}_{replicate}.fq.gz"
    log: "output/{sample}/logs/{sample}_trim.log"
    threads: config["threads"]["trim"]
    run:
        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 \
             5 --cores 2 --max_n 40 --gzip -o trim raw_data/{wildcards.sample}_1.fastq.gz \
             raw_data/{wildcards.sample}_2.fastq.gz 2> {log}")
        shell("for file in output/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc " + " ".join(expand("output/{sample}/trim/{sample}_{replicate}.fq.gz",sample=SAMPLE_LIST,replicate=["1", "2"])))

        #for sample in config["samples"]:
        #        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2" +
        #              " 5 --cores 2 --max_n 40 --gzip -o trim raw_data/" + sample +
        #              "_1.fastq.gz raw_data/" + sample + "_2.fastq.gz 2> {log}")
        #shell("for file in trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        #shell("fastqc " + " ".join(expand("trim/{sample}_{replicate}.fq.gz",sample=config["samples"],replicate=config["replicates"])))

rule sample:
    input: "output/{sample}/trim/{sample}_{replicate}.fq.gz"
    output: "output/{sample}/trim/sample/{sample}_{replicate}_1k.fq.gz"
    run:
        shell("zcat < output/{wildcards.sample}/trim/{wildcards.sample}_1.fq.gz | seqkit \
                sample -n 1000 -o output/{wildcards.sample}/trim/sample/{wildcards.sample}_1_1k.fq.gz")
        shell("zcat < output/{wildcards.sample}/trim/{wildcards.sample}_2.fq.gz | seqkit \
                sample -n 1000 -o trim/sample/{wildcards.sample}_2_1k.fq.gz")
        #for sample in config["samples"]:
        #    shell("zcat < trim/" + sample + "_1.fq.gz | seqkit sample -n 1000 -o trim/sample/" + sample + "_1_1k.fq.gz")
        #    shell("zcat < trim/" + sample + "_2.fq.gz | seqkit sample -n 1000 -o trim/sample/" + sample + "_2_1k.fq.gz")

rule align:
    input: "output/{sample}/trim/sample/{sample}_{replicate}_1k.fq.gz"
    output:"output/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
    log: "output/{sample}/logs/{sample}_align.log"
    threads: config["threads"]["align"]
    run:
        shell("STAR --runMode alignReads --runThreadN 1 \
                --readFilesCommand gunzip -c --outFilterMultimapNmax 10 --alignIntronMin 25 \
                --alignIntronMax 10000 --genomeDir " + config["genomeDir"]  + "--readFilesIn \
                output/{wildcards.sample}/trim/sample/{wildcards.sample}_1_1k.fq.gz \
                output/{wildcards.sample}/trim/sample/{wildcards.sample}_2_1k.fq.gz \
                --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + config["bamPath"] + "{wildcards.sample}.bam")
        #for sample in config["samples"]:
        #    shell("STAR --runMode alignReads --runThreadN 1" +
        #          " --readFilesCommand gunzip -c --outFilterMultimapNmax 10 --alignIntronMin 25" +
        #          " --alignIntronMax 10000 --genomeDir " + config["genomeDir"] +
        #          " --readFilesIn trim/sample/" + sample + "_1_1k.fq.gz trim/sample/" + sample + "_2_1k.fq.gz" +
        #          " --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + config["bamPath"] + sample + ".bam")
